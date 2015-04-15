'''
Get SV output, run post processing
'''

import os
import string
import numpy as np
import itertools
import ipdb
import pysam
import pandas as pd
import pandasql
import bamtools
import sys
import sqlite3
import subprocess

tr = 10 #threshold by how much read has to overlap breakpoint
sc_len = 25 #soft-clip threshold by which we call split reads
window = 500
max_cn=15

read_dtype = [('query_name', 'S25'), ('chrom', 'S10'), ('ref_start', int), ('ref_end', int), \
              ('align_start', int), ('align_end', int), ('len', int), ('is_reverse', np.bool)]

def read_to_array(x,bamf):
    chrom = bamf.getrname(x.reference_id)
    try:
        read = np.array((x.query_name,chrom,x.reference_start,x.reference_end,x.query_alignment_start,
                         x.query_alignment_end,x.query_length,np.bool(x.is_reverse)),dtype=read_dtype)
        return read
    except TypeError:
        print 'Warning: record %s contains invalid attributes' % x.query_name
        return np.empty(len(read_dtype),dtype=read_dtype)

def is_soft_clipped(r):
    return r['align_start'] != 0 or (r['len'] + r['ref_start'] != r['ref_end'])

def is_normal_across_break(r,pos):
   return r['ref_start'] < (pos - tr) and r['ref_end'] > (pos + tr) and \
       r['align_start'] < (tr*2) and (r['align_end'] + r['ref_start'] - r['ref_end'] < (tr*2))

def is_normal_spanning(r,m,pos):
    if not (is_soft_clipped(r) or is_soft_clipped(m)):
        if (not r['is_reverse'] and m['is_reverse']) or (r['is_reverse'] and not m['is_reverse']):
            return r['ref_end'] < (pos + tr) and m['ref_start'] > (pos - tr)
    return False

def is_supporting_split_read(r,pos):
    """
    return whether read is a supporting split read
    doesn't yet check whether the soft-clip aligns
    to the other side - should it, or is that implicit?
    """
    #TODO: check that split read has valid insert size
    if r['align_start'] < (tr/2): #a "soft" threshold if it is soft-clipped at the other end
        return r['ref_end'] > (pos - tr) and r['ref_end'] < (pos + tr) and \
            (r['len'] - r['align_end'] >= sc_len)
    else:
        return r['ref_start'] > (pos - tr) and r['ref_start'] < (pos + tr) and \
            (r['align_start'] >= sc_len)

def get_sc_bases(r,pos):
    """
    Return the number of soft-clipped bases
    """
    if r['align_start'] < (tr/2):
        return r['len'] - r['align_end']
    else:
        return r['align_start']
    

def get_bp_dist(x,bp_pos):
    if x['is_reverse']: 
        return (x['ref_end'] - bp_pos)
    else: 
        return (bp_pos - x['ref_start'])

def is_supporting_spanning_pair(r,m,bp1,bp2,inserts):
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    dir1 = bp1['dir']
    dir2 = bp2['dir']

    if is_soft_clipped(r) or is_soft_clipped(m):
        return False

#    #check if pair is facing the right way (only if it is a duplication)
#    if sv_type == 'DUP':
#        #if dup is the "correct way around" i.e. pos1 < pos2
#        if dir1=='-' and dir2=='+':
#            if not (not r['is_reverse'] and m['is_reverse']):
#                return False
#        #dup is wrong way round
#        elif dir1=='+' and dir2=='-': 
#            if not (r['is_reverse'] and not m['is_reverse']):
#                return False
#        elif dir1 == dir2:
#            return False #not a duplicaton?
    
    max_ins = 2*inserts[1]+inserts[0]
    #ensure this isn't just a regular old spanning pair    
    if r['chrom']==m['chrom']:
        if r['ref_start']<m['ref_start']:
            if m['ref_start']-r['ref_end'] < max_ins: return False
        else:
            if r['ref_start']-m['ref_end'] < max_ins: return False

    ins_dist1 = get_bp_dist(r,pos1)
    ins_dist2 = get_bp_dist(m,pos2)

    if ins_dist1<0 or ins_dist2<0:
        return False
    else:
        return (ins_dist1+ins_dist2) < max_ins

def get_loc_reads(bp,bamf):
    loc = '%s:%d:%d' % (bp['chrom'], bp['start'], bp['end'])
    
    loc_reads = np.empty([0,len(read_dtype)],dtype=read_dtype)    
    iter_loc = bamf.fetch(region=loc)
    for x in iter_loc:
        read = read_to_array(x,bamf) 
        loc_reads = np.append(loc_reads,read)
    
    loc_reads = np.sort(loc_reads,axis=0,order=['query_name','ref_start'])
    return loc_reads

def reads_to_sam(reads,bam,bp1,bp2,name):
    """
    For testing read assignemnts.
    Takes reads from array, matches them to bam 
    file reads by query name and outputs them to Sam
    """
    bamf = pysam.AlignmentFile(bam, "rb")
    loc1 = '%s:%d:%d' % (bp1['chrom'], bp1['start'], bp1['end'])
    loc2 = '%s:%d:%d' % (bp2['chrom'], bp2['start'], bp2['end'])
    iter_loc1 = bamf.fetch(region=loc1)
    iter_loc2 = bamf.fetch(region=loc2)
    
    loc1 = '%s-%d' % (bp1['chrom'], (bp1['start']+bp1['end'])/2)
    loc2 = '%s-%d' % (bp2['chrom'], (bp1['start']+bp1['end'])/2)
    sam_name = '%s_%s-%s' % (name,loc1,loc2)
    bam_out = pysam.AlignmentFile('%s.sam'%sam_name, "w", header=bamf.header)
    
    for x in iter_loc1:
        [bam_out.write(x) for r in reads if r['query_name']==x.query_name]
    for x in iter_loc2:
        [bam_out.write(x) for r in reads if r['query_name']==x.query_name]
    
    bamf.close()
    bam_out.close()
    
def windowed_norm_read_count(loc_reads,inserts):
    max_ins = 3*inserts[1]+inserts[0]
    cnorm = 0
    for idx,r in enumerate(loc_reads):
        if idx+1 >= len(loc_reads):
            break    
        r1 = np.array(loc_reads[idx],copy=True)
        r2 = np.array(loc_reads[idx+1],copy=True)
        if r1['query_name']!=r2['query_name'] or r1['chrom']!=r2['chrom']:
            continue
        ins_dist = r2['ref_end']-r1['ref_start']
        facing = not r1['is_reverse'] and r2['is_reverse']
        if not is_soft_clipped(r1) and not is_soft_clipped(r2) and facing and ins_dist > 0 and ins_dist < max_ins:
            cnorm = cnorm + 2
    return cnorm

def get_loc_counts(loc_reads,pos,rc,reproc,split,bp_num=1):
    for idx,x in enumerate(loc_reads):
        if idx+1 >= len(loc_reads):            
            break        
        r1 = loc_reads[idx]
        r2 = loc_reads[idx+1] if (idx+2)<=len(loc_reads) else None
        if is_normal_across_break(x,pos):
            split_norm = 'bp%d_split_norm'%bp_num
            rc[split_norm] = rc[split_norm]+1 
        elif is_supporting_split_read(x,pos):
            split = np.append(split,x)            
            split_supp = 'bp%d_split'%bp_num
            split_cnt = 'bp%d_sc_bases'%bp_num
            rc[split_supp] = rc[split_supp]+1 
            rc[split_cnt] = rc[split_cnt]+get_sc_bases(x,pos)
        elif r2!=None and r1['query_name']==r2['query_name'] and is_normal_spanning(r1,r2,pos):
            span_norm = 'bp%d_span_norm'%bp_num
            rc[span_norm] = rc[span_norm]+1 
        else:
            reproc = np.append(reproc,x) #may be spanning support or anomalous
    return rc, reproc, split

def get_sv_read_counts(bp1,bp2,bam,columns,inserts,max_dp):
    bamf = pysam.AlignmentFile(bam, "rb")
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    loc1_reads = get_loc_reads(bp1,bamf)    
    loc2_reads = get_loc_reads(bp2,bamf)    
    
    print('%d reads extracted'%len(loc1_reads))
    print('%d reads extracted'%len(loc2_reads))
    if len(loc1_reads)>max_dp or len(loc2_reads)>max_dp:
        return pd.Series
    bamf.close() 
    
    rc = pd.Series(np.zeros(len(columns)),index=columns,dtype='int')
    #rc['sv']='%s:%d-%s:%d'%(bp1['chrom'],pos1,bp2['chrom'],pos2)
    reproc = np.empty([0,len(read_dtype)],dtype=read_dtype)
    
    split = np.empty([0,len(read_dtype)],dtype=read_dtype)
    rc, reproc, split = get_loc_counts(loc1_reads,pos1,rc,reproc,split)
    rc, reproc, split = get_loc_counts(loc2_reads,pos2,rc,reproc,split,2)
    rc['bp1_win_norm'] = windowed_norm_read_count(loc1_reads,inserts)
    rc['bp2_win_norm'] = windowed_norm_read_count(loc2_reads,inserts)
     
    reproc = np.sort(reproc,axis=0,order=['query_name','ref_start'])
    span = np.empty([0,len(read_dtype)],dtype=read_dtype)
    for idx,x in enumerate(reproc):
        if idx+1 >= len(reproc):
            break        
        if reproc[idx+1]['query_name'] != reproc[idx]['query_name']:
            #not paired
            #rc['bp1_anomalous'] = rc['bp1_anomalous']+1 
            continue
        mate = np.array(reproc[idx+1],copy=True)
        r1 = np.array(x,copy=True)
        r2 = np.array(mate,copy=True)
        #if read corresponds to bp2 and mate to bp1
        if (bp1['chrom']!=bp2['chrom'] and x['chrom']==bp2['chrom']) or \
            (pos1 > pos2 and bp1['chrom']==bp2['chrom']):
            r1 = mate
            r2 = np.array(x,copy=True)
        if is_supporting_spanning_pair(r1,r2,bp1,bp2,inserts):
            span = np.append(span,r1)            
            rc['spanning'] = rc['spanning']+1 
        #else:
        #    rc['bp1_anomalous'] = rc['bp1_anomalous']+1
    #reads_to_sam(span,bam,bp1,bp2,'span')
    #reads_to_sam(split,bam,bp1,bp2,'split')
    return rc

def proc_header(header,columns):
    hd = pd.read_csv(header,delimiter='=',header=None)
    hdt = pd.DataFrame(hd[1].values,index=hd[0].values).transpose()
    keyfields = ['bp1_chr','bp1_pos','bp1_dir','bp2_chr','bp2_pos','bp2_dir','classification']
    if not np.all(np.sort(map(str,columns))==np.sort(hd[1].values)):
        print('Headers in cfg and input file do not match! Exiting')
        sys.exit()
    try:
        return hdt[keyfields].values[0]
    except KeyError: 
        print('Key fields incorrect. Set key fields as bp1_chr, bp1_pos, bp1_dir, bp2_chr,bp2_pos, bp2_dir and classification')
        sys.exit()

def run(svin,bam,out,header,db_out,mean_dp): 
    inserts = bamtools.estimateInsertSizeDistribution(bam)
    rlen = bamtools.estimateTagSize(bam)
    max_dp = ((mean_dp*(window*2))/rlen)*max_cn

    bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S2')]
    sv_dtype = [('bp1_chr','S20'),('bp1_pos',int),('bp1_dir','S5'),('bp2_chr','S20'), \
                ('bp2_pos',int),('bp2_dir','S5'),('classification','S100')]               
    #TODO: make this flexible
    svs = pd.read_csv(svin,delimiter='\t',dtype=sv_dtype)

    columns = ['bp1_split_norm','bp1_span_norm','bp1_win_norm','bp1_split','bp1_sc_bases', \
               'bp2_split_norm','bp2_span_norm','bp2_win_norm','bp2_split','bp2_sc_bases','spanning']
    #rinfo = pd.DataFrame(columns=columns,index=range(len(svs.index)),dtype='float')
   
    bp1_chr,bp1_pos,bp1_dir,bp2_chr,bp2_pos,bp2_dir,classification = proc_header(header,svs.columns)

    for idx,row in svs.iterrows():
        sv_str = '%s:%d|%s:%d'%(row[bp1_chr],row[bp1_pos],row[bp2_chr],row[bp2_pos])
        bp1 = np.array((row[bp1_chr],row[bp1_pos]-window,row[bp1_pos]+window,row[bp1_dir]),dtype=bp_dtype)
        bp2 = np.array((row[bp2_chr],row[bp2_pos]-window,row[bp2_pos]+window,row[bp2_dir]),dtype=bp_dtype)
        sv_rc  = get_sv_read_counts(bp1,bp2,bam,columns,inserts,max_dp)
        if bool(sv_rc.empty):
            print('Skipping %s, read depth too high!'%sv_str)
            continue
        newrow = pd.DataFrame(row.append(sv_rc)).transpose()
        newrow['norm1']=sv_rc['bp1_split_norm']+sv_rc['bp1_span_norm']
        newrow['norm2']=sv_rc['bp2_split_norm']+sv_rc['bp2_span_norm']
        newrow['support']=sv_rc['bp1_split']+sv_rc['bp2_split']+sv_rc['spanning']
        newrow['vaf1']=newrow['support']/(newrow['support']+newrow['norm1'])
        newrow['vaf2']=newrow['support']/(newrow['support']+newrow['norm2'])
        #rinfo.loc[idx] = sv_rc
        con = sqlite3.connect(db_out)
        newrow.to_sql('sv_info',con,if_exists='append',index=False)
        con.close()
        print('processed %s'%sv_str)
    con = sqlite3.connect(db_out)
    pd.read_sql('select * from sv_info',con).to_csv(out,sep="\t",index=False)
    con.close()
    #con.close()
    #subprocess.call(['samtools','view','-H',bam],stdout=open('head.sam','w'))
    #subprocess.call(['cat','head.sam','span_*.sam'],stdout=open('spanning_all.sam','w'),shell=True)
    #subprocess.call(['samtools','view','-hbS','spanning_all.sam'],stdout=open('spanning_all.bam','w'))
    #subprocess.call(['samtools','sort','spanning_all.bam','spanning_all_sort'])
    #subprocess.call(['samtools','index','spanning_all_sort.bam'])

#    rinfo = svs.join(rinfo)
#    #calculate VAFs
#    vafs = pd.DataFrame(columns=['bp1_norm','bp2_norm','support','vaf1','vaf2'],index=range(len(svs.index)))
#    for idx,ri in rinfo.iterrows():
#        support = (float(ri['bp1_split'])+float(ri['bp2_split'])+float(ri['spanning']))
#        norm1 = float(ri['bp1_split_norm'])+float(ri['bp1_span_norm'])
#        norm2 = float(ri['bp2_split_norm'])+float(ri['bp2_span_norm'])
#        vafs['bp1_norm'][idx] = norm1
#        vafs['bp2_norm'][idx] = norm2
#        vafs['support'][idx] = support
#        vafs['vaf1'][idx] = support / (support + norm1)
#        vafs['vaf2'][idx] = support / (support + norm2)
#    rinfo = rinfo.join(vafs)
#    print rinfo
#    rinfo.to_csv(out,sep="\t",index=False)    

