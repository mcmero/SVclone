'''
Takes an SV matrix and attempts to build a phylogeny
'''
import numpy as np
import os
import string
import subprocess
import random
import ipdb
import csv
import random
import graph as gp
import phylo_sv as sv
import itertools
from itertools import combinations
from subprocess import call
from numpy import loadtxt
from numpy import zeros
from scipy import stats

def get_allocs(class_mat,svmat):
    allocs = []
    for sv in svmat:
        match_locs = np.all(sv==class_mat,axis=1)
        x = np.where(match_locs)[0]
        allocs = np.append(allocs,x)
    return allocs

def unique2D(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]

def get_mp(nshared,svrot=True):
    '''
    Construct M' matrix
    Obtain binary encoding for each feature, sort by binary value, 
    then rotate and reverse the resulting matrix
    @param nshared nxm matrix of samples (cols) and rows (SVs)
    @param svrot is true, the matrix is of form rows = svs, cols = samples
    '''
    if not svrot:
        nshared = np.rot90(nshared)
    nshared = unique2D(nshared)
    ncol = len(nshared[0])    
    binrs = []
    for ns in nshared:
        binr = '0b'+''.join(map(str,ns[:ncol]))
        binrs.append(int(binr,2))
    # remove any zeroes
    if 0 in binrs:
        nshared = np.delete(nshared,binrs.index(0),axis=0)
        binrs.remove(0)
    ods = np.argsort(binrs)[::-1]
    nshared = nshared[ods] #reorder the array
    mp = np.rot90(nshared[:,:ncol])[::-1] #rotate/reverse
    return mp

def get_k1(mp,nodes):
    '''
    @param mp the m prime matrix
    @param nodes corresponding to the matrix
    Generate k1 from m' matrix
    Allows for checking of perfect phylogeny
    '''
    k1 = np.empty([0,len(mp[0])],dtype='|S15')
    for m in mp:
        chars = nodes[m!=0]#characters in row
        mrow = np.zeros(len(mp[0]),dtype='|S15')
        for idx,c in enumerate(chars):
            mrow[idx] = c
        if len(chars)<len(mp[0]): mrow[len(chars)]='%'
        k1 = np.append(k1,[mrow],axis=0)
    return k1

def tree_exists(k1,nodes):
    '''
    Determine whether perfect phylogeny exists from a k1 matrix    
    '''
    nlocs = []
    for n in nodes:
        npos = set([])
        for k in k1:
            [npos.add(loc) for loc in list(np.where(k==n)[0])]
        nlocs.append(npos)    
    nlocs = np.array([len(locs)>1 for locs in nlocs])
    if np.any(nlocs):
        #print 'No phylogeny found!'
        return False
    
    print 'Success! Found phylogeny!'
    return True

def get_cons(x,cons):
    '''
    get all pairs in cons that contain element x
    '''
    ncons = []
    for cn in cons:
        if x in cn: 
            ncons.append(cn)
    return ncons

def get_k(k,q,con):
    '''
    return the first K vector with E[K] >= 1
    (the K vector contains all elements of a chain of connected 
    vertices with more than just one connection between them)
    '''
    if not q:
        return k
    else:
        q1 = q.pop()
        k.add(q1) 
        cn = get_cons(q1,con)
        cn = set([c for x in cn for c in x])
        for c in cn:
            if c not in k:
                q.add(c)
        return get_k(k,q,con)

def build_graph(m,s,c):
    '''
    take the _M_ matrix and return the corresponding connection graph    
    '''
    g = gp.Graph()
    for si in s:
        g.addVertex(si)
    for ci in c:
        g.addVertex(ci)

    for i in range(len(m)):
        nc = c[np.where(m[i]==1)]
        for ci in nc:
            g.addEdge(s[i],ci)
    return g

def get_edge_pairs(g):
    '''
    get all the pairwise connections of a graph
    '''
    con = []
    for v in g:
        for w in v.getConnections():
            con.append([v.getId(),w.getId()])
    return con

def remove_edges_with(k,con):
    '''
    removes any edge pairs containing any of the elements in k
    '''
    n_con = []
    for i in range(len(con)):
        if not any([ki in con[i] for ki in k]): n_con.append(con[i])
    return n_con

def solve_incomplete(m,s,c):
    '''
    attempt to solve phylogeny assuming row j is 'incomplete'
    method from Pe'er 2004
    @param m3 input matrix _M_
    @param j row to consider 'incomplete'
    @param s samples
    @param c character features (SVs)
    '''
   
    #remove any entries with no 0s
    m3 = m.copy()
    m_inf = m.copy()
    ct = c.copy()

    mi = np.empty(0,dtype='int')
    for i in range(len(m3[0])):
        mcol = m3[:,i:i+1]
        if not np.any(mcol==0):
            mi = np.append(mi,i) #col indicies to delete 
    m3 = np.delete(m3,mi,axis=1)
    ct = np.delete(ct,mi)
    
    if len(m3[0])==0:
        print 'Matrix M has no 0s!'
        return []
    
    #print 'm input'
    #print m3
    
    t = [] #initialise tree t
    g = build_graph(m3,s,ct)
    con = get_edge_pairs(g)
    q = set(con[0]) #pick the first elements as k
    k = get_k(set(),q,con[1:])

    u = []
    while len(con) > 1:
        # get a new k if E[K] <= 1
        while len(k) < 3:
            con = remove_edges_with(k,con)
            if len(con) > 1:
                con = con[1:]
                k = get_k(set(),set(con[0]),con)
            else: 
                break
        
        sp = set(s).intersection(k)
        cp = set(c).intersection(k)
       
        si = [s.index(sx) for sx in sp]
        m_tmp = m3[si]
        
        #if all values in col c for s in sp 
        #contain no 0s, add ci to u
        u = []
        for ci in cp:
            ctd = np.where(ci==ct)[0]
            mcol = m_tmp[:,ctd:ctd+1]
            if np.all(mcol!=0):
                cm = np.where(ci==c)[0]
                u.append(ci)
        if not u:
            break
      
        # set any incompletes to 1 if associated with sp
        for ci in cp:
            ctd = np.where(ci==ct)[0]
            mcol = m_tmp[:,ctd:ctd+1]
            if np.any(mcol==-1):
                mcol[mcol==-1] = 1
                cm = np.where(ci==c)[0]
                m_inf[si,cm:cm+1] = mcol

        con = remove_edges_with(u,con)
        t.append(sp)
        if con:
            k = get_k(set(),set(con[0]),con) #get a new k

    # make any -1s in the matrix 0s
    for i in range(len(m_inf[0])):
        tcol = m_inf[:,i:i+1]
        tcol[tcol==-1] = 0
        m_inf[:,i:i+1] = tcol
  
    #print 'm post inferral' 
    #print m_inf
    
    if not u:
        #print 'M has no perfect phylogeny!'
        return [m_inf,[]]
    
    print 'tree T'
    print t
    return [m_inf,t]

def get_ns(nshared,row):
    '''
    return number shared corresponding to the row given
    '''
    ncol = len(nshared[0])
    mrow = np.all(nshared[:,:ncol-1]==np.rot90(row)[0],axis=1)
    return(nshared[np.where(mrow)[0],ncol-1])

def get_nshared_rnk(nshared):
    '''
    Gets the nshared matrix ranked by num shared descending
    '''
    ncol = len(nshared[0])
    #mrow = np.all(nshared[:,:ncol-1]==np.rot90(row))
    ns = nshared[:,ncol-1:ncol]
    ns = [x for n in ns for x in n]   
    return(nshared[np.argsort(ns)][::-1])

def get_coords(mp,val=0):
    '''
    get coordinates of val in matrix
    by default, gives cell locations of all zeroes
    @mp M' matrix
    @val integer value to match
    '''
    zs = list()
    for i in range(len(mp)):
        zs.extend([(i,j) if mp[i,j]==val else (-1,-1) for j in range(len(mp[0]))])
        zs = filter(lambda x: sum(x)>=0,zs)
    return zs

def get_conflicts(mp):
    '''
    return conflicting columns in matrix
    @mp M' matrix
    '''
    #print 'testing matrix for conflicts\n'
    #print mp
    nc = len(mp[0])
    cg = list(combinations(range(len(mp[0])),2)) #combinations to test for conflict
    conf_list = list()
    for c in cg:
        comp = list()
        for i in range(len(mp)):
            comp.append((mp[i,c[0]],mp[i,c[1]]))
        if (1,0) in comp and (1,1) in comp and (0,1) in comp:
            #print 'cols %s,%s in conflict' % (c[0]+1,c[1]+1)
            conf_list.append((c[1]+1,c[0]+1)) #x,y 1-indexed
    return conf_list

def output_dot(pname,samples,nodes,cmat,ftab):
    '''
    Output drawn phylogenetic tree using dot
    '''
    tree = []    
    nums = ftab[1:,1]
    nums = [0 if m=='' else int(m) for m in ftab[:,1]]
    with open(tree_file,'r') as tf:
        for line in tf:
            tree.append(line.rstrip())
    with open('results/%s_trees.dot'%pname,'w') as df:
        df.write('digraph %s {\n'%pname)
        df.write('\tlabel="%s"\n'%("%s matched=%d unmatched=%s"%(pname,sum(nums[1:]),nums[0])))
        for t in tree:
            df.write('\t%s;\n'%t)
        df.write('\tgraph[size="7.75,10.25"]\n}\n')
        for idxs,s in enumerate(samples):
            df.write('digraph %s {\n\tlabel="%s";\n'%(s,s))
            for idxn,n in enumerate(nodes):
                if cmat[idxn,idxs]==1:
                    df.write('\t%s_%s[label="%s"];\n'%(s,n,str(nums[idxn+1])))
                else:
                    df.write('\t%s_%s[label="0"];\n'%(s,n))
            for t in tree:
                t = t.split('->')
                df.write('\t%s_%s -> %s_%s;\n' % (s,t[0].rstrip(),s,t[1].lstrip()))
            df.write('}\n') 
   
    subprocess.call(["dot","-Tps","results/%s_trees.dot"%pname],stdout=open(r'results/%s_trees.ps'%pname,'w')) 
