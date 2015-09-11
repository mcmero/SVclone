import numpy as np
import pandas as pd
import pymc as pm
import ipdb

from operator import methodcaller
from . import parameters as param

def get_weighted_cns(gtypes):
    gtypes_split = [map(methodcaller("split",","),x) for x in map(methodcaller("split","|"),gtypes)]    
    cn_vals = []
    for gtype in gtypes_split:
        cn_val = sum([int(eval(g[0]))+int(eval(g[1]))*float(eval(g[2])) if g!=[''] else 2 for g in gtype])
        cn_vals.append(cn_val)
    return np.array(cn_vals)/2

def normalise_wins_by_cn(df_flt):
    bp1_win = df_flt.bp1_win_norm.map(float).values
    bp2_win = df_flt.bp2_win_norm.map(float).values
    
    bp1_wcn = get_weighted_cns(df_flt.gtype1.values)
    bp2_wcn = get_weighted_cns(df_flt.gtype2.values)

    bp1_nonzero = np.logical_not(bp1_wcn==0)
    bp2_nonzero = np.logical_not(bp2_wcn==0)
    bp1_win[bp1_nonzero] = (bp1_win[bp1_nonzero]/bp1_wcn[bp1_nonzero])
    bp2_win[bp2_nonzero] = (bp2_win[bp2_nonzero]/bp2_wcn[bp2_nonzero])
    
    return bp1_win,bp2_win

def get_cn_mu_v(cn):
    cn_v = [0.,0.]
    mu_v = [0.,0.]

    c = cn.split(',')
    if len(c)<2:
        return tuple(cn_v),tuple(mu_v)

    c = map(float,c)
    cn_t = float(c[0]+c[1])
    cn_v[0] = float(cn_t)
    mu_v[0] = c[0]/cn_t if cn_t!=0 else 0.

    if c[0]!=c[1] and c[1]>0:
        cn_v[1] = cn_t
        mu_v[1] = c[1]/cn_t if cn_t!=0 else 0.

    return tuple(cn_v),tuple(mu_v)

#def get_allele_combos_tuple(c):
#    cn_r = [tuple([0.,0.]),tuple([0.,0.])]
#    cn_v = [tuple([0.,0.]),tuple([0.,0.])]
#    mu_v = [tuple([0.,0.]),tuple([0.,0.])]
#
#    if len(c) == 0 or c[0]=='':
#        return cn_r,cn_v,mu_v
#
#    cn_v[0],mu_v[0] = get_cn_mu_v(c[0])
#    cn1_tmp = map(float,c[0].split(',')) if len(c[0])>1 else c[0]
#    cnr1_tmp = cn1_tmp[0]+cn1_tmp[1]
#
#    if len(c) > 1:
#        cn_v[1],mu_v[1] = get_cn_mu_v(c[1])
#        cn2_tmp = map(float,c[1].split(','))
#        cnr2_tmp = cn2_tmp[0]+cn2_tmp[1]
#        cn_r[1],cn_r[0] = tuple([cnr1_tmp,cnr1_tmp]),tuple([cnr2_tmp,cnr2_tmp])
#    else:
#        cn_r[0],cn_r[1] = tuple([cnr1_tmp,cnr1_tmp]),tuple([cnr1_tmp,cnr1_tmp])
#
#    return tuple(cn_r),tuple(cn_v),tuple(mu_v)

#def get_sv_allele_combos_tuple(sv):
#    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
#    cnr_bp1,cnv_bp1,mu_bp1 = get_allele_combos(cn_tmp[0])
#    cnr_bp2,cnv_bp2,mu_bp2 = get_allele_combos(cn_tmp[1])
#
#    return pd.Series([tuple([cnr_bp1,cnr_bp2]),tuple([cnv_bp1,cnv_bp2]),tuple([mu_bp1,mu_bp2])])

def get_allele_combos(c):
    combos = []

    if len(c) == 0 or c[0]=='':
        return combos
    
    cn1_v,mu1_v = get_cn_mu_v(c[0])
    cn1 = map(float,c[0].split(',')) if len(c[0])>1 else c[0]
    cn1_ref_total  = sum(cn1[:2])
    cn1_r = tuple([cn1_ref_total, cn1_ref_total])
    combos.extend(zip(cn1_r,cn1_v,mu1_v))

    if len(c) > 1:
        cn2_v,mu2_v = get_cn_mu_v(c[1])
        cn2_tmp = map(float,c[1].split(','))
        cn2_ref_total = cn2_tmp[0]+cn2_tmp[1]
        cn2_r = tuple([cn2_ref_total, cn2_ref_total])
        combos[0],combos[1] = zip(cn2_r,cn1_v,mu1_v)
        combos.extend(zip(cn1_r,cn2_v,mu2_v))
    
    return filter_cns(combos)

def get_sv_allele_combos(sv):
    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
    combos_bp1 = get_allele_combos(cn_tmp[0])
    combos_bp2 = get_allele_combos(cn_tmp[1])

    return tuple([combos_bp1,combos_bp2])

def get_sv_vals(df,rlen):
    n = zip(np.array(df.norm1.values),np.array(df.norm2.values))
    s = np.array(df.bp1_split.values+df.bp2_split.values)
    d = np.array(df.spanning.values)
    #cn_r,cn_v,mu_v = [],[],[]
    
    combos = df.apply(get_sv_allele_combos,axis=1)
    #cn_r,cn_v,mu_v = combos[0],combos[1],combos[2]

#    for idx,sv in df.iterrows():
#        cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
#        cnr_bp1,cnv_bp1,mu_bp1 = get_allele_combos_tuple(cn_tmp[0])
#        cnr_bp2,cnv_bp2,mu_bp2 = get_allele_combos_tuple(cn_tmp[1])
#        cn_r.append(tuple([cnr_bp1,cnr_bp2]))
#        cn_v.append(tuple([cnv_bp1,cnv_bp2]))
#        mu_v.append(tuple([mu_bp1,mu_bp2]))

    sup = np.array(d+s,dtype=float)
    n_bp1,n_bp2 = [ni[0] for ni in n],[ni[1] for ni in n]
    Nvar = len(sup)
    
    # calculate the average depth using windowed counts
    win1 = (df.bp1_win_norm.values*rlen)/(param.window*2)
    win2 = (df.bp2_win_norm.values*rlen)/(param.window*2)
    av_cov = np.mean([win1,win2])
     
    sides = np.zeros(Nvar,dtype=int)
    sides[df.gtype1.values=='']=1 #no genotype for bp1, use bp2

    has_both_gts = np.logical_and(df.gtype1.values!='',df.gtype2.values!='')
    
    bp1_win, bp2_win = normalise_wins_by_cn(df)
    bp1_win, bp2_win = (bp2_win*rlen)/(param.window*2), (bp2_win*rlen)/(param.window*2)

    # for sides with gtypes for both sides, pick the side where the adjusted window count is closest to the average coverage
    av_cov = np.mean([bp1_win,bp2_win])
    dev_from_cov = np.array(zip(abs(bp1_win-av_cov),abs(bp2_win-av_cov)))

    # ALT: pick the side where the norm count is closest to the mean coverage
    #dev_from_cov = np.array(zip(abs(n_bp1-av_cov),abs(n_bp2-av_cov)))
    # OR: prefer sides with subclonal genotype data
    #gt1_sc = np.array(map(len,map(methodcaller("split","|"),df.gtype1.values)))>1
    #gt2_sc = np.array(map(len,map(methodcaller("split","|"),df.gtype1.values)))>1

    sides[has_both_gts] = [np.where(x==min(x))[0][0] for x in dev_from_cov[has_both_gts]]

    norm = [ni[si] for ni,si in zip(n,sides)]
    #norm = map(np.mean,n)
    dep = np.array(norm+sup,dtype=float)

    #return sup,dep,cn_r,cn_v,mu_v,sides,Nvar
    cn_states = [cn[side] for cn,side in zip(combos,sides)]
    return sup,dep,cn_states,Nvar

def get_snv_vals(df):
    n = df['ref'].values
    b = df['var'].values

    #cn_r,cn_v,mu_v = [],[],[]
    df = df.fillna('')
    df = df[df.chrom.values!='']

    def get_snv_allele_combos(snv):
        return get_allele_combos(snv.gtype.split('|'))
     
    combos = df.apply(get_snv_allele_combos,axis=1)
#    for idx,snv in df.iterrows():
#        cn_tmp = snv.gtype.split('|')
#        cnr,cnv,mu = get_allele_combos_tuple(cn_tmp)
#        cn_r.append([cnr])
#        cn_v.append([cnv])
#        mu_v.append([mu])
    return b,(n+b),combos,len(b)

def fit_and_sample(model, iters, burn, thin, use_map):
    #TODO: suppress warning about using fmin method
    if use_map:
        map_ = pm.MAP( model )
        map_.fit( method = 'fmin_powell' )
    mcmc = pm.MCMC( model )
    mcmc.sample( iters, burn=burn, thin=thin )
    return mcmc

def get_pv(phi,cn_r,cn_v,mu,pi):
    rem = 1.0 - phi
    rem = rem.clip(min=0)
    pn =  (1.0 - pi) * 2            #proportion of normal reads coming from normal cells
    pr =  pi * rem * cn_r           #proportion of normal reads coming from other clusters
    pv =  pi * phi * cn_v           #proportion of variant + normal reads coming from this cluster

    norm_const = pn + pr + pv
    pv = pv / norm_const

    return pv * mu

#def get_most_likely_cn(cn_r,cn_v,mu_v,si,di,phi,pi):
#    llik = []
#    vals = []
#    # for each subclone
#    for cn_rx, cn_vx, mu_vx in zip(cn_r,cn_v,mu_v):
#        # for each allelic state
#        for cn_ri, cn_vi, mu_vi in zip(cn_rx,cn_vx,mu_vx):
#            vals.append([cn_ri,cn_vi,mu_vi])
#            if cn_vi==0 or mu_vi==0:
#                llik.append(float('nan'))
#            else:
#                pv = get_pv(phi,cn_ri,cn_vi,mu_vi,pi)
#                temp = pm.binomial_like(si,di,pv)
#                llik.append(temp)
#    
#    vals = np.array(vals)
#    idx = np.where(np.array(llik)==np.nanmax(llik))[0]
#    if len(idx)==0:
#        #ipdb.set_trace()
#        # this shouldn't happen
#        return [0.,0.,0.]
#    else:
#        return vals[idx[0]]

def filter_cns(cn_states):
    cn_str = [','.join(map(str,cn)) for cn in cn_states if cn[2]!=0 and cn[1]!=0]
    cn_str = np.unique(np.array(cn_str))
    return [map(float,cn) for cn in map(methodcaller('split',','),cn_str)]

#def get_cn_states(cn_r,cn_v,mu_v,sides):
#    # get all unique allele state combos per SV
#    cnr = [cn[side] for cn,side in zip(cn_r,sides)]
#    cnv = [cn[side] for cn,side in zip(cn_v,sides)]
#    muv = [mu[side] for mu,side in zip(mu_v,sides)]
#    cn_states = [[[cr_i[0][0],cv_i[0][0],mv_i[0][0]],[cr_i[0][1],cv_i[0][1],mv_i[0][1]], \
#                  [cr_i[1][0],cv_i[1][0],mv_i[1][0]],[cr_i[1][1],cv_i[1][1],mv_i[1][1]]] for cr_i,cv_i,mv_i in zip(cnr,cnv,muv)]    
#    return [filter_cns(state) for state in cn_states]

def calc_lik(combo,si,di,phi_i,pi):
    pvs = [ get_pv(phi_i,c[0],c[1],c[2],pi) for c in combo ]
    lls = [ pm.binomial_like(si,di,pvs[i]) for i,c in enumerate(combo)]
    return (pvs, lls)

def get_most_likely_cn_states(cn_states,s,d,phi,pi):
    '''
    Obtain the copy-number states which maximise the binomial likelihood
    of observing the supporting read depths at each variant location
    '''
    cn_ll = [ calc_lik(cn_states[i],s[i],d[i],phi[i],pi)  for i in range(len(cn_states)) ]
    most_likely_pv = [ cn_lik[0][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]] for i,cn_lik in enumerate(cn_ll)]
    most_likely_cn = [ cn_states[i][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]] for i,cn_lik in enumerate(cn_ll)]
    #TODO: currently using precision fudge factor to get around 
    #0 probability errors when pv = 1 - look into this bug more
    return most_likely_cn, np.array(most_likely_pv)-0.00000001

def cluster(sv_df,snv_df,pi,rlen,insert,ploidy,iters,burn,thin,beta,use_map,Ndp=param.clus_limit):
    '''
    clustering model using Dirichlet Process
    '''
    pl = ploidy
    sup, dep, cn_states = [], [], []
    #sides = [] #only relevant for SVs
    cn_states = []
    Nvar = 0
    
    if len(sv_df)>0 and len(snv_df)>0:
        print("Coclustering SVs & SNVs...")
        sup,dep,cn_states,Nvar = get_snv_vals(snv_df)
        sv_sup,sv_dep,sv_cn_states,sv_Nvar = get_sv_vals(sv_df,rlen)
        
        sup = np.append(sup,sv_sup)
        dep = np.append(dep,sv_dep)
        cn_states = np.append(cn_states,sv_cn_states)
        Nvar = Nvar + sv_Nvar
    elif len(snv_df)>0:
        print("Clustering SNVs...")
        sup,dep,cn_states,Nvar = get_snv_vals(snv_df)
        #sides = np.zeros(Nvar,dtype=int)
    else:
        print("Clustering SVs...")
        sup,dep,cn_states,Nvar = get_sv_vals(sv_df,rlen)

    sens = 1.0 / ((pi/float(pl))*np.mean(dep))
    alpha = pm.Gamma('alpha',0.9,1/0.9,value=2)
    #beta = pm.Gamma('beta',param.beta_shape,param.beta_rate)
    #beta = pm.Gamma('beta',1,10**(-7))
    #beta = pm.Uniform('beta', 0.01, 1, value=0.1)
    #print("Beta value:%f"%beta)

    h = pm.Beta('h', alpha=1, beta=alpha, size=Ndp)
    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1.0-h[:i]) for i,u in enumerate(h)]
        value /= np.sum(value)
        #value[-1] = 1.0-sum(value[:-1])
        return value

    z = pm.Categorical('z', p=p, size=Nvar, value=np.zeros(Nvar))
    #phi_init = (np.mean(sup/dep)/pi)*2
    phi_k = pm.Uniform('phi_k', lower=sens, upper=2, size=Ndp)#, value=[phi_init]*Ndp)
    
    @pm.deterministic
    def p_var(z=z,phi_k=phi_k):
#        ml_cn = [get_most_likely_cn(cn_ri[side],cn_vi[side],mu_vi[side],si,di,phi,pi) \
#                for cn_ri,cn_vi,mu_vi,si,di,phi,side in zip(cn_r,cn_v,mu_v,sup,dep,phi_k[z],sides)]
#        cn_rn = [m[0] for m in ml_cn]
#        cn_vn = [m[1] for m in ml_cn]
#        mu_vn = [m[2] for m in ml_cn]
#
#        return  get_pv(phi_k[z],cn_rn,cn_vn,mu_vn,pi)
        most_lik_cn_states, pvs = get_most_likely_cn_states(cn_states,sup,dep,phi_k[z],pi)
        return pvs
    
    cbinom = pm.Binomial('cbinom', dep, p_var, observed=True, value=sup)

    model = pm.Model([alpha,h,p,phi_k,z,p_var,cbinom])
    mcmc = fit_and_sample(model,iters,burn,thin,use_map)
    return mcmc
