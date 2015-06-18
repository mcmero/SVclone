import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import IPython
import pymc as pm
import colorsys
import ipdb

from IPython.core.pylabtools import figsize
from . import parameters as param

def gen_new_colours(N):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def plot_cluster_hist(clusters,assignments,df,pl,pi):
    fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(12.5,8))
    RGB_tuples = gen_new_colours(len(clusters))
    
    for idx,clus in enumerate(clusters):
        t1 = df[np.array(assignments)==clus]
        s1 = t1.support.values
        n1 = map(np.mean,zip(t1.norm1.values,t1.norm2.values))
        axes[0].set_title("Clusters post-cluster merge: Cell fractions (raw VAFs purity-ploidy-adjusted)")
        axes[0].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        axes[1].set_title("Raw VAFs")
        axes[1].hist(s1/(n1+s1),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
    
def plot_clusters(center_trace,clusters,assignments,df,pl,pi):
    fig, axes = plt.subplots(3, 1, sharex=False, sharey=False, figsize=(12.5,11))

    RGB_tuples = gen_new_colours(len(clusters))
    
    axes[0].set_ylim([0,1])
    axes[0].set_title("Trace of $\phi_k$")
    
    for idx,clus in enumerate(clusters): 
        axes[0].plot(center_trace[:, clus], label="trace of center %d" % clus, c=RGB_tuples[idx], lw=1)
        # print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))
    
    leg = axes[0].legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

    for idx,clus in enumerate(clusters):
        t1 = df[np.array(assignments)==clus]
        s1 = t1.support.values
        n1 = map(np.mean,zip(t1.norm1.values,t1.norm2.values))
        axes[1].set_title("Unmerged clusters: Cell fractions (raw VAFs purity-ploidy-adjusted)")
        axes[1].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        #print 'mean cluster VAF: %f' % (np.mean(s1/(n1+s1)/pi))
        axes[2].set_title("Raw VAFs")
        axes[2].hist(s1/(n1+s1),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

def get_cn_mu_v(cn):
    cn_v = [0.,0.]
    mu_v = [0.,0.]

    c = cn.split(',')
    if len(c)<2:
        return tuple(cn_v),tuple(mu_v)
    
    c = map(float,c)
    cn_t = c[0]+c[1]
    cn_v[0] = cn_t
    mu_v[0] = float(c[0])/(c[0]+c[1])

    if c[0]!=c[1] and c[1]>0: 
        cn_v[1] = cn_t
        mu_v[1] = float(c[1])/(c[0]+c[1])
    
    return tuple(cn_v),tuple(mu_v)

def get_allele_combos_tuple(c):
    cn_r = [tuple([0.,0.]),tuple([0.,0.])]
    cn_v = [tuple([0.,0.]),tuple([0.,0.])]
    mu_v = [tuple([0.5,0.5]),tuple([0.5,0.5])]
    
    if len(c) < 1:
        return cn_r,cn_v,mu_v

    cn_v[0],mu_v[0] = get_cn_mu_v(c[0])
    cn1_tmp = map(float,c[0].split(',')) if len(c[0])>1 else c[0]
    if len(c) > 1:
        cn_v[1],mu_v[1] = get_cn_mu_v(c[1])
        cn2_tmp = map(float,c[1].split(','))
        
        cnr1_tmp = cn1_tmp[0]+cn1_tmp[1]
        cnr2_tmp = cn2_tmp[0]+cn2_tmp[1]
        cn_r[1],cn_r[0] = tuple([cnr1_tmp,cnr1_tmp]),tuple([cnr2_tmp,cnr2_tmp])
    else:
        cn_r[0],cn_r[1] = tuple([2.,2.]),tuple([2.,2.])
    return tuple(cn_r),tuple(cn_v),tuple(mu_v)

def get_read_vals(df):
    #n = np.array(df.norm_mean)
    n = [np.array(df.norm1.values),np.array(df.norm2.values)]
    s = np.array(df.bp1_split.values+df.bp2_split.values)
    d = np.array(df.spanning.values)
    cn_r,cn_v,mu_v = [],[],[]

    for idx,sv in df.iterrows():
        cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
        cnr_bp1,cnv_bp1,mu_bp1 = get_allele_combos_tuple(cn_tmp[0])
        cnr_bp2,cnv_bp2,mu_bp2 = get_allele_combos_tuple(cn_tmp[1])
        cn_r.append(tuple([cnr_bp1,cnr_bp2]))
        cn_v.append(tuple([cnv_bp1,cnv_bp2]))
        mu_v.append(tuple([mu_bp1,mu_bp2]))

    return n,d,s,cn_r,cn_v,mu_v

def fit_and_sample(model, iters, burn, thin):
    # map_ = pm.MAP(model)
    # map_.fit()
    mcmc = pm.MCMC(model)
    mcmc.sample(iters, burn=burn, thin=thin)
    return mcmc

def cluster(df,pi,rlen,insert,ploidy,iters,burn,thin,beta,Ndp=param.clus_limit):
    '''
    inital clustering using Dirichlet Process
    '''    
    pl = ploidy
    n,d,s,cn_r,cn_v,mu_v = get_read_vals(df)
    #dep = [np.array(n[0]+d+s,dtype=int),np.array(n[1]+d+s,dtype=int)]
    n_max = np.array(map(max,zip(n[0],n[1])))
    sup = np.array(d+s,dtype=int)
    dep = np.array(n_max+sup,dtype=int)
    Nsv = len(sup)
    sens = 1.0 / ((pi/pl)*np.average(dep))
    #print(sens)
    
    beta = pm.Gamma('beta',0.9,1/0.9) 
    #beta = pm.Gamma('beta',param.beta_shape,param.beta_rate) 
    #beta = pm.Gamma('beta',1,10**(-7)) 
    #beta = pm.Uniform('beta', 0.01, 1, value=0.1)
    #print("Beta value:%f"%beta)
    
    h = pm.Beta('h', alpha=1, beta=beta, size=Ndp)
    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1.0-h[:i]) for i,u in enumerate(h)]
        value /= np.sum(value)
        #value[-1] = 1.0-sum(value[:-1])
        return value

    z = pm.Categorical('z', p=p, size=Nsv, value=np.zeros(Nsv))
    #phi_init = (np.mean((s+d)/(n+s+d))/pi)*2
    phi_k = pm.Uniform('phi_k', lower=sens, upper=1, size=Ndp)#, value=[phi_init]*Ndp)

    def get_pv(phi,cn_r,cn_v,mu):
        pn =  (1.0 - pi) * 2            #proportion of normal reads coming from normal cells
        pr =  pi * (1.0 - phi) * cn_r   #proportion of normal reads coming from other clusters
        pv =  pi * phi * cn_v           #proportion of variant reads coming from this cluster
        
        norm_const = pn + pr + pv
        pv = pv / norm_const  

        return pv * mu

    def get_most_likely_cn(cn_r,cn_v,mu_v,si,di,phi):
        llik = []
        vals = []
        # for each side of the break
        for side in range(2):
            # for each subclone
            for cn_rx, cn_vx, mu_vx in zip(cn_r[side],cn_v[side],mu_v[side]): 
                # for each allelic state
                for cn_ri, cn_vi, mu_vi in zip(cn_rx,cn_vx,mu_vx):                    
                    vals.append([cn_ri,cn_vi,mu_vi])
                    if cn_vi==0 or mu_vi==0: 
                        llik.append(float('nan'))
                    else:
                        pv = get_pv(phi,cn_ri,cn_vi,mu_vi)
                        temp = pm.binomial_like(si,di,pv)
                        llik.append(temp)
        
        idx = np.where(np.array(llik)==max(llik))[0]
        if len(idx)==0:
            return [2.,2.,0.5]
        elif len(idx) >1: 
            return vals[idx[0]]
        else:
            return vals[idx]

    @pm.deterministic
    def p_var(z=z,phi_k=phi_k):
        ml_cn = [get_most_likely_cn(cn_ri,cn_vi,mu_vi,si,di,phi) for cn_ri,cn_vi,mu_vi,si,di,phi in zip(cn_r,cn_v,mu_v,sup,dep,phi_k[z])]   
        cn_rn = [m[0] for m in ml_cn]
        cn_vn = [m[1] for m in ml_cn]
        mu_vn = [m[2] for m in ml_cn]    
        #ipdb.set_trace()
        
        return get_pv(phi_k[z],cn_rn,cn_vn,mu_vn)

#    @pm.deterministic
#    def p_var(z=z,phi_k=phi_k):
#        pn =  (1.0 - pi) * 2                #proportion of normal reads coming from normal cells
#        pr =  pi * (1.0 - phi_k[z]) * pl    #proportion of normal reads coming from other clusters
#        pv =  pi * phi_k[z] * (1.0/pl)      #proportion of variant reads coming from this cluster
#        pvn = pi * phi_k[z] * ((pl-1.0)/pl) #proportion of normal reads coming from this cluster
#        
#        norm_const = pn + pr + pv + pvn
#        pv = pv / norm_const    
#        
#        return pv
    #ipdb.set_trace() 
    cbinom = pm.Binomial('cbinom', dep, p_var, observed=True, value=sup)

    model = pm.Model([beta,h,p,phi_k,z,p_var,cbinom])
    mcmc = fit_and_sample(model,iters,burn,thin)
    return mcmc

