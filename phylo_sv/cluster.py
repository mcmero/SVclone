import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import IPython
import pymc as pm
import colorsys
import ipdb

from IPython.core.pylabtools import figsize
from . import parameters as param

def plot_clusters(center_trace,clusters,assignments,df,pl,pi):
    #TODO: incorporate merged cluster traces
    fig, axes = plt.subplots(3, 1, sharex=False, sharey=False, figsize=(12.5,11))

    N = len(clusters)
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    
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
        axes[1].set_title("Cell fractions (raw VAFs purity-ploidy-adjusted)")
        axes[1].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        #print 'mean cluster VAF: %f' % (np.mean(s1/(n1+s1)/pi))
        axes[2].set_title("Raw VAFs")
        axes[2].hist(s1/(n1+s1),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

def get_read_vals(df):
    n = np.array(df.norm_mean)
    s = np.array(df.bp1_split.values+df.bp2_split.values)
    d = np.array(df.spanning.values)
    return n,d,s

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
    n,d,s = get_read_vals(df)
    dep = np.array(n+d+s,dtype=int)
    sup = np.array(d+s,dtype=int)
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

    @pm.deterministic
    def p_var(z=z,phi_k=phi_k):
        pn =  (1.0 - pi) * 2                #proportion of normal reads coming from normal cells
        pr =  pi * (1.0 - phi_k[z]) * pl    #proportion of normal reads coming from other clusters
        pv =  pi * phi_k[z] * (1.0/pl)      #proportion of variant reads coming from this cluster
        pvn = pi * phi_k[z] * ((pl-1.0)/pl) #proportion of normal reads coming from this cluster
        
        norm_const = pn + pr + pv + pvn
        pv = pv / norm_const    
        
        return pv
        
    cbinom = pm.Binomial('cbinom', dep, p_var, observed=True, value=sup)

    model = pm.Model([beta,h,p,phi_k,z,p_var,cbinom])
    mcmc = fit_and_sample(model,iters,burn,thin)
    return mcmc

