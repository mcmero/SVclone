import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import IPython
import pymc as pm
import colorsys
import ipdb

from IPython.core.pylabtools import figsize
from . import parameters as param

def plot_clusters(center_trace,npoints,clusters):
    figsize(12.5, 9)
    plt.subplot(311)
    lw = 1

    N = len(clusters)
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    
    for idx,clus in enumerate(clusters): 
        plt.plot(center_trace[:, clus], label="trace of center 0", c=RGB_tuples[idx], lw=lw)
        # print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))
    
    plt.title("Trace of $\phi_k$")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

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

#def recluster(df,pi,rlen,insert,ploidy,iters,burn,thin,dump=True):
#    '''
#    reclusters group without Dirichlet (assuming only one group exists)
#    '''
#    print("Reclustering with %d SVs"%len(df))
#    n,d,s = get_read_vals(df)
#
#    phi_k = pm.Uniform('phi_k', lower=0, upper=1)
#
#    @pm.deterministic
#    def smu(phi_k=phi_k):
#        return (rlen/(rlen+0.5*insert))*((phi_k/ploidy)/pi)
#
#    @pm.deterministic
#    def dmu(phi_k=phi_k):
#        return (insert/(2*rlen+insert))*((phi_k/ploidy)/pi)
#
#    @pm.deterministic
#    def nmu(phi_k=phi_k):
#        return (1 - pi) + (pi * (1-phi_k))
#                
#    sp = pm.Poisson('sp', mu=smu, observed=True, value=s)
#    dp = pm.Poisson('dp', mu=dmu, observed=True, value=d)
#    normp = pm.Poisson('normp', mu=nmu, observed=True, value=n)
#
#    model = pm.Model([dp,sp,normp,phi_k])
#    mcmc = fit_and_sample(model,iters,burn,thin)
#   
#    return mcmc.trace("phi_k")[:]
#    # center_trace = center_trace[len(center_trace)/4:]
#    # phi = np.mean(center_trace[:])
#    # return center_trace

def cluster(df,pi,rlen,insert,ploidy,iters,burn,thin,beta,Ndp=param.clus_limit):
    '''
    inital clustering using Dirichlet Process
    '''    
    pl = ploidy
    #pl = 2
    n,d,s = get_read_vals(df)
    dep = np.array(n+d+s,dtype=int)
    sup = np.array(d+s,dtype=int)
    Nsv = len(sup)

    beta = pm.Gamma('beta',0.8,1/0.8) 
    #beta = pm.Uniform('beta', 0.01, 1, value=0.1)
    #print("Beta value:%f"%beta)

    h = pm.Beta('h', alpha=1, beta=beta, size=Ndp)
    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1-h[:i]) for i,u in enumerate(h)]
        #value /= np.sum(value)
        value[-1] = 1-sum(value[:-1])
        return value

    z = pm.Categorical('z', p=p, size=Nsv, value=np.zeros(Nsv))
    #phi_init = (np.mean((s+d)/(n+s+d))/pi)*2
    phi_k = pm.Uniform('phi_k', lower=0, upper=1, size=Ndp)#, value=[phi_init]*Ndp)

    @pm.deterministic
    def mu(z=z,phi_k=phi_k):
        pn = (1 - pi) * 2 #proportion of normal reads coming from normal cells
        pr = pi * (1 - phi_k[z]) * pl #proportion of normal reads coming from other clusters
        pv = pi * phi_k[z] / pl #proportion of variant reads coming from this cluster
        pvn = pi * phi_k[z] * (pl-1) / pl #proportion of normal reads coming from this cluster
        
        norm_const = pn + pr + pv + pvn
        pv = pv / norm_const    
        
        return (pv)

    cbinom = pm.Binomial('cbinom', dep, mu, observed=True, value=sup)

    model = pm.Model([h,p,phi_k,z,mu,cbinom])
    mcmc = fit_and_sample(model,iters,burn,thin)
    return mcmc

