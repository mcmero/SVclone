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
    print("Clustering with %d SVs"%len(df))
    n,d,s = get_read_vals(df)

    beta = pm.Gamma('beta', 8, 1/0.015,value=0.2)
    #print("Beta value:%f"%beta)
    h = pm.Beta('h', alpha=1, beta=beta, size=Ndp)
    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1-h[:i]) for i,u in enumerate(h)]
        value /= np.sum(value)   
        return value

    z = pm.Categorical('z', p=p, size=len(d), value=np.zeros(len(d)))
    #phi_init = np.mean((s+d)/(n+s+d))/pi
    #phi_k = pm.Gamma('phi_k', alpha=3, beta=1/0.15, size=Ndp)#, value=[phi_init]*Ndp)
    phi_k = pm.Uniform('phi_k', lower=0.01, upper=1, size=Ndp)#, value=[phi_init]*Ndp)

    @pm.deterministic
    def smu(z=z, phi_k=phi_k):
        #return (rlen/(rlen+0.5*insert))*((phi_k[z]/ploidy)*pi)
        return ( (rlen/(rlen+0.5*insert))*(phi_k[z]*pi) )#*(ploidy/2)

    @pm.deterministic
    def dmu(z=z, phi_k=phi_k):
        #return (insert/(2*rlen+insert))*((phi_k[z]/ploidy)*pi)
        return ( (insert/(2*rlen+insert))*(phi_k[z]*pi) )#*(ploidy/2)

    @pm.deterministic
    def nmu(z=z, phi_k=phi_k):
        return  ( 2*(1 - pi) + (ploidy-1)*pi*(phi_k[z]) + ploidy*pi*(1.05-phi_k[z]) )#*(ploidy/2)
        #return (1 - pi) + (pi * (1-(phi_k[z]/ploidy)))
        #return (1 - pi) + (pi * (1-(phi_k[z])))
                
    sp = pm.Poisson('sp', mu=smu, observed=True, value=s)
    dp = pm.Poisson('dp', mu=dmu, observed=True, value=d)
    normp = pm.Poisson('normp', mu=nmu, observed=True, value=n)

    model = pm.Model([phi_k,beta,h,p,z,normp,smu,dmu])
    mcmc = fit_and_sample(model,iters,burn,thin)
    return mcmc

