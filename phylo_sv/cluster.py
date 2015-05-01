import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import IPython
import pymc as pm
import scipy.optimize as op
from scipy.stats import itemfreq
from IPython.core.pylabtools import figsize

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def plot_clusters(mcmc,npoints):
    z_trace = mcmc.trace('z')[:]
    # assign points to highest probabiliy cluster
    cm = [np.bincount(z_trace[:,i]) for i in range(npoints)]
    cx = [index_max(c) for c in cm]
    cf = np.bincount(cx)
    ii = np.nonzero(cf)[0]

    # cluster distribution
    cd = zip(ii,cf[ii])
    print(cd)

    figsize(12.5, 9)
    plt.subplot(311)
    lw = 1
    center_trace = mcmc.trace("phi_k")[:]

    clus = np.where(cf==max(cf))[0]
    clus_sc = cf.argsort()[::-1][1]

    plt.plot(center_trace[:, clus], label="trace of center 0", c="#348ABD", lw=lw)
    plt.plot(center_trace[:, clus_sc], label="trace of center 0", c="#A60628", lw=lw)
    plt.title("Trace of $\phi_k$")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

    print('phi_1 ~ %f' % np.mean(center_trace[2000:,clus]))
    print('phi_2 ~ %f' % np.mean(center_trace[2000:,clus_sc]))

def cluster(df,pi,rlen,insert):
    n1 = df.norm1.values
    n2 = df.norm2.values
    n = np.array(df.norm_mean)
    d = np.array(df.bp1_split.values+df.bp2_split.values)
    s = np.array(df.spanning.values)

    Ndp = 25
    beta = pm.Uniform('beta', lower=0.01, upper=1)
    h = pm.Beta('h', alpha=1, beta=beta, size=Ndp)

    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1-h[:i]) for i,u in enumerate(h)]
        value /= np.sum(value)   
        return value

    z = pm.Categorical('z', p=p, size=len(d), value=np.zeros(len(d)))
    phi_k = pm.Uniform('phi_k', lower=0, upper=pi, size=Ndp)

    @pm.deterministic
    def smu(z=z, phi_k=phi_k):
        return (rlen/(rlen+0.5*insert))*(phi_k[z]*pi)

    @pm.deterministic
    def dmu(z=z, phi_k=phi_k):
        return (insert/(2*rlen+insert))*(phi_k[z]*pi)

    @pm.deterministic
    def nmu(z=z, phi_k=phi_k):
        return (1-phi_k[z])+(0.5*(1-pi))
                
    sp = pm.Poisson('sp', mu=smu, observed=True, value=s)
    dp = pm.Poisson('dp', mu=dmu, observed=True, value=d)
    normp = pm.Poisson('normp', mu=nmu, observed=True, value=n)

    model = pm.Model([p,z,dp,sp,normp,phi_k])

    mcmc = pm.MCMC(model)
    mcmc.sample(100000, burn=5000)
    plot_clusters(mcmc,len(d))

    
