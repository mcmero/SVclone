import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import IPython
import pymc as pm
import scipy.optimize as op
import colorsys
import ipdb
import sys
from ete2 import Tree
from scipy.stats import itemfreq
from IPython.core.pylabtools import figsize

sc_thresh = 0.05 #subclone size threshold

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def plot_clusters(center_trace,npoints,clusters):
    figsize(12.5, 9)
    plt.subplot(311)
    lw = 1

    N = len(clusters)
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    
    for idx,clus in enumerate(clusters): 
        plt.plot(center_trace[:, clus], label="trace of center 0", c=RGB_tuples[idx], lw=lw)
        print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))
    
    plt.title("Trace of $\phi_k$")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

def cluster(df,pi,rlen,insert,Ndp=35,mcmc_samples=100000,burn=5000):
    n1 = df.norm1.values
    n2 = df.norm2.values
    n = np.array(df.norm_mean)
    d = np.array(df.bp1_split.values+df.bp2_split.values)
    s = np.array(df.spanning.values)

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
    mcmc.sample(mcmc_samples, burn=burn)
    return mcmc

def infer_tree(df,pi,rlen,insert,out):
    mcmc = cluster(df,pi,rlen,insert)
    npoints = len(df.spanning.values)

    # assign points to highest probabiliy cluster
    z_trace = mcmc.trace('z')[:]
    cm = [np.bincount(z_trace[:,i]) for i in range(npoints)]
    cx = [index_max(c) for c in cm]
    cf = np.bincount(cx)
    ii = np.nonzero(cf)[0]

    # cluster distribution, filter clones below threshold
    cd = zip(ii,cf[ii])
    print("Initial clustering")
    print(cd)
    flt_clones = sum(cf)*sc_thresh<cf
    ii = ii[flt_clones]
    cf = cf[flt_clones]
    
    clusters = np.array([np.where(np.array(cx)==i)[0] for i in ii])
    if len(clusters) < 1:
        print("Could not converge on any major clusters! Exiting.")
        sys.exit
    
    center_trace = mcmc.trace("phi_k")[:]
    plot_clusters(center_trace,npoints,ii)
    
    s_order = np.argsort(cf)[::-1]
    clus_ids = ii[s_order]
    clusters = clusters[s_order] 
    phis = [np.mean(center_trace[2000:,cid]) for cid in clus_ids]
    clus_idx = range(len(clus_ids))

    print("Filtered clusters")
    print(zip(clus_ids,cf))

    if (phis[0]*2.)>(1+sc_thresh):
        print("Major clonal cluster VAF exceeds 0.5 - either germline breaks were not " + \
              "adequately filtered, or provided copy-numbers are inaccurate. Please double-check and rerun")
        sys.exit
    
    # 1 cluster - no tree building possible
    if len(clusters)==1:
        print("Single sample contains 1 clonal cluster, no tree building attempted")
        t = Tree()
        t.add_child(name=clus_idx[0])
        with open('%s.txt'%out,'w') as outf:
            outf.write("id\tphi\tnum_children\tnum_svs\tsvs\n")            
            svs = ['%s:%s|%s:%s' % x for x in zip(svs.bp1_chr,svs.bp1_pos,svs.bp2_chr,svs.bp2_pos)]
            outf.write('%d\t%f\t%d\t%d\t%s'%(clus_idx[0],phis[0],0,len(clusters[0]),','.join(svs)))
            outf.write(t.get_ascii(show_internal=True))
    else:
        # can we use the mutual-exclusion principle?
        if phis[0]<(sum(phis[1:])+sc_thresh):
            # algorithm to deal with this
            print("Need to implement pigeon-hole principle")
        else:
            # we can't be sure whether the phylogeny is branching or linear
            t = Tree()
            with open('%s.txt'%out,'w') as outf:
                outf.write("id\tphi\tnum_children\tnum_svs\tsvs\n")
                lin_phis = [phi-phi[idx+1] for idx,phi in enumerate(phis) if idx+1>=len(phis)]
                lin_phis.extend(phis[len(phis)-1])
                
                for idx,c in enumerate(clusters):
                    t.add_child(name=clus_idx[idx])
                    svs = ['%s:%s|%s:%s' % x for x in zip(svs.bp1_chr,svs.bp1_pos,svs.bp2_chr,svs.bp2_pos)]
                    outf.write('%d\t%f\t%d\t%d\t%s'%(clus_idx[idx],lin_phis[idx],idx,len(clusters[idx]),','.join(svs)))

                outf.write(t.get_ascii(show_internal=True))
