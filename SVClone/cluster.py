import numpy as np
import pandas as pd
import pymc as pm
import ipdb

from operator import methodcaller
from . import parameters as param

#def get_cn_mu_v(cn):
#    cn_v = [0.,0.]
#    mu_v = [0.,0.]
#
#    c = cn.split(',')
#    if len(c)<2:
#        return tuple(cn_v),tuple(mu_v)
#    if c[0]>1 or c[1]>1:
#        ipdb.set_trace()
#    
#    c = map(float,c)
#    cn_t = float(c[0]+c[1])
#    cn_v[0] = float(cn_t)
#    mu_v[0] = c[0]/cn_t if cn_t!=0 else 0.
#
#    if c[0]!=c[1] and c[1]>0:
#        cn_v[1] = cn_t
#        mu_v[1] = c[1]/cn_t if cn_t!=0 else 0.
#    
#    return tuple(cn_v),tuple(mu_v)

def add_copynumber_combos(combos, var_maj, var_min, ref_cn):
    '''
    ref_cn = total reference (non-variant) copy-number
    var_total = total variant copy-number 
    mu_v = copies containing variant / var_total
    
    possible copynumber states for variant are:
        - major / var total
        - minor / var_total, and 
        - 1 / var_total (if neither alleles are 1)
    '''
    var_total = float(var_maj + var_min)
    if var_total == 0.:
        mu_v = 0.
        #combos.append([ref_cn, var_total, mu_v])
    #elif var_maj == 1.:
    #    mu_v = 1. / var_total
    #    combos.append([ref_cn, var_total, mu_v])
    #elif var_maj == var_min:
    #    mu_v = var_maj / var_total
    #    combos.append([ref_cn, var_total, mu_v])
    else:
        for i in range(1,int(var_maj+1)):
            mu_v = 1.0*i / var_total
            combos.append([ref_cn, var_total, mu_v])

    # add combos for the minor allele being the variant allele
    #if var_maj != var_min:
        # do the exact same thing as above, swapping
        # the major allele for the minor
    #    mu_v = var_min / var_total if var_total != 0 else 0.
    #    combos.append([ref_cn, var_total, mu_v])

        #if var_min > 1:
        #    mu_v = 1. / var_total
        #    combos.append([ref_cn, var_total, mu_v])
    #print(combos)
    return combos

def get_allele_combos(cn):
    combos = []

    if len(cn) == 0 or cn[0]=='':
        return combos

    if len(cn)>1:
        # split subclonal copy-numbers
        c1 = map(float,cn[0].split(','))
        c2 = map(float,cn[1].split(','))
        major1, minor1, total1 = c1[0], c1[1], c1[0]+c1[1]
        major2, minor2, total2 = c2[0], c2[1], c2[0]+c2[1]
        
        # generate copynumbers for each subclone being the potential variant pop 
        combos = add_copynumber_combos(combos, major1, minor1, total2)
        combos = add_copynumber_combos(combos, major2, minor2, total1)
    else:
        c = map(float,cn[0].split(','))
        major, minor = c[0], c[1]
        combos = add_copynumber_combos(combos, major, minor, major + minor)

    return filter_cns(combos)

def get_sv_allele_combos(sv):
    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
    combos_bp1 = get_allele_combos(cn_tmp[0])
    combos_bp2 = get_allele_combos(cn_tmp[1])

    return tuple([combos_bp1,combos_bp2])

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
    pn =  (1.0 - pi) * 2.0            #proportion of normal reads coming from normal cells
    #pr =  pi * rem * cn_r           #proportion of normal reads coming from other (reference) clusters
    pv =  pi * cn_v           #proportion of variant + normal reads coming from this (the variant) cluster
    norm_const = pn + pv
    pv = pv / norm_const

    return phi * pv * mu

def filter_cns(cn_states):
    cn_str = [','.join(map(str,cn)) for cn in cn_states if cn[2]!=0 and cn[1]!=0]
    cn_str = np.unique(np.array(cn_str))
    return [map(float,cn) for cn in map(methodcaller('split',','),cn_str)]

def calc_lik(combo,si,di,phi_i,pi):
    pvs = [ get_pv(phi_i,c[0],c[1],c[2],pi) for c in combo ]
    lls = [ pm.binomial_like(si,di,pvs[i]) for i,c in enumerate(combo)]
    return (pvs, lls)

def get_most_likely_cn_states(cn_states,s,d,phi,pi):
    '''
    Obtain the copy-number states which maximise the binomial likelihood
    of observing the supporting read depths at each variant location
    '''
    
    def get_most_likely_pv(cn_lik):
        if len(cn_lik[0])>0:
            return cn_lik[0][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]]
        else:
            return 0.0000001

    def get_most_likely_cn(cn_lik,i):
        if len(cn_lik[0])>0:
            return cn_states[i][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]] 
        else:
            return [float('nan'), float('nan'), float('nan')]
    
    cn_ll = [ calc_lik(cn_states[i],s[i],d[i],phi[i],pi) for i in range(len(cn_states)) ]
    most_likely_pv = [ get_most_likely_pv(cn_lik) for i,cn_lik in enumerate(cn_ll)]
    most_likely_cn = [ get_most_likely_cn(cn_lik,i) for i,cn_lik in enumerate(cn_ll)]

    #TODO: currently using precision fudge factor to get around 
    #0 probability errors when pv = 1 - look into this bug more
    return most_likely_cn, np.array(most_likely_pv)-0.00000001

def cluster(sup,dep,cn_states,Nvar,sparams,cparams,Ndp=param.clus_limit):
    '''
    clustering model using Dirichlet Process
    '''
    sens = 1.0 / ((sparams['pi']/float(sparams['ploidy']))*np.mean(dep))
    beta_a, beta_b, beta_init = map(lambda x: float(eval(x)), cparams['beta'].split(','))
    alpha = pm.Gamma('alpha',beta_a,beta_b)#,value=beta_init)
    print("Dirichlet concentration gamma values: alpha = %f, beta= %f, init = %f" % (beta_a, beta_b, beta_init))

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
        most_lik_cn_states, pvs = get_most_likely_cn_states(cn_states,sup,dep,phi_k[z],sparams['pi'])
        return pvs
    
    cbinom = pm.Binomial('cbinom', dep, p_var, observed=True, value=sup)

    model = pm.Model([alpha,h,p,phi_k,z,p_var,cbinom])
    mcmc = fit_and_sample(model,cparams['n_iter'],cparams['burn'],cparams['thin'],cparams['use_map'])
    return mcmc
