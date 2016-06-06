import numpy as np
import pandas as pd
import pymc3 as pm
import pymc as pm2
import math
import ipdb
import theano
import theano.tensor as t

from theano.compile.ops import as_op
from scipy import stats, optimize

def add_copynumber_combos(combos, var_maj, var_min, ref_cn, frac):
    '''
    ref_cn = total reference (non-variant) copy-number
    var_total = total variant copy-number 
    mu_v = copies containing variant / var_total
    
    possible copynumber states for variant are:
        - (1 .. major) / var total
    '''
    var_total = float(var_maj + var_min)
    if var_total == 0.:
        mu_v = 0.
    else:
        for i in range(1, int(var_maj+1)):
            mu_v = 1.0*i / var_total
            combos.append([ref_cn, var_total, mu_v, frac])

    return combos

def get_allele_combos(cn):
    combos = []

    if len(cn) == 0 or cn[0]=='':
        return combos

    if len(cn)>1:
        # split subclonal copy-numbers
        c1 = [float(c) for c in cn[0].split(',')]
        c2 = [float(c) for c in cn[1].split(',')]
        major1, minor1, total1, frac1 = c1[0], c1[1], c1[0]+c1[1], c1[2]
        major2, minor2, total2, frac2 = c2[0], c2[1], c2[0]+c2[1], c2[2]

        # generate copynumbers for each subclone being the potential variant pop 
        combos = add_copynumber_combos(combos, major1, minor1, total2, frac1)
        combos = add_copynumber_combos(combos, major2, minor2, total1, frac2)
    else:
        c = [float(c) for c in cn[0].split(',')]
        major, minor, frac = c[0], c[1], c[2]
        combos = add_copynumber_combos(combos, major, minor, major + minor, frac)

    return filter_cns(combos)

def get_sv_allele_combos(sv):
    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
    combos_bp1 = get_allele_combos(cn_tmp[0])
    combos_bp2 = get_allele_combos(cn_tmp[1])

    return tuple([combos_bp1,combos_bp2])

def get_pv(phi, cn_r, cn_v, mu, cn_f, pi):
    pn = (1.0 - pi) * 2.0       #proportion of normal reads coming from normal cells
    pr = pi * cn_r * (1. - cn_f) if cn_f < 1 else 0 # incorporate ref population CNV fraction if present
    pv = pi * cn_v * cn_f       #proportion of variant + normal reads coming from this (the variant) cluster
    norm_const = pn + pv + pr
    pv = pv / norm_const

    return phi * pv * mu

def filter_cns(cn_states):
    cn_str = [','.join(map(str,cn)) for cn in cn_states if cn[2]!=0 and cn[1]!=0]
    cn_str = np.unique(np.array(cn_str))
    cn_str = [x.split(',') for x in cn_str]
    return [list(map(float,cn)) for cn in cn_str]

def calc_lik(combo,si,di,phi_i,pi):
    pvs = [get_pv(phi_i,c[0],c[1],c[2],c[3],pi) for c in combo]
    lls = [pm2.binomial_like(si,di,pvs[i]) for i,c in enumerate(combo)]
    #currently using precision fudge factor to get 
    #around 0 probability errors when pv = 1
    #TODO: investigate look this bug more
    return (pvs, np.array(lls)-0.00000001)

def get_probs(var_states,s,d,phi,pi):
    #probs = [(1/float(len(var_states)))]*len(var_states)
    llik = calc_lik(var_states,s,d,phi,pi)[1]
    probs = get_probs_from_llik(llik)
    probs = ','.join(map(lambda x: str(round(x,4)),probs))
    return probs

def get_probs_from_llik(cn_lik):
    probs = np.array([1.])
    if len(cn_lik) > 1:        
        probs = [math.exp(x) for x in cn_lik]
        probs = np.array(probs)/sum(probs)
    return probs

def get_most_likely_cn_states(cn_states, s, d, phi, pi):
    '''
    Obtain the copy-number states which maximise the binomial likelihood
    of observing the supporting read depths at each variant location
    '''
    
    def get_most_likely_pv(cn_lik):
        if len(cn_lik[0])>0:
            return cn_lik[0][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]]
        else:
            return 0.0000001

    def get_most_likely_cn(cn_states, cn_lik, i):
        '''
        use the most likely phi state, unless p < cutoff when compared to the 
        most likely clonal (phi=1) case (log likelihood ratio test) 
        - in this case, pick the most CN state with the highest clonal likelihood
        '''
        pval_cutoff = 0.01
        cn_lik_clonal, cn_lik_phi = cn_lik

        #reinstitute hack - uncomment below
        #cn_lik_phi = cn_lik_clonal
        if len(cn_lik_clonal)==0:
            return [float('nan'), float('nan'), float('nan')]
   
        # log likelihood ratio test; null hypothesis = likelihood under phi
        LLR   = 2 * (np.nanmax(cn_lik_clonal) - np.nanmax(cn_lik_phi))
        p_val = stats.chisqprob(LLR,1) if not np.isnan(LLR) else 1
        if p_val < pval_cutoff:
            return cn_states[i][np.where(np.nanmax(cn_lik_clonal)==cn_lik_clonal)[0][0]] 
        else:
            return cn_states[i][np.where(np.nanmax(cn_lik_phi)==cn_lik_phi)[0][0]] 
    
    cn_ll_clonal = [calc_lik(cn_states[i],s[i],d[i],np.array(1),pi)[1] for i in range(len(cn_states))]
    cn_ll_phi = [calc_lik(cn_states[i],s[i],d[i],phi[i],pi)[1] for i in range(len(cn_states))]
    cn_ll_combined = zip(cn_ll_clonal, cn_ll_phi)
    
    most_likely_cn = [get_most_likely_cn(cn_states,cn_lik,i) for i,cn_lik in enumerate(cn_ll_combined)]
    cn_ll = [calc_lik(cn_states[i],s[i],d[i],phi[i],pi) for i in range(len(most_likely_cn))]
    most_likely_pv = [get_most_likely_pv(cn_lik) for i,cn_lik in enumerate(cn_ll)]

    return most_likely_cn, most_likely_pv

def cluster(sup,dep,cn_states,Nvar,sparams,cparams,phi_limit):
    '''
    clustering model using Dirichlet Process
    '''

    Ndp, iters, pi = cparams['clus_limit'], cparams['n_iter'], sparams['pi']
    Ndp_vals = [x for x in range(Ndp)]
    map_ = cparams['use_map']
    sens = 1.0 / ((sparams['pi']/float(sparams['ploidy']))*np.mean(dep))
    beta_a, beta_b = map(lambda x: float(eval(x)), cparams['beta'].split(','))
    print("Dirichlet concentration gamma values: alpha = %f, beta= %f" % (beta_a, beta_b))

    # to avoid out of bounds error
    Ndp = Nvar

    # pymc3 prefers lists vs. nparrays for some instances? can't say for sure
    sup = [int(si) for si in sup]
    cn_states = list(cn_states)

    with pm.Model() as model:

        alpha = pm.Gamma('alpha', alpha=beta_a, beta=beta_b, testval=beta_a/beta_b)
        h = pm.Beta('h', alpha=1, beta=alpha, shape=Ndp, testval=1/(1+alpha))
        phi_k = pm.Uniform('phi_k', lower=sens, upper=phi_limit, shape=Ndp)

        @as_op(itypes=[t.dvector], otypes=[t.dvector])
        def stick_breaking(h=h):
            value = [u * np.prod(1.0 - h[:i]) for i, u in enumerate(h)]
            value /= np.sum(value)
            if np.sum(value)!=1:
                value[len(value)-1] += 1-np.sum(value)#floaty roundy error
            return value

        stick_breaking.grad = lambda *x: x[0] #dummy gradient (otherwise fit function fails)
        p = stick_breaking(h)
        z = pm.Categorical('z', p=p, testval=0, shape=Nvar)

        @theano.compile.ops.as_op(itypes=[t.lvector, t.dvector], otypes=[t.dvector])
        def p_var(z=z, phi_k=phi_k):
            most_lik_cn_states, pvs = \
                    get_most_likely_cn_states(cn_states, sup, dep, phi_k[z], pi)
            return np.array(pvs)

        p_var.grad = lambda *x: [t.cast(x[0][0], dtype='float64'), x[0][1]]
        pv = p_var(z, phi_k)
        #cbinom = pm.Binomial('cbinom', dep, pv, shape=Nvar, observed=sup)

        bb_scale = 10
        bb_init = dep/bb_scale/0.001
        bb_beta = pm.Gamma('bb_beta', alpha=dep/bb_scale, beta=0.001, shape=Nvar, testval=bb_init)
        cbbinom = pm.BetaBinomial('cbbinom', alpha=-(bb_beta*pv)/(pv-1), beta=bb_beta, n=dep, observed=sup)

    with model:

        if map_:
            start = pm.find_MAP(fmin=optimize.fmin_cg)

        step1 = pm.Metropolis(vars=[alpha, h, p, phi_k])
        step2 = pm.ElemwiseCategorical(vars=[z], values=Ndp_vals)
        step3 = pm.Metropolis(vars=[pv, bb_beta, cbbinom])
        #step3 = pm.Metropolis(vars=[pv, cbinom])

        if map_:
            trace = pm.sample(iters, [step1, step2, step3], start)
        else:
            trace = pm.sample(iters, [step1, step2, step3])

    return trace, model
