import numpy as np
import pandas as pd
import pymc as pm
import math
import ipdb
import collections

from scipy import stats
from operator import methodcaller

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
        c1 = map(float,cn[0].split(','))
        c2 = map(float,cn[1].split(','))
        major1, minor1, total1, frac1 = c1[0], c1[1], c1[0]+c1[1], c1[2]
        major2, minor2, total2, frac2 = c2[0], c2[1], c2[0]+c2[1], c2[2]

        # generate copynumbers for each subclone being the potential variant pop 
        combos = add_copynumber_combos(combos, major1, minor1, total2, frac1)
        combos = add_copynumber_combos(combos, major2, minor2, total1, frac2)
    else:
        c = map(float,cn[0].split(','))
        major, minor, frac = c[0], c[1], c[2]
        combos = add_copynumber_combos(combos, major, minor, major + minor, frac)

    return filter_cns(combos)

def get_sv_allele_combos(sv):
    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
    combos_bp1 = get_allele_combos(cn_tmp[0])
    combos_bp2 = get_allele_combos(cn_tmp[1])

    return tuple([combos_bp1,combos_bp2])

def fit_and_sample(model, iters, burn, thin, use_map):
    if use_map:
        map_ = pm.MAP( model )
        map_.fit(method = 'fmin_cg')

    mcmc = pm.MCMC( model )
    #burn-in and thinning now done in post processing
    #mcmc.sample( iters, burn=burn, thin=thin )
    mcmc.sample( iters )

    if use_map:
        return mcmc, map_
    else:
        return mcmc, None

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
    return [map(float,cn) for cn in map(methodcaller('split',','),cn_str)]

def calc_lik_with_clonal(combo, si, di, phi_i, pi):
    pvs = np.array([get_pv(phi_i, c[0], c[1], c[2], c[3], pi) for c in combo])
    lls = np.array([pm.binomial_like(si, di, pvs[i]) for i,c in enumerate(combo)])-0.00000001
    # also calculate with clonal phi
    pvs_cl = np.array([get_pv(np.array(1), c[0], c[1], c[2], c[3], pi) for c in combo])
    lls_cl = np.array([pm.binomial_like(si, di, pvs[i]) for i,c in enumerate(combo)])-0.00000001
    #lls currently uses precision fudge factor to get
    #around 0 probability errors when pv = 1
    #TODO: investigate look this bug more
    return np.array([[pvs, lls], [pvs_cl, lls_cl]])

def calc_lik(combo, si, di, phi_i, pi):
    pvs = np.array([get_pv(phi_i, c[0], c[1], c[2], c[3], pi) for c in combo])
    lls = np.array([pm.binomial_like(si, di, pvs[i]) for i,c in enumerate(combo)])-0.00000001
    return np.array([pvs, lls])

def get_probs(var_states,s,d,phi,pi):
    llik = calc_lik(var_states,s,d,phi,pi)[1]
    probs = get_probs_from_llik(llik)
    probs = ','.join(map(lambda x: str(round(x,4)),probs))
    return probs

def get_probs_from_llik(cn_lik):
    probs = np.array([1.])
    if len(cn_lik)>1:        
        probs = map(math.exp,cn_lik)
        probs = np.array(probs)/sum(probs)
    return probs

def index_of_max(lik_list):
    result = collections.defaultdict(list)
    for val, idx in enumerate(lik_list):
        result[idx] = val
    return result[np.nanmax(lik_list)]

def get_most_likely_pv(cn_lik):
    if np.all(np.isnan(cn_lik[0])):
        return 0.0000001
    elif len(cn_lik[0]) > 0:
        return cn_lik[0][index_of_max(cn_lik[1])]
        #return cn_lik[0][np.where(np.nanmax(cn_lik[1])==cn_lik[1])[0][0]]
    else:
        return 0.0000001

def get_most_likely_cn(cn_states, cn_lik, i, pval_cutoff):
    '''
    use the most likely phi state, unless p < cutoff when compared to the
    most likely clonal (phi=1) case (log likelihood ratio test)
    - in this case, pick the most CN state with the highest clonal likelihood
    '''
    cn_lik_clonal, cn_lik_phi = cn_lik

    #reinstitute hack - uncomment below
    #cn_lik_phi = cn_lik_clonal
    if len(cn_lik_clonal)==0:
        return [float('nan'), float('nan'), float('nan')]

    # log likelihood ratio test; null hypothesis = likelihood under phi
    LLR   = 2 * (np.nanmax(cn_lik_clonal) - np.nanmax(cn_lik_phi))
    p_val = stats.chisqprob(LLR, 1) if not np.isnan(LLR) else 1

    if np.all(np.isnan(cn_lik_phi)):
        return cn_states[i][0]
    elif p_val < pval_cutoff:
        return cn_states[i][index_of_max(cn_lik_clonal)]
    else:
        return cn_states[i][index_of_max(cn_lik_phi)]

def get_most_likely_cn_states(cn_states, s, d, phi, pi, pval_cutoff):
    '''
    Obtain the copy-number states which maximise the binomial likelihood
    of observing the supporting read depths at each variant location
    '''
    #cn_ll_clonal = [calc_lik(cn_states[i],s[i],d[i],np.array(1),pi)[1] for i in range(len(cn_states))]
    cn_ll_combined = [calc_lik_with_clonal(cn_states[i],s[i],d[i],phi[i],pi)[1] for i in range(len(cn_states))]
    #cn_ll_combined = zip(cn_ll_clonal,cn_ll_phi)

    most_likely_cn = [get_most_likely_cn(cn_states,cn_lik,i,pval_cutoff) for i, cn_lik in enumerate(cn_ll_combined)]
    cn_ll = [calc_lik(cn_states[i],s[i],d[i],phi[i],pi) for i in range(len(most_likely_cn))]
    most_likely_pv = [get_most_likely_pv(cn_lik) for cn_lik in cn_ll]

    return most_likely_cn, most_likely_pv

def cluster(sup,dep,cn_states,Nvar,sparams,cparams,phi_limit):
    '''
    clustering model using Dirichlet Process
    '''
    Ndp = cparams['clus_limit']
    purity, ploidy = sparams['pi'], sparams['ploidy']
    sens = 1.0 / ((purity/ float(ploidy)) * np.mean(dep))

    beta_a, beta_b = map(lambda x: float(eval(x)), cparams['beta'].split(','))
    beta_init = float(beta_a) / beta_b
    print("Dirichlet concentration gamma values: alpha = %f; beta= %f; init = %f" % (beta_a, beta_b, beta_init))

    alpha = pm.Gamma('alpha', beta_a, beta_b, value = beta_init)
    pval_cutoff = cparams['clonal_cnv_pval']

    h = pm.Beta('h', alpha=1, beta=alpha, size=Ndp)
    @pm.deterministic
    def p(h=h):
        value = [u*np.prod(1.0-h[:i]) for i,u in enumerate(h)]
        value /= np.sum(value)
        return value

    z = pm.Categorical('z', p = p, size = Nvar, value = np.zeros(Nvar))
    phi_init = max(sens, min(np.mean(sup) / np.mean(dep) / purity * 2, 1))
    phi_k = pm.Uniform('phi_k', lower = sens, upper = phi_limit, size = Ndp)#, value=[phi_init]*Ndp)
    #print('phi lower limit: %f; phi upper limit: %f; phi init: %f' % (sens, phi_limit, phi_init))
    
    @pm.deterministic
    def p_var(z=z,phi_k=phi_k):
        if np.any(z < 0):
            z = [0 for x in z]
            # ^ some fmin optimization methods initialise this array with -ve numbers
        most_lik_cn_states, pvs = \
                get_most_likely_cn_states(cn_states, sup, dep, phi_k[z], purity, pval_cutoff)
        return pvs

    cbinom = pm.Binomial('cbinom', dep, p_var, observed=True, value=sup)
    model = pm.Model([alpha, h, p, phi_k, z, p_var, cbinom])
    mcmc, map_ = fit_and_sample(model,cparams['n_iter'],cparams['burn'],cparams['thin'],cparams['use_map'])

    return mcmc, map_
