import numpy as np
import pandas as pd
import pymc3 as pm
import pymc as pm2

import math
import ipdb
import collections
import theano
import theano.tensor as t

from theano.compile.ops import as_op
from scipy import stats, optimize
from sklearn.cluster import KMeans
from operator import methodcaller

def add_copynumber_combos(combos, var_maj, var_min, ref_cn, frac, cparams):
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
    elif cparams['restrict_cnss']:
        combos.append([ref_cn, var_total, 1. / var_total, frac])
        combos.append([ref_cn, var_total, float(var_maj) / var_total, frac])
        combos.append([ref_cn, var_total, float(var_min) / var_total, frac])
    else:
        for i in range(1, int(var_maj+1)):
            mu_v = 1.0*i / var_total
            combos.append([ref_cn, var_total, mu_v, frac])

    return combos

def get_allele_combos(cn, cparams):
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
        combos = add_copynumber_combos(combos, major1, minor1, total2, frac1, cparams)
        combos = add_copynumber_combos(combos, major2, minor2, total1, frac2, cparams)
    else:
        c = [float(c) for c in cn[0].split(',')]
        major, minor, frac = c[0], c[1], c[2]
        combos = add_copynumber_combos(combos, major, minor, major + minor, frac, cparams)

    return filter_cns(combos)

def get_sv_allele_combos(sv, cparams):
    cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
    combos_bp1 = get_allele_combos(cn_tmp[0], cparams)
    combos_bp2 = get_allele_combos(cn_tmp[1], cparams)

    return tuple([combos_bp1,combos_bp2])

def fit_and_sample(model, iters, burn, thin, use_map):
    if use_map:
        map_ = pm.MAP( model )
        map_.fit(method = 'fmin_cg')

    mcmc = pm.MCMC( model )
    #burn-in and thinning now done in post processing
    mcmc.sample( iters, burn=burn, thin=thin )

    if use_map:
        return mcmc, map_
    else:
        return mcmc, None

def get_pv(phi, combo, pi, ni):
    cn_r, cn_v, mu, cn_f = combo

    pn = (1.0 - pi) * ni        #proportion of normal reads coming from normal cells
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

def calc_lik_with_clonal(combo, si, di, phi_i, pi, ni):

    # calculate with given phi
    # (lls currently uses precision fudge factor to get around 0 probability errors when pv = 1)
    pvs = np.array([get_pv(phi_i, c, pi, ni) for c in combo])
    lls = np.array([pm2.binomial_like(si, di, pvs[i]) for i,c in enumerate(combo)])-0.00000001

    # calculate with clonal phi
    pvs_cl = np.array([get_pv(np.array(1), c, pi, ni) for c in combo])
    lls_cl = np.array([pm2.binomial_like(si, di, pvs_cl[i]) for i,c in enumerate(combo)])-0.00000001

    return np.array([[pvs, lls], [pvs_cl, lls_cl]])

def calc_lik(combo, si, di, phi_i, pi, ni):
    pvs = np.array([get_pv(phi_i, c, pi, ni) for c in combo])
    lls = np.array([pm2.binomial_like(si, di, pvs[i]) for i,c in enumerate(combo)])-0.00000001
    return np.array([pvs, lls])

def get_probs(var_states,s,d,phi,pi,norm):
    llik = calc_lik(var_states,s,d,phi,pi,norm)[1]
    probs = get_probs_from_llik(llik)
    probs = ','.join(map(lambda x: str(round(x,4)),probs))
    return probs

def get_probs_from_llik(cn_lik):
    probs = np.array([1.])
    if len(cn_lik) > 1:
        probs = [math.exp(x) for x in cn_lik]
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
    else:
        return 0.0000001

def get_most_likely_cn(combo, cn_lik, pval_cutoff):
    '''
    use the most likely phi state, unless p < cutoff when compared to the
    most likely clonal (phi=1) case (log likelihood ratio test)
    - in this case, pick the most CN state with the highest clonal likelihood
    '''
    cn_lik_phi, cn_lik_clonal = cn_lik
    ll_phi, ll_clonal = cn_lik_phi[1], cn_lik_clonal[1]

    empty_result = [float('nan'), float('nan'), float('nan'), float('nan')]
    if len(combo) == 0:
        return empty_result
    elif len(combo) == 1:
        return combo[0]
    elif np.all(np.isnan(ll_phi)) and np.all(np.isnan(ll_clonal)):
        return empty_result
    elif np.all(np.isnan(ll_phi)):
        return combo[index_of_max(ll_clonal)]
    elif np.all(ll_phi == ll_clonal) or pval_cutoff == 0:
        return combo[index_of_max(ll_phi)]

    # log likelihood ratio test; null hypothesis = likelihood under phi
    # use clonal if best clonal solution significantly better than worst phi solution
    #LLR   = 2 * (np.nanmax(ll_clonal) - np.nanmax(ll_phi))
    LLR   = 2 * (np.nanmax(ll_clonal) - np.nanmin(ll_phi))
    p_val = stats.chisqprob(LLR, 1) if not np.isnan(LLR) else 1

    if p_val < pval_cutoff:
        return combo[index_of_max(ll_clonal)]
    else:
        return combo[index_of_max(ll_phi)]

def get_most_likely_cn_states(cn_states, s, d, phi, pi, pval_cutoff, norm):
    '''
    Obtain the copy-number states which maximise the binomial likelihood
    of observing the supporting read depths at each variant location
    '''
    cn_ll_combined = [calc_lik_with_clonal(cn_states[i],s[i],d[i],phi[i],pi,norm[i]) for i in range(len(cn_states))]
    most_likely_cn = [get_most_likely_cn(cn_states[i],cn_lik,pval_cutoff) for i, cn_lik in enumerate(cn_ll_combined)]

    cn_ll = [calc_lik(cn_states[i],s[i],d[i],phi[i],pi,norm[i]) for i in range(len(most_likely_cn))]
    most_likely_pv = [get_most_likely_pv(cn_lik) for cn_lik in cn_ll]

    return most_likely_cn, most_likely_pv


def get_initialisation(nclus_init, Ndp, sparams, sup, dep, norm, cn_states, sens, phi_limit, pval_cutoff):

    purity, ploidy, mean_cov = sparams['pi'], sparams['ploidy'], sparams['mean_cov']
    mlcn, mlpv = get_most_likely_cn_states(cn_states, sup, dep, np.ones(len(sup)), purity, pval_cutoff, norm)
    data = (sup / dep) * (1 / np.array(mlpv))
    data = np.array([d if d < phi_limit else phi_limit for d in data])
    data = np.array([d if d > sens else sens for d in data])

    if nclus_init < 1:
        # derive sensible N cluster initialisation, based on mean chromosomal copy depth
        nrpcc = mean_cov * (purity / (purity * ploidy + (1 - purity) * 2))
        nclus_init = int(round(nrpcc / 5)) # min resolution is about 5 reads between cluster means

    kme = KMeans(nclus_init)
    kme.fit(data.reshape(-1,1))

    print('Cluster initialisation set to %d' % nclus_init)
    z_init = kme.labels_
    phi_init = [c[0] if c[0] < phi_limit else phi_limit for c in kme.cluster_centers_]

    # fill the rest of the phi slots (up to max clusters) randomly
    if (Ndp - len(phi_init)) > 0:
        phi_fill = np.random.rand(Ndp - len(phi_init)) * phi_limit
        phi_fill = np.array([sens if x < sens else x for x in phi_fill])
        phi_init = np.concatenate([phi_init, phi_fill])

#    # for testing initalisation
#    import colorsys
#    import matplotlib.pyplot as plt
#    def gen_new_colours(N):
#        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
#        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
#        return RGB_tuples
#    RGB_tuples = np.array(gen_new_colours(Ndp))
#
#    for lab in np.unique(kme.labels_):
#        dat = data[np.where(kme.labels_==lab)[0]]
#        plt.hist(dat,bins=np.array(range(0,200,2))/100.,alpha=0.75,color=RGB_tuples[lab])
#    plt.savefig('test')

    return z_init, phi_init

def cluster(sup,dep,cn_states,Nvar,sparams,cparams,phi_limit,norm,recluster=False):
    '''
    clustering model using Dirichlet Process
    '''
    Ndp = cparams['clus_limit'] if not recluster else 1
    n_iter = cparams['n_iter'] if not recluster else cparams['merge_iter']
    burn = cparams['burn'] if not recluster else cparams['merge_burn']
    thin, use_map = cparams['thin'], cparams['use_map']
    use_map = False if recluster else use_map
    nclus_init = cparams['nclus_init']

    purity, ploidy = sparams['pi'], sparams['ploidy']
    fixed_alpha, gamma_a, gamma_b = cparams['fixed_alpha'], cparams['alpha'], cparams['beta']
    sens = 1.0 / ((purity/ float(ploidy)) * np.mean(dep))
    pval_cutoff = cparams['clonal_cnv_pval']
    print('phi lower limit: %f; phi upper limit: %f' % (sens, phi_limit))

    alpha, beta = cparams['alpha'], cparams['beta']
    print("Dirichlet concentration gamma values: alpha = %f, beta= %f" % (alpha, beta))

    a_s, b_s = map(lambda x: float(eval(x)), cparams['precision'].split(','))
    print("Precision values: a_s = %f, a_b = %f" % (a_s, b_s))

    # to avoid out of bounds error
    Ndp = Nvar

    # pymc3 prefers lists vs. nparrays? can't say for sure
    sup = [int(si) for si in sup]
    cn_states = list(cn_states)

#    z_init = np.zeros(Nvar, dtype=np.int)
#    phi_init = np.random.rand(Ndp) * phi_limit
#
#    # use smart initialisation if nclus_init specified
#    if not nclus_init.lower() in ("no", "false", "f"):
#        try:
#            nclus_init = nclus_init if not recluster else 1
#            nclus_init = int(nclus_init)
#            nclus_init = Ndp if nclus_init > Ndp else nclus_init
#            if nclus_init == 1:
#                phi_init[0] = 1.
#            else:
#                z_init, phi_init = get_initialisation(nclus_init, Ndp, sparams, sup, dep, norm,
#                                                      cn_states, sens, phi_limit, pval_cutoff)
#        except ValueError:
#            pass

    model = pm.Model()
    with model:
        alpha = pm.Gamma('alpha', alpha=alpha, beta=beta, testval=alpha/beta)
        h = pm.Beta('h', alpha=1, beta=alpha, shape=Ndp, testval=1/(1+alpha))
        phi_k = pm.Uniform('phi_k', lower=sens, upper=phi_limit, shape=Ndp)

        @as_op(itypes=[t.dvector], otypes=[t.dvector])
        def stick_breaking(h=h):
            value = [u * np.prod(1.0 - h[:i]) for i, u in enumerate(h)]
            value[-1] = 1-sum(value[:-1])
            return np.array(value)

        stick_breaking.grad = lambda *x: x[0] #dummy gradient (otherwise fit function fails)
        p = stick_breaking(h)
        z = pm.Categorical('z', p=p, testval=0, shape=Nvar)

        @theano.compile.ops.as_op(itypes=[t.lvector, t.dvector], otypes=[t.dvector])
        def p_var(z=z, phi_k=phi_k):
            most_lik_cn_states, pvs = \
                    get_most_likely_cn_states(cn_states, sup, dep, phi_k[z], purity, pval_cutoff, norm)
            return np.array(pvs)

        p_var.grad = lambda *x: [t.cast(x[0][0], dtype='float64'), x[0][1]]
        pv = p_var(z, phi_k)
        #cbinom = pm.Binomial('cbinom', dep, pv, shape=Nvar, observed=sup)

        bb_scale = a_s
        bb_init  = dep/bb_scale/b_s
        bb_beta  = pm.Gamma('bb_beta', alpha=dep/bb_scale, beta=b_s, shape=Nvar, testval=bb_init)
        cbbinom  = pm.BetaBinomial('cbbinom', alpha=-(bb_beta*pv)/(pv-1), beta=bb_beta, n=dep, observed=sup)

        if use_map:
            start = pm.find_MAP(fmin=optimize.fmin_cg)

        Ndp_vals = [x for x in range(Ndp)]

        step1 = pm.Metropolis(vars=[alpha, h, p, phi_k])
        #step2 = pm.CategoricalGibbsMetropolis(vars = [z], values=Ndp_vals)
        step2 = pm.ElemwiseCategorical(vars=[z], values=Ndp_vals)
        step3 = pm.Metropolis(vars=[pv, bb_beta, cbbinom])
        #step3 = pm.Metropolis(vars=[pv, cbinom])

        if use_map:
            trace = pm.sample(n_iter, [step1, step2, step3], start)
        else:
            trace = pm.sample(n_iter, [step1, step2, step3])

    return trace, model
