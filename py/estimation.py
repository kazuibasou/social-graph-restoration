from collections import defaultdict

# Unbiased estimators based on re-weighted random walk

def size_estimator(samplinglist):
    # Estimating sizes of social networks via biased sampling, WWW, 2011.

    r = len(samplinglist)
    m = int(round(r*0.025)) # parameter
    phi = 0.0
    psi = 0.0
    for k in range(0,r-m):
        v_k = samplinglist[k].index
        d_k = len(samplinglist[k].nlist)
        for l in range(k+m, r):
            v_l = samplinglist[l].index
            d_l = len(samplinglist[l].nlist)
            if v_k == v_l:
                phi += 2
            psi += (float(d_k)/d_l) + (float(d_l)/d_k)

    return float(psi)/phi

def average_degree_estimator(samplinglist):
    # Walking in Facebook: A case study of unbiased sampling of osns, INFOCOM, 2010.

    est = 0.0
    for data in samplinglist:
        est += float(1)/len(data.nlist)

    return float(len(samplinglist))/est

def degree_distribution_estimator(samplinglist):
    # Walking in Facebook: A case study of unbiased sampling of osns, INFOCOM, 2010.

    est = defaultdict(float)
    x = 0
    for data in samplinglist:
        d = len(data.nlist)
        est[d] += float(1)/d
        x += float(1)/d

    for d in est:
        est[d] = float(est[d])/x

    return est

def jdd_estimator_induced_edges(samplinglist, est_n, est_aved):
    # 2.5 k-graphs: from sampling to generation, INFOCOM, 2013.
    # The original estimator does not correctly converge to the real value.
    # This is a modified unbiased estimator.

    phi = defaultdict(lambda: defaultdict(float))
    r = len(samplinglist)
    m = int(round(r * 0.025))

    for i in range(0,r-m):
        v = samplinglist[i].index
        k = len(samplinglist[i].nlist)
        for j in range(i+m,r):
            l = len(samplinglist[j].nlist)
            c = samplinglist[j].nlist.count(v)
            value = float(c)/(k*l)
            phi[k][l] += value
            phi[l][k] += value

    est = defaultdict(lambda: defaultdict(float))
    num_sample = (r-m)*(r-m+1)
    sum_d = est_n*est_aved
    for k in phi:
        for l in phi[k]:
            value = float(phi[k][l])/num_sample
            value *= sum_d
            est[k][l] = value

    return est

def jdd_estimator_traversed_edges(samplinglist):
    # 2.5 k-graphs: from sampling to generation, INFOCOM, 2013.
    # The original estimator does not correctly converge to a real value.
    # This is a modified unbiased estimator.

    est = defaultdict(lambda: defaultdict(float))
    r = len(samplinglist)

    for i in range(0, r-1):
        k = len(samplinglist[i].nlist)
        l = len(samplinglist[i+1].nlist)
        est[k][l] += 1
        est[l][k] += 1

    for k in est:
        for l in est[k]:
            est[k][l] = float(est[k][l])/(2*(r-1))

    return est

def JDD_estimator_hybrid(samplinglist, est_n, est_aved):
    # 2.5 k-graphs: from sampling to generation, INFOCOM, 2013.
    # The original estimator does not correctly converge to a real value.
    # This is a modified unbiased estimator.

    est_jdd_ie = jdd_estimator_induced_edges(samplinglist, est_n, est_aved)
    est_jdd_te = jdd_estimator_traversed_edges(samplinglist)

    est_jdd = defaultdict(lambda: defaultdict(float))

    for k in est_jdd_ie:
        for l in est_jdd_ie[k]:
            if (k + l) >= 2 * est_aved:
                est_jdd[k][l] = est_jdd_ie[k][l]
                est_jdd[l][k] = est_jdd_ie[l][k]

    for k in est_jdd_te:
        for l in est_jdd_te[k]:
            if (k+l) < 2*est_aved:
                est_jdd[k][l] = est_jdd_te[k][l]
                est_jdd[l][k] = est_jdd_te[l][k]

    return est_jdd

def degree_dependent_clustering_coefficient_estimator(samplinglist):
    # Estimating clustering coefficients and size of social networks via random walk, WWW, 2013.

    phi = defaultdict(float)
    psi = defaultdict(float)
    r = len(samplinglist)

    for i in range(0, r):
        d = len(samplinglist[i].nlist)
        psi[d] += float(1)/d

    for i in range(1, r-1):
        d = len(samplinglist[i].nlist)
        if d == 0 or d == 1:
            continue

        s = samplinglist[i-1].index
        v = samplinglist[i].index
        t = samplinglist[i+1].index

        if s != v and v != t and t != s:
            c = samplinglist[i+1].nlist.count(s)
            phi[d] += float(c)/(d-1)

    for d in phi:
        phi[d] = float(phi[d])/(r-2)

    for d in psi:
        psi[d] = float(psi[d])/r

    est_ddcc = defaultdict(float)
    for d in psi:
        est_ddcc[d] = float(phi[d])/psi[d]

    return est_ddcc