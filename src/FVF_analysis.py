import numpy as np

def quality_factor_result(model, val, vec, target_T=None, target_Q=None, target_theta_zeros=0, target_region='equator',
                          eq_cutoff=0.5, target_symmetric=True,oscillate=False,
                          wt_T=1., wt_Q=1., wt_zero=1., wt_region=1., wt_sym=1.):
    QFsq = 0.
    numQF = 0.
    if target_T:
        QFsq += (quality_T(model, val, target_T))**2
        numQF +=1
    if target_Q:
        QFsq += (quality_Q(model, val, target_Q))**2
        numQF += 1
    if target_theta_zeros:
        QFsq += (quality_theta_zeros(model, vec, target_theta_zeros, oscillate=oscillate))
        numQF +=1
    if target_region:
        QFsq += (quality_region(model, vec, target_region, eq_cutoff))**2
        numQF += 1
    if target_symmetric:
        QFsq += (quality_symmetric(model, vec))**2
        numQF += 1
    return (QFsq/numQF)**0.5


def quality_Q(model, val, target_Q):
    Q = np.abs(val.imag / (2 * val.real))
    return np.abs(Q-target_Q)/target_Q

def quality_T(model, val, target_T):
    T = (2 * np.pi / val.imag) * model.t_star / (24. * 3600. * 365.25)
    return np.abs(T-target_T)/target_T

def quality_theta_zeros(model, vec, target_theta_zeros, oscillate=False):
    zeros = get_theta_zero_crossings(model, vec, oscillate=oscillate)
    return np.abs(zeros-target_theta_zeros)/max(1, target_theta_zeros)

def quality_region(model, vec, target_region, eq_cutoff, var='ur'):
    var_out = model.get_variable(vec, var)
    split = eq_cutoff
    noneq_power = abs(np.concatenate((var_out[:, :int((model.Nl-1)*(0.5-split/2.))],
                                         var_out[:, int((model.Nl-1)*(0.5+split/2.)):]),
                                         axis=1)).sum()
    eq_power = abs(var_out[:, int((model.Nl-1)*(0.5-split/2.)):int((model.Nl-1)*(0.5+split/2.))]).sum()
    if target_region=='equator':
        return eq_power/(noneq_power+eq_power)
    else:
        return noneq_power/(noneq_power+eq_power)

def quality_symmetric(model, vec, var='ur'):
    var_out = model.get_variable(vec, var)
    north_power = np.abs(var_out[:, model.Nl//2:]).sum()
    south_power = np.abs(var_out[:, :model.Nl//2]).sum()
    return 1.-np.abs(north_power-south_power)/np.abs(north_power+south_power)

def apply_d2(model, vec):
    try:
        model.d2Mat
    except:
        model.make_d2Mat()
    return model.d2Mat.tocsr() * vec


def apply_dth(model, vec):
    try:
        model.dthMat
    except:
        model.make_dthMat()
    return model.dthMat.tocsr()*vec

def get_max_d2_norm(model, vec, var=None):
    if var:
        d2_var = model.get_variable(apply_d2(vec), var)
        var_out = model.get_variable(vec, var)
        return abs(d2_var).max()/abs(var_out).max()
    else:
        maxes = []
        d2_vec = apply_d2(model, vec)
        for var in model.model_variables:
            d2_var = model.get_variable(d2_vec, var)
            var_out = model.get_variable(vec, var)
            maxes.append(abs(d2_var).max()/abs(var_out).max())
    return max(maxes)

def get_max_dth_norm(model, vec, var=None):
    if var:
        dth_var = model.get_variable(apply_dth(vec), var)
        var_out = model.get_variable(vec, var)
        return abs(dth_var).max()/abs(var_out).max()
    else:
        maxes = []
        dth_vec = apply_dth(model, vec)
        for var in model.model_variables:
            dth_var = model.get_variable(dth_vec, var)
            var_out = model.get_variable(vec, var)
            maxes.append(abs(dth_var).max()/abs(var_out).max())
    return max(maxes)

def get_equator_power_excess(model, vec, var='ur', split=0.5):
    var_out = model.get_variable(vec, var)
    var_noneq_power = abs(np.concatenate((var_out[:, :int((model.Nl-1)*(0.5-split/2.))],
                                         var_out[:, int((model.Nl-1)*(0.5+split/2.)):]),
                                         axis=1)).sum()
    var_eq_power = abs(var_out[:, int((model.Nl-1)*(0.5-split/2.)):int((model.Nl-1)*(0.5+split/2.))]).sum()
    return var_eq_power-var_noneq_power

def shift_longitude(model, vec, phi):
    return vec*np.exp(1j*model.m_values[0]*phi)

def shift_vec_real(model, vec, var='ur'):
    ''' shift given vector's phase so that given variable (default ur) is
    dominantly real'''
    v = model.get_variable(vec, var)
    angs = np.angle(v) % np.pi
    abs_v = np.abs(v)
    avg_ang = np.average(angs, weights=abs_v) # shift phase angle
    # shift case to deal with vectors that are already dominantly real
    var_ang = np.average((angs - avg_ang)**2, weights=abs_v)
    if var_ang > 0.5:
        shift = np.exp(0.5j*np.pi)
        v_s = v*shift
        angs_s = np.angle(v_s) % np.pi
        avg_ang = np.average(angs_s, weights=abs_v)
        return vec*np.exp(-1j*(avg_ang-0.5*np.pi))
    else:
        return vec*np.exp(-1j*avg_ang)

def get_theta_zero_crossings(model, vec, var='uth', oscillate=False):
    z = model.get_variable(vec, var)
    if oscillate:
        z[:,::2] = -z[:,::2]
    ind = np.argmax(np.mean(np.abs(z),axis=1))
    signs = np.sign(z[ind,:])
    return np.where(signs[1:] != signs[:-1])[0]

def get_Q(model, val):
    return np.abs(val.imag/(2*val.real))

def filter_by_theta_zeros(model, vals, vecs, zeros_wanted, var='uth', verbose=False, oscillate=False):
    if type(zeros_wanted) is not list:
        zeros_wanted = list(zeros_wanted)
    if verbose:
        print('zeros wanted: {0}'.format(zeros_wanted))
    filtered_vals = []
    filtered_vecs = []
    for ind,(val, vec) in enumerate(zip(vals, vecs)):
        zc = get_theta_zero_crossings(model, vec, var=var, oscillate=oscillate)
        if verbose:
            print("{0}: val= {1}, zc = {2}".format(ind, val, zc))
        if len(zc)-1 in zeros_wanted:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_dth(model, vals, vecs, max_dth):
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        if get_max_dth_norm(model, vec) < max_dth:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_d2(model, vals, vecs, max_d2):
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        if get_max_d2_norm(model, vec) < max_d2:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_Q(model, vals, vecs, min_Q):
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        if get_Q(model, val) > min_Q:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_equator_power(model, vals, vecs, equator_fraction=0.5, var='ur'):
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        if get_equator_power_excess(model, vec, var=var,
                                    split=equator_fraction) > 0.:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_period(model, vals, vecs, minT, maxT):
    '''
    filter results by wave period in years 
    '''
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        Period = (2*np.pi/val.imag)*model.t_star/(24.*3600.*365.25)
        if Period > minT and Period < maxT:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def get_period(model, val):
    return 2*np.pi/val*model.t_star/(24.*3600.*365.25)