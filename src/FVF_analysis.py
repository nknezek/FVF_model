import numpy as np


def filter_by_misfit(model, vals, vecs, num_to_keep, target_T=None, target_Q=None, target_r_order=0, target_th_order=0,
                     target_region='equator', th_ord_var ='uph', r_ord_var ='uph', eq_var='uph',
                     eq_cutoff=0.5, target_symmetric=True,
                     wt_T=1., wt_Q=1., wt_th_order=1., wt_r_order=1., wt_region=1., wt_sym=1., wt_r_sm=1., wt_th_sm=1.):
    mf = np.zeros(len(vals))
    for i,(val, vec) in enumerate(zip(vals, vecs)):
        mf[i] = misfit_result(model, val, vec, target_T=target_T, target_Q=target_Q, target_r_order=target_r_order,
                                   target_th_order=target_th_order, target_region=target_region,
                                   eq_cutoff=eq_cutoff, target_symmetric=target_symmetric,
                                th_ord_var=th_ord_var, r_ord_var=r_ord_var, eq_var=eq_var,
                                   wt_T=wt_T, wt_Q=wt_Q, wt_r_order=wt_r_order, wt_th_order=wt_th_order,
                                    wt_region=wt_region, wt_sym=wt_sym, wt_r_sm=wt_r_sm, wt_th_sm=wt_th_sm)
    sorted_ind = np.argsort(mf)
    fvals = []
    fvecs = []
    for i in range(num_to_keep):
        fvals.append(vals[sorted_ind[i]])
        fvecs.append(vecs[sorted_ind[i]])
    return fvals, fvecs

def filter_by_rth_zeros(model, vals, vecs):
    max_th_zeros = model.Nl//2+1
    max_r_zeros = model.Nk//2+1
    fvals = []
    fvecs = []
    for val, vec in zip(vals, vecs):
        if get_theta_zero_crossings(model, vec) <= max_th_zeros:
            if get_r_zero_crossings(model, vec) <= max_r_zeros:
                fvals.append(val)
                fvecs.append(vec)
    return fvals, fvecs


def misfit_result(model, val, vec, target_T=None, target_Q=None, target_r_order=0, target_th_order=0, target_region='equator',
                          eq_cutoff=0.5, target_symmetric=True,
                        th_ord_var ='uph', r_ord_var ='uph', eq_var='uph',
                          wt_T=1., wt_Q=1., wt_r_order=1., wt_th_order=1., wt_region=1., wt_sym=1., wt_r_sm=1, wt_th_sm=1.):
    mfsq = 0.
    nummf = 0.
    if target_T is not None:
        mf_T = misfit_T(model, val, target_T)
        # print("\n T mf = {0}".format(mf_T))
        mfsq += (mf_T*wt_T)**2
        nummf +=1
    if target_Q is not None:
        mf_Q = misfit_Q(model, val, target_Q)
        # print("Q mf = {0}".format(mf_Q))
        mfsq += (mf_Q*wt_Q)**2
        nummf += 1
    if target_th_order is not None:
        mf_th_ord = misfit_th_order(model, vec, target_th_order,  var=th_ord_var)
        # print("th ord mf = {0}".format(mf_th_ord))
        mfsq += (mf_th_ord*wt_th_order)**2
        nummf +=1
    if target_r_order is not None:
        mf_r_ord = misfit_r_order(model, vec, target_r_order,  var=r_ord_var)
        # print("r ord mf = {0}".format(mf_r_ord))
        mfsq += (mf_r_ord*wt_r_order)**2
        nummf +=1
    if target_region is not None:
        mf_region = misfit_region(model, vec, target_region, eq_cutoff, var=eq_var)
        # print("region mf = {0}".format(mf_region))
        mfsq += (mf_region*wt_region)**2
        nummf += 1
    if target_symmetric is True:
        mf_sym = misfit_symmetric(model, vec, var=eq_var)
        # print("sym mf = {0}".format(mf_sym))
        mfsq += (mf_sym*wt_sym)**2
        nummf += 1
    if wt_th_sm > 0.:
        mf_th_sm = misfit_smoothness_th(model, vec, var=th_ord_var)
        # print("th sm mf = {0}".format(mf_th_sm))
        mfsq += (mf_th_sm*wt_th_sm)**2
        nummf +=1
    if wt_r_sm > 0.:
        mf_r_sm = misfit_smoothness_r(model, vec, var=r_ord_var)
        # print("r sm mf = {0}".format(mf_r_sm))
        mfsq += (mf_r_sm*wt_r_sm)**2
        nummf +=1
    return (mfsq/nummf)**0.5


def misfit_Q(model, val, target_Q):
    Q = np.abs(val.imag / (2 * val.real))
    return np.abs(Q-target_Q)/max(1.,target_Q)

def misfit_T(model, val, target_T):
    T = (2 * np.pi / val.imag) * model.t_star / (24. * 3600. * 365.25)
    return np.abs(T-target_T)/target_T

def misfit_th_order(model, vec, target_th_order,  var='uph'):
    zeros = get_theta_zero_crossings(model, vec,  var=var)
    return np.abs(zeros-target_th_order)/max(1, target_th_order)/4

def misfit_r_order(model, vec, target_r_order,  var='uph'):
    zeros = get_r_zero_crossings(model, vec,  var=var)
    return np.abs(zeros-target_r_order)/max(1, target_r_order)/4

def misfit_region(model, vec, target_region, eq_cutoff, var='uph'):
    var_out = model.get_variable(vec, var)
    split = eq_cutoff
    noneq_power = abs(np.concatenate((var_out[:, :int((model.Nl-1)*(0.5-split/2.))],
                                         var_out[:, int((model.Nl-1)*(0.5+split/2.)):]),
                                         axis=1)).sum()
    eq_power = abs(var_out[:, int((model.Nl-1)*(0.5-split/2.)):int((model.Nl-1)*(0.5+split/2.))]).sum()
    if target_region=='equator':
        return noneq_power/(noneq_power+eq_power)
    else:
        return eq_power/(noneq_power+eq_power)

def misfit_symmetric(model, vec, var='uph'):
    var_out = model.get_variable(vec, var)
    north_power = np.mean(np.abs(var_out[:, model.Nl//2:]))
    south_power = np.mean(np.abs(var_out[:, :model.Nl//2]))
    return np.abs((north_power-south_power)/(north_power+south_power))*10

def misfit_smoothness_th(model, vec, var='uph'):
    y = (model.get_variable(vec, var)).real
    return np.mean(np.abs(y[:,2:]+y[:,:-2]-2*y[:,1:-1]))/np.mean(np.abs(y))

def misfit_smoothness_r(model, vec, var='uph'):
    y = (model.get_variable(vec, var)).real
    return np.mean(np.abs(y[2:,:]+y[:-2,:]-2*y[1:-1,:]))/np.mean(np.abs(y))


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

def get_theta_zero_crossings(model, vec, var='uth', cutoff=0.075):
    z = model.get_variable(vec, var)
    ind = np.argmax(np.mean(np.abs(z),axis=1))
    zi = z[ind,1:-1]
    zi[np.abs(zi) < cutoff * np.max(np.abs(zi))] = 0.
    signs = np.sign(zi)
    stripped_signs = signs[np.nonzero(signs)]
    zeros = np.where(stripped_signs[1:] != stripped_signs[:-1])[0]
    return len(zeros)

def get_r_zero_crossings(model, vec, var='uph', cutoff=0.075):
    z = model.get_variable(vec, var)
    ind = np.argmax(np.mean(np.abs(z),axis=0))
    zi = z[1:-1,ind]
    zi[np.abs(zi) < cutoff * np.max(np.abs(zi))] = 0.
    signs = np.sign(zi)
    stripped_signs = signs[np.nonzero(signs)]
    zeros = np.where(stripped_signs[1:] != stripped_signs[:-1])[0]
    return len(zeros)

def get_Q(model, val):
    return np.abs(val.imag/(2*val.real))

def filter_by_theta_zeros(model, vals, vecs, zeros_wanted, var='uth', verbose=False):
    if type(zeros_wanted) is not list:
        zeros_wanted = list(zeros_wanted)
    if verbose:
        print('zeros wanted: {0}'.format(zeros_wanted))
    filtered_vals = []
    filtered_vecs = []
    for ind,(val, vec) in enumerate(zip(vals, vecs)):
        zc = get_theta_zero_crossings(model, vec, var=var)
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