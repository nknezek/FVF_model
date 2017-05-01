import numpy as np

def filter_results(model, vals, vecs, filter_dict):
    fvals = vals
    fvecs = vecs
    for k in filter_dict.keys():
        fvals, fvecs = filter_by_type(k, filter_dict[k], model, fvals, fvecs)
    return fvals, fvecs

def filter_by_type(type_of_filter, parameters, model, vals, vecs):
    if type_of_filter is 'Q':
        return filter_by_Q(vals, vecs, parameters)
    if type_of_filter is 'order_r':
        return filter_by_order_r(model, vals, vecs, parameters)
    if type_of_filter is 'order_th':
        return filter_by_order_th(model, vals, vecs, parameters)

def filter_by_order_th(model, vals, vecs, parameters):
    filtered_vals = []
    filtered_vecs = []
    for ind,(val, vec) in enumerate(zip(vals, vecs)):
        zc = get_order_th(model, vec, var=parameters['var'])
        keep = True
        if 'maximum' in parameters.keys():
            if zc > parameters['maximum']:
                keep = False
        if 'minimum' in parameters.keys():
            if zc < parameters['minimum']:
                keep = False
        if keep:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_order_r(model, vals, vecs, parameters):
    filtered_vals = []
    filtered_vecs = []
    for ind,(val, vec) in enumerate(zip(vals, vecs)):
        zc = get_order_r(model, vec, var=parameters['var'])
        keep = True
        if 'maximum' in parameters.keys():
            if zc > parameters['maximum']:
                keep = False
        if 'minimum' in parameters.keys():
            if zc < parameters['minimum']:
                keep = False
        if keep:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def filter_by_Q(vals, vecs, parameters):
    filtered_vals = []
    filtered_vecs = []
    for ind, (val, vec) in enumerate(zip(vals, vecs)):
        Q = get_Q(val)
        keep = True
        if 'maximum' in parameters.keys():
            if Q > parameters['maximum']:
                keep = False
        if 'minimum' in parameters.keys():
            if Q < parameters['minimum']:
                keep = False
        if keep:
            filtered_vals.append(val)
            filtered_vecs.append(vec)
    return filtered_vals, filtered_vecs

def sort_by_total_misfit(model, vals, vecs, misfit_dict):
    mf = np.zeros(len(vals))
    for i,(val, vec) in enumerate(zip(vals, vecs)):
        mf[i] = misfit_rms(model, val, vec, misfit_dict)
    sorted_ind = np.argsort(mf)
    svals = []
    svecs = []
    for i in range(len(vals)):
        svals.append(vals[sorted_ind[i]])
        svecs.append(vecs[sorted_ind[i]])
    return svals, svecs

def misfit_rms(model, val, vec, misfit_dict):
    mfsq = 0.
    for k in misfit_dict.keys():
        mfsq += (misfit_by_type(k, misfit_dict[k], model, val, vec)*misfit_dict[k]['weight'])**2
    return (mfsq/len(misfit_dict))**0.5

def misfit_by_type(type_of_misfit, parameters, model, val, vec):
    if type_of_misfit is 'Q':
        return misfit_Q(val, parameters)
    elif type_of_misfit is 'T':
        return misfit_T(model, val, parameters)
    elif type_of_misfit is 'order_r':
        return misfit_order_r(model, vec, parameters)
    elif type_of_misfit is 'order_th':
        return misfit_order_th(model, vec, parameters)
    elif type_of_misfit is 'region':
        return misfit_region(model, vec, parameters)
    elif type_of_misfit is 'symmetry':
        return misfit_symmetry(model, vec, parameters)
    elif type_of_misfit is 'smoothness_r':
        return misfit_smoothness_r(model, vec, parameters)
    elif type_of_misfit is 'smoothness_th':
        return misfit_smoothness_th(model, vec, parameters)
    elif type_of_misfit is 'power_in_layers':
        return misfit_power_in_layers(model, vec, parameters)

def misfit_Q(val, parameters):
    '''
    return misfit of how far away the quality factor is from the target quality factor

    :param val: eigenvalue
    :param parameters: dict containing 'target', the desired Q value
    :return: misfit (range 0 to 1)
    '''
    Q = np.abs(val.imag / (2 * val.real))
    return np.abs(Q - parameters['target']) / max(Q, parameters['target'])

def misfit_T(model, val,parameters):
    T = (2 * np.pi / np.abs(val.imag)) * model.t_star / (24. * 3600. * 365.25)
    return np.abs(T-parameters['target'])/max(parameters['target'], T)

def misfit_order_th(model, vec, parameters):
    zeros = get_order_th(model, vec,  var=parameters['var'])
    return saturation(np.abs(zeros-parameters['target']), c=2)

def misfit_order_r(model, vec, parameters):
    zeros = get_order_r(model, vec,  var=parameters['var'])
    return saturation(np.abs(zeros-parameters['target']), c=2)

def misfit_region(model, vec, parameters):
    var_out = model.get_variable(vec, parameters['var'])
    split = parameters['split']
    noneq_power = abs(np.concatenate((var_out[:, :int((model.Nl-1)*(0.5-split/2.))],
                                         var_out[:, int((model.Nl-1)*(0.5+split/2.)):]),
                                         axis=1)).sum()
    eq_power = abs(var_out[:, int((model.Nl-1)*(0.5-split/2.)):int((model.Nl-1)*(0.5+split/2.))]).sum()
    if parameters['target']=='equator':
        return noneq_power/(noneq_power+eq_power)
    else:
        return eq_power/(noneq_power+eq_power)

def misfit_symmetry(model, vec, parameters):
    var_out = model.get_variable(vec, parameters['var'])
    north_power = np.mean(np.abs(var_out[:, model.Nl//2:]))
    south_power = np.mean(np.abs(var_out[:, :model.Nl//2]))
    if parameters['target'] is 'symmetric':
        return np.abs((north_power - south_power)/(north_power+south_power))
    elif parameters['target'] is 'asymmetric':
        return 1.-np.abs((north_power + south_power) / (north_power + south_power))

def misfit_smoothness_r(model, vec, parameters):
    y = (model.get_variable(vec, parameters['var'])).real
    return saturation(np.mean(np.abs(y[2:,:]+y[:-2,:]-2*y[1:-1,:]))/np.mean(np.abs(y)), c=1)

def misfit_smoothness_th(model, vec, parameters):
    y = (model.get_variable(vec, parameters['var'])).real
    return saturation(np.mean(np.abs(y[:,2:]+y[:,:-2]-2*y[:,1:-1]))/np.mean(np.abs(y)), c=1)

def misfit_power_in_layers(model, vec, parameters):
    var = parameters['var']
    split_index = parameters['split_index']
    layer_want = parameters['layer_wanted']
    var_out = model.get_variable(vec, var)
    layer_power0 = np.mean(np.abs(var_out[:split_index,:]))
    layer_power1 = np.mean(np.abs(var_out[split_index:,:]))
    if layer_want == 0:
        return 0.5-0.5*(layer_power0-layer_power1)/(layer_power0+layer_power1)
    else:
        return 0.5+0.5*(layer_power0-layer_power1)/(layer_power0+layer_power1)

def saturation(x, c=1):
    return x/(x+c)

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

def normalize_vec_physical_units(model, vec, var, normalization_value, velocity_type='mm/s'):
    var_data = model.get_variable(vec, var)
    var_data_physical = convert_var_to_physical_units(model, var, var_data, velocity_type=velocity_type)
    return vec * normalization_value/np.max(np.abs(var_data_physical))

def normalize_vec(vec, normalization_value):
    return vec * normalization_value/np.max(np.abs(vec))

def convert_var_to_physical_units(model, var, var_data, velocity_type='mm/s'):
    if var in ['ur', 'uth', 'uph']:
        if velocity_type == 'mm/s':
            var_data = var_data * model.u_star
        elif velocity_type == 'km/yr':
            var_data = var_data * model.u_star * 31556.926
        else:
            raise TypeError('velocity type is not understood')
    elif var in ['br', 'bth', 'bph']:
        var_data = var_data * model.B_star
    elif var == 'r_disp':
        var_data = var_data * model.r_star
    elif var == 'p':
        var_data = var_data * model.P_star
    return var_data

def get_order_th(model, vec, var='uth', cutoff=0.075):
    z = model.get_variable(vec, var)
    ind = np.argmax(np.mean(np.abs(z),axis=1))
    zi = z[ind,1:-1]
    zi[np.abs(zi) < cutoff * np.max(np.abs(zi))] = 0.
    signs = np.sign(zi)
    stripped_signs = signs[np.nonzero(signs)]
    zeros = np.where(stripped_signs[1:] != stripped_signs[:-1])[0]
    return len(zeros)

def get_order_r(model, vec, var='uph', cutoff=0.075):
    z = model.get_variable(vec, var)
    ind = np.argmax(np.mean(np.abs(z),axis=0))
    zi = z[1:-1,ind]
    zi[np.abs(zi) < cutoff * np.max(np.abs(zi))] = 0.
    signs = np.sign(zi)
    stripped_signs = signs[np.nonzero(signs)]
    zeros = np.where(stripped_signs[1:] != stripped_signs[:-1])[0]
    return len(zeros)

def get_Q(val):
    return np.abs(val.imag/(2*val.real))

def get_period(model, val):
    return (2*np.pi/val.imag)*model.t_star/(24.*3600.*365.25)

def get_decay(model, val):
    return (2 * np.pi / val.real) * model.t_star / (24. * 3600. * 365.25)