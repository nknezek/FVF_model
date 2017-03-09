#! /usr/bin/env python3

import sys, os
import FVF_plotlib as fplt
import itertools as it
import numpy as np
import importlib

#%% Import configuration file
default_config = "cfg_make_general"
sys.path.append('../config')
try:
    config_file = os.path.splitext(sys.argv[1])[0]
    cfg = importlib.import_module(config_file)
    print("used config file from command line {0}".format(config_file))
except:
    try:
        config_file = default_config
        cfg = importlib.import_module(config_file)
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" +
              "\n\n ALERT!\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"+
              "Could not load config file from command-line argument, so used default config file " + default_config+"\n")
    except:
        raise ImportError("could not find a configuration file")

#%% Store constant parameters
flog = cfg.flog
R = cfg.R
Omega = cfg.Omega
rho = cfg.rho
mu_0 = cfg.mu_0
g = cfg.g
model_variables = cfg.model_variables
dir_suf = cfg.dir_suf
ep = cfg.ep

# create list of all combinations of iteratable parameters
iter_param_names = ['m', 'Nk', 'Nl', 'h', 'nu', 'eta', 'dCyr', 'B_type', 'Bd', 'Br', 'Bth', 'const', 'Bmax', 'Bmin', 'sin_exp', 'Uphi', 'buoy_type', 'buoy_ratio', 'fvf', 'Bnoise']
iter_params = {}
for name in iter_param_names:
    value = eval('cfg.'+name)
    if type(value) is not list:
        value = [value]
    iter_params[name] = value
varNames = sorted(iter_params)
combinations = [ dict(zip(varNames, prod)) for prod in it.product(*(iter_params[varName] for varName in varNames))]

for c in combinations:
    # Store Parameters for this model run
    m = c['m']
    Nk = c['Nk']
    Nl = c['Nl']
    h = c['h']
    nu = c['nu']
    eta = c['eta']
    dCyr = c['dCyr']
    B_type = c['B_type']
    Bd = c['Bd']
    Br = c['Br']
    Bth = c['Bth']
    const = c['const']
    Bmax = c['Bmax']
    Bmin = c['Bmin']
    sin_exp = c['sin_exp']
    buoy_type = c['buoy_type']
    buoy_ratio = c['buoy_ratio']
    fvf = c['fvf']
    Uphi = c['Uphi']
    Bnoise = c['Bnoise']

    # Directory name to save model
    dir_name = ('../data/k'+str(Nk) + '_l' + str(Nl) +
                '_m{1:.0f}_nu{2:.0e}_{3:.0f}km_{7}N{4:.0f}_{6}B{5:.0f}{8}/'.format(dCyr, m, nu, h/1e3, buoy_ratio*100., Br*1e5, B_type, buoy_type, dir_suf))
    filemodel = 'model.p' # name of model in directory
    fileA = 'A' # name of A matrix data
    fileB = 'B' # name of M matrix data

    #%% Set up Model
    #===============================================================================
    model_parameters = {'Nk': Nk, 'Nl': Nl, 'm': m}
    physical_constants = {'R': R, 'Omega': Omega, 'rho': rho,
                          'h': h, 'nu': nu, 'eta': eta,
                          'mu_0': mu_0, 'g': g}
    model = fvf.Model(model_variables, model_parameters, physical_constants)
    model.set_B_by_type(B_type=B_type, Bd=Bd, Br=Br, Bth=Bth, const=const, Bmax=Bmax, Bmin=Bmin, sin_exp=sin_exp, noise=Bnoise)
    model.set_buoy_by_type(buoy_type=buoy_type, buoy_ratio=buoy_ratio)
    if type(dCyr) == (float or int):
        model.set_CC_skin_depth(dCyr)
    else:
        raise TypeError('dCyr not right type')
    model.set_Uphi(Uphi)
    model.make_operators()

    flog.ensure_dir(dir_name)

    print('done setting up model')

    # %% Save Model info
    #==============================================================================
    fplt.plot_buoy_struct(model, dir_name=dir_name)
    print('plotted buoyancy structure')
    fplt.plot_B(model, dir_name=dir_name)
    print('plotted background magnetic field structure')
    fplt.plot_Uphi(model, dir_name=dir_name)
    print('plotted background Uphi structure')

    logger = flog.setup_custom_logger(dir_name=dir_name, filename='model.log')
    logger.info('\n' +
    "Model Information:\n" +
    "from config file: {0}".format(config_file) + '\n\n' +
    'm = ' + str(model.m) + '\n' +
    'Nk = ' + str(model.Nk) + '\n' +
    'Nl = ' + str(model.Nl) + '\n' +
    'R = ' + str(model.R) + '\n' +
    'h = ' + str(model.h) + '\n' +
    'Omega = ' + str(model.Omega) + '\n' +
    'rho = ' + str(model.rho) + '\n' +
    'nu = ' + str(model.nu) + '\n' +
    'eta = ' + str(model.eta) + '\n' +
    'mu_0 = ' + str(model.mu_0) + '\n' +
    'g = ' + str(model.g) + '\n' +
    'dCyr = ' + str(dCyr) + '\n' +
    'B_Type = ' + str(B_type) + '\n' +
    'Bd = ' + str(Bd) + '\n' +
    'Br = ' + str(model.Br.max()) + ' to ' + str(model.Br.min()) + '\n' +
    'Bth = ' + str(model.Bth.max()) + ' to ' + str(model.Bth.min()) + '\n' +
    'Uph = ' + str(model.Uphi.max()) + ' to ' + str(model.Uphi.min()) + '\n' +
    'buoy_type = ' + str(buoy_type) + '\n' +
    'buoy_ratio = ' + str(buoy_ratio) +'\n' +
    'model variables = ' + str(model.model_variables) + '\n'
    )
    print('model will be saved in ' + str(dir_name))

    #%% Make matricies used for later analysis
    #==============================================================================
    model.make_Bobs()
    print('created Bobs matrix')
#    model.make_d2Mat()
#    print('created D3sq matrix')
#    model.make_dthMat()
#    print('created dth matrix')

    #%% Save Model Information
    #==============================================================================
    model.save_model(dir_name + filemodel)
    print('saved model to ' + str(dir_name))

    #%% Create Matrices
    #===============================================================================
    model.make_B()
    print('created B matrix')
    epB = np.min(np.abs(model.B.data[np.nonzero(model.B.data)]))*ep
    model.save_mat_PETSc(dir_name+fileB+'.dat', model.B.toPETSc(epsilon=epB))
    print('saved PETSc B matrix ' + str(dir_name))
    model.make_A()
    print('created A matrix')
    epA = np.min(np.abs(model.A.data[np.nonzero(model.A.data)]))*ep
    model.save_mat_PETSc(dir_name+fileA+str(dCyr)+'.dat', model.A.toPETSc(epsilon=epA))
    print('saved PETSc A matrix for dCyr = {0} to '.format(dCyr) + str(dir_name))

    # import pickle; pickle.dump(model.A.todense(), open('A_py3matrix.p','wb')); pickle.dump(model.B.todense(), open('B_py3matrix.p','wb'));