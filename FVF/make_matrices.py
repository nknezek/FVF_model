#! /usr/bin/env python3

import sys, os
import itertools as it
import numpy as np
import importlib
import multiprocess as mp
from . import notify as fvn
import datetime

#%% Import configuration file
default_config = "cfg_make_general"
sys.path.append('../config')
try:
    config_file = os.path.splitext(sys.argv[1])[0]
    cfg = importlib.import_module(config_file)
    if cfg.verbose:
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
R = cfg.R
Omega = cfg.Omega
rho = cfg.rho
mu_0 = cfg.mu_0
g = cfg.g
model_variables = cfg.model_variables
dir_suf = cfg.dir_suf
ep = cfg.ep

# create list of all combinations of iteratable parameters
iter_param_names = ['m', 'Nk', 'Nl', 'h', 'nu', 'eta', 'skin_depth_period',
                    'B_type', 'Bd', 'Br', 'Bth','Bph','Brconst', 'Brnoise','Brmult', 'Bthconst', 'Bthnoise', 'Bthmult','use_Bth',
                    'Vphi', 'buoyancy_type', 'N', 'model_type']
iter_params = {}
for name in iter_param_names:
    value = eval('cfg.'+name)
    if type(value) is not list:
        value = [value]
    iter_params[name] = value
varNames = sorted(iter_params)
combinations = [ dict(zip(varNames, prod)) for prod in it.product(*(iter_params[varName] for varName in varNames))]

def append_iter_num_to_combinations(combinations):
    for i,c in enumerate(combinations):
        c['iter_num']=i+1
        c['total_iter'] = len(combinations)

append_iter_num_to_combinations(combinations)
def make_matrix(c):
    '''
    :param c: dictionary of parameters
    :return:
    '''
    import FVF_loglib as flog
    import FVF_plotlib as fplt
    from . import utilities as util

    exec('import {0} as fvf'.format(c['model_type']))

    print('working on combination {0}/{1}'.format(c['iter_num'], c['total_iter']))

    # Store Parameters for this model run into local variables to make things easier
    m = c['m']
    Nk = c['Nk']
    Nl = c['Nl']
    h = c['h']
    nu = c['nu']
    eta = c['eta']
    skin_depth_period = c['skin_depth_period']
    buoyancy_type = c['buoyancy_type']
    N = c['N']
    Vphi = c['Vphi']

    # Directory name to save model
    dir_name = futil.get_directory_name(c, dir_suf=cfg.dir_suf)
    filemodel = 'model.p' # name of model in directory
    fileA = 'A' # name of A matrix data
    fileB = 'B' # name of M matrix data

    #%% Set up Model
    #===============================================================================
    model_parameters = {'Nk': Nk, 'Nl': Nl, 'm': m}
    physical_constants = {'R': R, 'Omega': Omega, 'rho': rho,
                          'h': h, 'nu': nu, 'eta': eta,
                          'mu_0': mu_0, 'g': g}
    # have to call fvf from the locals() dict direction, as for some reason exec during multiprocess doesn't update the local namespace
    model = locals()['fvf'].Model(model_variables, model_parameters, physical_constants)
    model.set_B_by_type(c['B_type'], Bd=c['Bd'], Br=c['Br'], Bth=c['Bth'], Bph=c['Bph'], use_Bth=c['use_Bth'],
                      Brconst=c['Brconst'], Brnoise=c['Brnoise'], Brmult=c['Brmult'], Bthconst=c['Bthconst'], Bthnoise=c['Bthnoise'], Bthmult=c['Bthmult'])
    model.set_buoyancy_by_type(buoyancy_type=buoyancy_type, N=N)
    if type(skin_depth_period) == (float or int):
        model.set_mag_skin_depth(skin_depth_period)
    else:
        raise TypeError('skin_depth_period not right type')
    model.set_Vphi(Vphi)
    model.make_operators()

    futil.ensure_dir(dir_name)

    if cfg.verbose:
        print('done setting up model')

    # %% Save Model info
    #==============================================================================
    fplt.plot_buoyancy_struct(model, dir_name=dir_name)
    if cfg.verbose:
        print('plotted buoyancy structure')
    fplt.plot_B(model, dir_name=dir_name)
    if cfg.verbose:
        print('plotted background magnetic field structure')
    fplt.plot_Vphi(model, dir_name=dir_name)
    if cfg.verbose:
        print('plotted background Vphi structure')

    logger = flog.setup_custom_logger(dir_name=dir_name, filename='model.log', verbose=cfg.verbose)
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
    'skin_depth_period = ' + str(skin_depth_period) + '\n' +
    'B_Type = ' + str(c['B_type']) + '\n' +
    'Bd = ' + str(c['Bd']) + '\n' +
    'Br = ' + str(model.Br.max()) + ' to ' + str(model.Br.min()) + '\n' +
    'Bth = ' + str(model.Bth.max()) + ' to ' + str(model.Bth.min()) + '\n' +
    'Uph = ' + str(model.Vphi.max()) + ' to ' + str(model.Vphi.min()) + '\n' +
    'buoyancy_type = ' + str(buoyancy_type) + '\n' +
    'N = ' + str(N) +'\n' +
    'model variables = ' + str(model.model_variables) + '\n'
    )
    if cfg.verbose:
        print('model will be saved in ' + str(dir_name))

    #%% Make matricies used for later analysis
    #==============================================================================
    model.make_Bobs()
    if cfg.verbose:
        print('created Bobs matrix')

    #%% Save Model Information
    #==============================================================================
    model.save_model(dir_name + filemodel)
    if cfg.verbose:
        print('saved model to ' + str(dir_name))

    #%% Create Matrices
    #===============================================================================
    model.make_B()
    if cfg.verbose:
        print('created B matrix')
    epB = np.min(np.abs(model.B.data[np.nonzero(model.B.data)]))*ep
    model.save_mat_PETSc(dir_name+fileB+'.dat', model.B.toPETSc(epsilon=epB))
    if cfg.verbose:
        print('saved PETSc B matrix ' + str(dir_name))
    model.make_A()
    if cfg.verbose:
        print('created A matrix')
    epA = np.min(np.abs(model.A.data[np.nonzero(model.A.data)]))*ep
    model.save_mat_PETSc(dir_name+fileA+str(skin_depth_period)+'.dat', model.A.toPETSc(epsilon=epA))
    if cfg.verbose:
        print('saved PETSc A matrix for skin_depth_period = {0} to '.format(skin_depth_period) + str(dir_name))
    print('done with combination {0}/{1}'.format(c['iter_num'], c['total_iter']))

    return

if __name__ == '__main__':
    if cfg.num_threads is None:
        procs = mp.cpu_count()
    else:
        procs = cfg.num_threads
    p = mp.Pool(processes=procs)
    p.map(make_matrix, combinations)
    time = datetime.datetime.today().ctime()
    message = 'finished making {0} matrices using {1} at {2}'.format(len(combinations), config_file, time)
    print(message)
    if cfg.notify_me_by_text:
        cli = fvn.MessageClient()
        cli.send_message(message)