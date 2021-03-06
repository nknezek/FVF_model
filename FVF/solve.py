#! /usr/bin/env python3

import multiprocess as mp
import numpy as np
import itertools as it
from datetime import datetime
import sys, os
import importlib
import shutil
from . import notify as fvn

# Import configuration file
default_config = "cfg_solve_general"
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
              "\n\n ALERT!\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+
              "\n\nCould not load config file from command-line argument, so used default config file " + default_config+"\n")
    except:
        raise ImportError("Could not import config file")


# Import constant parameters from config file
num_solutions_to_calculate = cfg.num_solutions_to_calculate
filemodel = cfg.filemodel
fileA = cfg.fileA
fileB = cfg.fileB
savefile = cfg.savefile
tol = cfg.tol

# Iterate over parameters that can vary
iter_param_names = ['data_dir', 'T_list']
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

# Store main output directory
out_dir_base = '../output/{0}_{1}/'.format(datetime.today().strftime("%Y-%m-%d_%H-%M-%S"), config_file.split('_')[-1])


def solve_for_combo(c):
    import time
    time_start = -time.time()
    import FVF_loglib as flog
    import FVF_plotlib as fplt
    from . import analyze as fana
    from . import utilities as util
    import petsc4py
    petsc4py.init()
    import slepc4py
    slepc4py.init()
    import dill

    print('working on combination {0}/{1}'.format(c['iter_num'], c['total_iter']))

    data_dir = c['data_dir']
    T = c['T_list']

    # Set up directory to store solution data
    try:
        out_dir = futil.get_out_dir(out_dir_base, data_dir, len(cfg.data_dir), T, len(cfg.T_list))
        futil.ensure_dir(out_dir)
    except:
        print('problem setting out directory')
        return

    # Set up logger
    try:
        logger = flog.setup_custom_logger(dir_name=out_dir, filename='run.log', verbose=cfg.verbose)
    except:
        print('problem creating logger')
        return

    try:
        # Store config file for later reference
        logger.info('used config file {0}.py'.format(config_file))
        futil.store_config_file(config_file, out_dir_base)

        logger.info('Main output directory set to {0}'.format(out_dir_base))
        logger.info('Output subdirectory set to {0}'.format(out_dir))

        # Convert Time in years to model frequency
        t_star = (23.9345*3600)/(2*np.pi)
        Target_j = 2*np.pi/(T*365.25*24*3600/t_star)*1j
        Target = Target_j + Target_j*1j/(2*cfg.target_Q)

        # Find which CC matrix to use
        dCyr_list = futil.find_available_skin_depths(data_dir)
        dCyr_use = futil.find_closest_CC(T, dCyr_list)
        logger.info('{0} dCyr used'.format(dCyr_use))
    except:
        problem = "problem storing config file or finding correct magnetic skin depth"
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))

        return

    # %% Load Matrices and model from files
    #==============================================================================
    try:
        viewer = petsc4py.PETSc.Viewer().createBinary(data_dir+fileA+str(dCyr_use)+'.dat', 'r')
        A = petsc4py.PETSc.Mat().load(viewer)
        viewer = petsc4py.PETSc.Viewer().createBinary(data_dir+fileB+'.dat', 'r')
        B = petsc4py.PETSc.Mat().load(viewer)
        try:
            model = dill.load(open(data_dir+filemodel,'rb'))
        except:
            model = dill.load(open(data_dir+filemodel,'rb'),encoding='latin1')
        logger.info('A'+str(dCyr_use)+' matrix used')
        logger.info('matrices and model loaded into memory from ' + data_dir)
    except:
        problem = "Problem loading matrices from file."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
        return

    # %% Set up SLEPc Solver
    #==============================================================================
    try:
        EPS = slepc4py.SLEPc.EPS().create()
        EPS.setDimensions(num_solutions_to_calculate, petsc4py.PETSc.DECIDE)
        EPS.setOperators(A, B)
        EPS.setType(EPS.Type.KRYLOVSCHUR)
        EPS.setProblemType(slepc4py.SLEPc.EPS.ProblemType.PGNHEP)
        EPS.setTarget(Target)
        EPS.setWhichEigenpairs(EPS.Which.TARGET_MAGNITUDE)
        EPS.setTolerances(tol)
        EPS.setFromOptions()
        ST = EPS.getST()
        ST.setType(slepc4py.SLEPc.ST.Type.SINVERT)
        logger.info('Solver set up, Target Period = {0:.1f}, Number of solutions to calculate = {1}'.format(T, num_solutions_to_calculate))
    except:
        problem = "Problem setting up SLEPc Solver."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
        return

    # %% Solve Problem
    #==============================================================================
    try:
        EPS.solve()
        logger.info('problem solved')
    except:
        problem = "Problem Solving Eigenvalues."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
        return
    try:
        # Save Computed Solutions
        conv = EPS.getConverged()
        logger.info('{0} eigenvalues converged'.format(conv))
        vals = []
        vecs = []
        for ind in range(conv):
            vs, ws = petsc4py.PETSc.Mat.getVecs(A)
            v = EPS.getEigenpair(ind, ws)
            vals.append(v)
            vecs.append(ws.getArray())
        Periods = (2*np.pi/np.array([x.imag for x in vals]))*model.t_star/(24.*3600.*365.25)
        Period_max = Periods.max()
        Period_min = Periods.min()
        logger.info('min Period = {0:.1f}yrs, max Period = {1:.1f}yrs'.format(Period_min, Period_max))
    except:
        problem = "Problem Computing Eigenvalues."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
        return

    #%% Filter Solutions
    #==============================================================================
    try:
        logger.info('Filtering Eigenvalues:')

        # filter results to keep only those that satisfy requirements specified in filter_dict
        fvals, fvecs = fana.filter_results(model, vals, vecs, cfg.filter_dict)

        # Sort by fit to given parameter choices
        svals, svecs = fana.sort_by_total_misfit(model, fvals, fvecs, cfg.sort_dict)

    except:
        try:
            problem = "Problem Saving Eigenvalues."
            logger.error(problem, exc_info=1)
            print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
            svals = vals
            svecs = vecs
        except:
            svals = []
            svecs = []
            problem = "Problem Saving Eigenvalues."
            logger.error(problem, exc_info=1)
            print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))

    # %% Save Filtered Eigenvectors
    #==============================================================================
    try:
        if savefile:
            dill.dump({'vals': svals, 'vecs': svecs, 'model':model},open(out_dir + savefile, 'wb'))
            logger.info('saved {0:d} vals and vecs saved to '.format(len(svals)) + out_dir + savefile)
    except:
        problem = "Problem Saving Eigenvalues."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))
        return
    # %% Plot Filtered Eigenvectors
    #==============================================================================
    try:
        logger.info('Plotting:')

        for ind in range(min(cfg.num_solutions_to_plot,len(svals))):
            val = svals[ind]
            vec = fana.shift_vec_real(model, svecs[ind], var='vth')
            vec = fana.normalize_vec(vec, 10)
            Period = fana.get_period(model, val)
            Q = fana.get_Q(val)
            r_ord = fana.get_order_r(model, vec)
            th_ord = fana.get_order_th(model, vec)
            if abs(Period) < 1.0:
                title = ('{0:03d} k={3}, l={4}, m={5}, T={1:.2f}dys, Q={2:.2f}'.format(ind, Period*365.25, Q, r_ord, th_ord, model.m)
                         + '\nH={:.0f}km, B={:.2f}mT, N={:.2f} pvbc'.format(model.h/1e3, np.mean(model.Br)*model.B_star*1e3, np.mean(model.N)))
            else:
                title = ('{0:03d} k={3}, l={4}, m={5}, T={1:.2f}yrs, Q={2:.2f}'.format(ind, Period, Q, r_ord, th_ord, model.m)
                         + '\nH={:.0f}km, B={:.2f}mT, N={:.2f} pvbc'.format(model.h/1e3, np.mean(model.Br)*model.B_star*1e3, np.mean(model.N)))
            if (model.Nk > 1):
                # fplt.plot_fast_solution(model, vec, title=title, dir_name=out_dir)
                fplt.plot_full_solution(model, val, vec, title=title, dir_name=out_dir, layer_boundary=3380)
            else:
                if('br' in model.model_variables):
                    fplt.plot_1D(model, vec, val, ind, dir_name=out_dir, title=title)
                else:
                    fplt.plot_1D_noB(model, vec, val, ind, dir_name=out_dir, title=title)
            logger.info('\t plotted ind={0}, T={1:.2f}yrs (eig={2:.2e})'.format(ind, Period, val))
        logger.info('run complete')
    except:
        problem = "Problem Plotting Eigenvalues."
        logger.error(problem, exc_info=1)
        print('combination {0}/{1} broke, {2}'.format(c['iter_num'], c['total_iter'], problem))

        return
    dtime = (time_start + time.time())/60.
    print('done with combination {0}/{1} in {2:.2f}min'.format(c['iter_num'], c['total_iter'], dtime))

if __name__ == '__main__':
    try:
        if cfg.num_threads is None:
            procs = mp.cpu_count()
        else:
            procs = cfg.num_threads
        mainlogging = mp.log_to_stderr()
        p = mp.Pool(processes=procs)
        p.map(solve_for_combo, combinations)
        time = datetime.today().ctime()
        message = 'finished run of {0} parameter sets using {1} at {2}'.format(len(combinations), config_file, time)
        print(message)
        if cfg.notify_me_by_text:
            cli = fvn.MessageClient()
            cli.send_message(message)
    except:
        message = "something went wrong with the multi-threading"
        print(message)
        if cfg.notify_me_by_text:
            cli = fvn.MessageClient()
            cli.send_message(message)

