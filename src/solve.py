#! /usr/bin/env python3

import multiprocess as mp
import numpy as np
import itertools as it
from datetime import datetime
import sys, os
import importlib
import shutil
import FVF_notify as fvn


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
wt_T=cfg.wt_T
nev = cfg.nev
target_Q = cfg.target_Q
wt_Q=cfg.wt_Q
target_th_order = cfg.target_th_order
wt_th_order=cfg.wt_th_order
target_r_order = cfg.target_r_order
wt_r_order=cfg.wt_r_order
num_to_keep = cfg.num_to_keep
wt_region=cfg.wt_region
wt_sym=cfg.wt_sym
eq_split = cfg.eq_split

eq_var = cfg.eq_var
real_var = cfg.real_var

filemodel = cfg.filemodel
fileA = cfg.fileA
fileB = cfg.fileB
savefile = cfg.savefile

plot_robinson = cfg.plot_robinson
plot_B_obs = cfg.plot_B_obs
plot_vel = cfg.plot_vel

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
    import FVF_loglib as flog
    import FVF_plotlib as fplt
    import FVF_analysis as fana
    import FVF_utilities as futil
    import petsc4py
    petsc4py.init()
    import slepc4py
    slepc4py.init()
    import dill

    print('working on combination {0}/{1}'.format(c['iter_num'], c['total_iter']))

    data_dir = c['data_dir']
    T = c['T_list']

    # Set up directory to store solution data
    out_dir = futil.get_out_dir(out_dir_base, data_dir, len(cfg.data_dir), T, len(cfg.T_list))
    futil.ensure_dir(out_dir)

    # Set up logger
    logger = flog.setup_custom_logger(dir_name=out_dir, filename='run.log', verbose=cfg.verbose)

    # Store config file for later reference
    logger.info('used config file {0}.py'.format(config_file))
    futil.store_config_file(config_file, out_dir_base)

    logger.info('Main output directory set to {0}'.format(out_dir_base))
    logger.info('Output subdirectory set to {0}'.format(out_dir))

    # Convert Time in years to model frequency
    t_star = (23.9345*3600)/(2*np.pi)
    Target_j = 2*np.pi/(T*365.25*24*3600/t_star)*1j
    Target = Target_j + Target_j*1j/(2*target_Q)

    # Find which CC matrix to use
    dCyr_list = futil.find_available_skin_depths(data_dir)
    dCyr_use = futil.find_closest_CC(T, dCyr_list)
    logger.info('{0} dCyr used'.format(dCyr_use))


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
        logger.error( "Could not load matrices from file", exc_info=1)


    # %% Set up SLEPc Solver
    #==============================================================================
    try:
        EPS = slepc4py.SLEPc.EPS().create()
        EPS.setDimensions(nev, petsc4py.PETSc.DECIDE)
        EPS.setOperators(A, B)
        EPS.setType(EPS.Type.KRYLOVSCHUR)
        EPS.setProblemType(slepc4py.SLEPc.EPS.ProblemType.PGNHEP)
        EPS.setTarget(Target)
        EPS.setWhichEigenpairs(EPS.Which.TARGET_MAGNITUDE)
        EPS.setTolerances(tol)
        EPS.setFromOptions()
        ST = EPS.getST()
        ST.setType(slepc4py.SLEPc.ST.Type.SINVERT)
        logger.info('solver set up, Period = {0:.1f}, nev = {1}'.format(T, nev))
        logger.info('eigenvalue target = {0:.1e}'.format(Target))
    except:
        logger.error( "Could not set up SLEPc solver ", exc_info=1)

    # %% Solve Problem
    #==============================================================================
    try:
        EPS.solve()
        logger.info('problem solved')
    except:
        logger.error("Could not solve problem.")
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
        logger.error("Could not get converged eigenvalues.", exc_info=1)

    #%% Filter Solutions
    #==============================================================================
    try:
        logger.info('Filtering Eigenvalues:')
        #%% Convert Eigenvectors to uph real
        vecs_tmp = []
        for vec in vecs:
            vecs_tmp.append(fana.shift_vec_real(model, vec, var=real_var))
        vecs = vecs_tmp
        logger.info('\t{0} eigenvectors shifted so that {1} is real to plot'.format(len(vecs), real_var))

        # filter by number of r/th zero crossings
        fvals, fvecs = fana.filter_by_rth_zeros(model, vals, vecs)
        # print('filtered out {0} eigenvectors because of too many zeros'.format(len(vals)-len(fvals)))
        # Filter by fit to given parameter choices
        fvals, fvecs = fana.filter_by_misfit(model, fvals, fvecs, num_to_keep, target_T=T, target_Q=target_Q,
                                            target_th_order=target_th_order, target_r_order=target_r_order,
                                            target_region=cfg.target_region, eq_cutoff=0.5, target_symmetric=True,
                                             wt_T=wt_T, wt_Q=wt_Q, wt_r_order=wt_r_order,
                                             wt_th_order=wt_th_order, wt_region=wt_region, wt_sym=wt_sym,
                                             wt_r_sm=cfg.wt_r_sm, wt_th_sm=cfg.wt_th_sm,
                                             th_ord_var=cfg.th_ord_var, r_ord_var=cfg.r_ord_var, eq_var=cfg.eq_var)

    except:
        logger.error("Problem Filtering Eigenvalues.", exc_info=1)

    # %% Save Filtered Eigenvectors
    #==============================================================================
    try:
        if savefile:
            dill.dump({'vals': fvals, 'vecs': fvecs, 'model':model},open(out_dir + savefile, 'wb'))
            logger.info('vals and vecs saved to ' + out_dir + savefile)
    except:
        logger.error("Problem Saving Filtered Eigenvalues.", exc_info=1)

    # %% Plot Filtered Eigenvectors
    #==============================================================================
    try:
        logger.info('Plotting:')
        for ind, (val, vec) in enumerate(zip(fvals,fvecs)):
            Period = (2*np.pi/val.imag)*model.t_star/(24.*3600.*365.25)
            Decay = (2*np.pi/val.real)*model.t_star/(24.*3600.*365.25)
            Q = abs(val.imag/(2*val.real))
            r_ord = fana.get_r_zero_crossings(model, vec, var=cfg.r_ord_var)
            th_ord = fana.get_theta_zero_crossings(model, vec, var=cfg.th_ord_var)

            if abs(Period) < 1.0:
                title = ('m={5}, l={4}, k={3}, T{1:.2f}dys_Q{2:.2f}_{0:.0f}'.format(ind, Period*365.25, Q, r_ord, th_ord, model.m))
            else:
                title = ('m={5}, l={4}, k={3}, T{1:.2f}yrs_Q{2:.2f}_{0:.0f}'.format(ind, Period, Q, r_ord, th_ord, model.m))
            if plot_vel:
                if (model.Nk > 1):
                    fplt.plot_pcolormesh_rth(model, val, vec, dir_name=out_dir, title=title, physical_units=True)
                else:
                    if('br' in model.model_variables):
                        fplt.plot_1D(model, vec, val, ind, dir_name=out_dir, title=title)
                    else:
                        fplt.plot_1D_noB(model, vec, val, ind, dir_name=out_dir, title=title)
            if plot_robinson:
                fplt.plot_robinson(model, vec, model.m,  dir_name=out_dir, title=str(ind)+'_T{0:.2f}yrs_'.format(Period)+'Divergence')
            if plot_B_obs:
                fplt.plot_B_obs(model, vec, model.m,  dir_name=out_dir, title=title+'_Bperturb')
            logger.info('\t plotted ind={0}, T={1:.2f}yrs (eig={2:.2e})'.format(ind, Period, val))
        logger.info('run complete')
    except:
        logger.error("Problem Plotting Eigenvalues.", exc_info=1)
    print('done with combination {0}/{1}'.format(c['iter_num'], c['total_iter']))

if __name__ == '__main__':
    procs = mp.cpu_count()
    p = mp.Pool(processes=procs)
    p.map(solve_for_combo, combinations)
    time = datetime.today().ctime()
    message = 'finished solvimg {0} parameter sets using {1} at {2}'.format(len(combinations), config_file, time)
    print(message)
    if cfg.notify_me_by_text:
        cli = fvn.MessageClient()
        cli.send_message(message)