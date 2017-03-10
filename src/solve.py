#! /usr/bin/env python3

import FVF_loglib as flog
import FVF_plotlib as fplt
import FVF_analysis as fana
import pickle as pkl
import numpy as np
import itertools as it
from datetime import datetime
import sys, os
import importlib
import shutil
import slepc4py
slepc4py.init(sys.argv)
from petsc4py import PETSc
from slepc4py import SLEPc
opts = PETSc.Options()

# Import configuration file
default_config = "cfg_solve_general"
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
              "\n\n ALERT!\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+
              "\n\nCould not load config file from command-line argument, so used default config file " + default_config+"\n")
    except:
        raise ImportError("Could not import config file")

# function to find nearest skin-depth value for a wave period
def find_closest_CC(Target, dCyr_list):
    dist = [abs(x-abs(Target)) for x in dCyr_list]
    return dCyr_list[dist.index(min(dist))]


# Import constant parameters from config file
T_list = cfg.T_list
wt_T=cfg.wt_T
delta_T = cfg.delta_T

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
dCyr_list = cfg.dCyr_list
savefile = cfg.savefile

plot_robinson = cfg.plot_robinson
plot_B_obs = cfg.plot_B_obs
plot_vel = cfg.plot_vel

tol = cfg.tol

# Iterate over parameters that can vary
iter_param_names = ['data_dir', 'oscillate']
iter_params = {}
for name in iter_param_names:
    value = eval('cfg.'+name)
    if type(value) is not list:
        value = [value]
    iter_params[name] = value
varNames = sorted(iter_params)
combinations = [ dict(zip(varNames, prod)) for prod in it.product(*(iter_params[varName] for varName in varNames))]

# Store main output directory
out_dir_base = '../output/{0}_{1}/'.format(datetime.today().strftime("%Y-%m-%d_%H-%M-%S"), config_file[9:])
# Set up logger
logger = flog.setup_custom_logger(dir_name=out_dir_base, filename='run.log')

# Store config file for later reference
logger.info('used config file {0}.py'.format(config_file))
shutil.copyfile('../config/'+config_file + '.py', out_dir_base+config_file+'.py')
logger.info('Main output directory set to {0}'.format(out_dir_base))


for cnum, c in enumerate(combinations):
    data_dir = c['data_dir']
    oscillate = c['oscillate']

    # Get information on layer thickness, buoyancy frequency, and magnetic field from data directory
    data_dir_info = data_dir.strip('/').split('_')
    m = data_dir_info[2]
    H = data_dir_info[4]
    N = data_dir_info[5]
    B = data_dir_info[6]

    # Set up directory to store solution data
    if(len(combinations) > 1):
        out_dir = out_dir_base + '{0}/{1}/{2}/{3}/'.format(m, N, B, H)
    else:
        out_dir = out_dir_base
    logger.info('\n\nParameter Set {0}, m = {1}, H = {2}, B = {3}, N = {4} \n'.format(cnum, m, H, B, N))

    for tnum,T in enumerate(T_list):
        logger.info('\n\nT guess {0}, Target {1:.1f} yrs, DelT {2:.1f} yrs'.format(tnum, T, delta_T))

        # Convert Time in years to model frequency
        t_star = (23.9345*3600)/(2*np.pi)
        Target_j = 2*np.pi/(T*365.25*24*3600/t_star)*1j
        Target = Target_j + Target_j*1j/(2*target_Q)

        # Find which CC matrix to use
        dCyr_use = find_closest_CC(T, dCyr_list)
        logger.info('{0} dCyr used'.format(dCyr_use))

        # Set out_dir
        if(len(T_list)>1):
            out_dir_T = out_dir+'{0:.2f}/'.format(T)
        else:
            out_dir_T = out_dir
        flog.ensure_dir(out_dir_T)
        logger.info('out_dir set to {0}'.format(out_dir_T))

        # %% Load Matrices and model from files
        #==============================================================================
        try:
            viewer = PETSc.Viewer().createBinary(data_dir+fileA+str(dCyr_use)+'.dat', 'r')
            A = PETSc.Mat().load(viewer)
            viewer = PETSc.Viewer().createBinary(data_dir+fileB+'.dat', 'r')
            B = PETSc.Mat().load(viewer)
            try:
                model = pkl.load(open(data_dir+filemodel,'rb'))
            except:
                model = pkl.load(open(data_dir+filemodel,'rb'),encoding='latin1')
            logger.info('A'+str(dCyr_use)+' matrix used')
            logger.info('matrices and model loaded into memory from ' + data_dir)
        except:
            logger.error( "Could not load matrices from file", exc_info=1)
            break


        # %% Set up SLEPc Solver
        #==============================================================================
        try:
            EPS = SLEPc.EPS().create()
            EPS.setDimensions(nev, PETSc.DECIDE)
            EPS.setOperators(A, B)
            EPS.setType(EPS.Type.KRYLOVSCHUR)
            EPS.setProblemType(SLEPc.EPS.ProblemType.PGNHEP)
            EPS.setTarget(Target)
            EPS.setWhichEigenpairs(EPS.Which.TARGET_MAGNITUDE)
            EPS.setTolerances(tol)
            EPS.setFromOptions()
            ST = EPS.getST()
            ST.setType(SLEPc.ST.Type.SINVERT)
            logger.info('solver set up, Period = {0:.1f}, nev = {1}'.format(T, nev))
            logger.info('eigenvalue target = {0:.1e}'.format(Target))
        except:
            logger.error( "Could not set up SLEPc solver ", exc_info=1)
            break

        # %% Solve Problem
        #==============================================================================
        try:
            EPS.solve()
            logger.info('problem solved')
        except:
            logger.error("Could not solve problem.")
            break
        try:
            # Save Computed Solutions
            conv = EPS.getConverged()
            logger.info('{0} eigenvalues converged'.format(conv))
            vals = []
            vecs = []
            for ind in range(conv):
                vs, ws = PETSc.Mat.getVecs(A)
                v = EPS.getEigenpair(ind, ws)
                vals.append(v)
                vecs.append(ws.getArray())
            Periods = (2*np.pi/np.array([x.imag for x in vals]))*model.t_star/(24.*3600.*365.25)
            Period_max = Periods.max()
            Period_min = Periods.min()
            logger.info('min Period = {0:.1f}yrs, max Period = {1:.1f}yrs'.format(Period_min, Period_max))
        except:
            logger.error("Could not get converged eigenvalues.", exc_info=1)
            break

        #%% Filter Solutions
        #==============================================================================
        try:
            logger.info('Filtering Eigenvalues:')

            # Filter by fit to given parameter choices
            fvals, fvecs = fana.filter_by_misfit(model, vals, vecs, num_to_keep, target_T=T, target_Q=target_Q,
                                                target_th_order=target_th_order, target_r_order=target_r_order,
                                                target_region='equator', eq_cutoff=0.5, target_symmetric=True,
                                                oscillate=False, wt_T=wt_T, wt_Q=wt_Q, wt_r_order=wt_r_order,
                                                 wt_th_order=wt_th_order, wt_region=wt_region, wt_sym=wt_sym)

            #%% Convert Eigenvectors to ur real
            fvecs_tmp = []
            for vec in fvecs:
                fvecs_tmp.append(fana.shift_vec_real(model, vec, var=real_var))
            fvecs = fvecs_tmp
            logger.info('\t{0} eigenvectors shifted so that {1} is real to plot'.format(len(fvecs), real_var))
        except:
            logger.error("Problem Filtering Eigenvalues.", exc_info=1)
            break

        # %% Save Filtered Eigenvectors
        #==============================================================================
        try:
            if savefile:
                pkl.dump({'vals': fvals, 'vecs': fvecs, 'model':model},open(out_dir_T + savefile, 'wb'))
                logger.info('vals and vecs saved to ' + out_dir_T + savefile)
        except:
            logger.error("Problem Saving Filtered Eigenvalues.", exc_info=1)
            break

        # %% Plot Filtered Eigenvectors
        #==============================================================================
        try:
            logger.info('Plotting:')
            for ind, (val, vec) in enumerate(zip(fvals,fvecs)):
                Period = (2*np.pi/val.imag)*model.t_star/(24.*3600.*365.25)
                Decay = (2*np.pi/val.real)*model.t_star/(24.*3600.*365.25)
                Q = abs(val.imag/(2*val.real))
                if abs(Period) < 1.0:
                    title = ('T{1:.2f}dys_Q{2:.2f}_{0:.0f}'.format(ind, Period*365.25, Q))
                else:
                    title = ('T{1:.2f}yrs_Q{2:.2f}_{0:.0f}'.format(ind, Period, Q))
                if plot_vel:
                    if (model.Nk > 1):
                        fplt.plot_pcolormesh_rth(model, val, vec, dir_name=out_dir_T, title=title, physical_units=True, oscillate_values=oscillate)
                    else:
                        if('br' in model.model_variables):
                            fplt.plot_1D(model, vec, val, ind, dir_name=out_dir_T, title=title)
                        else:
                            fplt.plot_1D_noB(model, vec, val, ind, dir_name=out_dir_T, title=title)
                if plot_robinson:
                    fplt.plot_robinson(model, vec, model.m, oscillate=oscillate, dir_name=out_dir_T, title=str(ind)+'_T{0:.2f}yrs_'.format(Period)+'Divergence')
                if plot_B_obs:
                    fplt.plot_B_obs(model, vec, model.m, oscillate=oscillate, dir_name=out_dir_T, title=title+'_Bperturb')
                logger.info('\t plotted ind={0}, T={1:.2f}yrs (eig={2:.2e})'.format(ind, Period, val))
            logger.info('run complete')
        except:
            logger.error("Problem Plotting Eigenvalues.", exc_info=1)
            break
