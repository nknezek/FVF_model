#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [20, -20, 46, -46, 91, -91, 181, -181]

nev = 100
num_to_keep = 30

eq_split = 0.5
target_r_order = 2
target_th_order = 2
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 1.
wt_r_order = 1e1
wt_th_order = 1.
wt_region = 1.
wt_sym = 1.
wt_r_sm = 1e2 # smooth variation in r-direciton
wt_th_sm = 1e2 # smooth variation in th-direction

eq_var = 'uth'
r_ord_var = 'uth'
th_ord_var = 'uth'
real_var = 'uth'

filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

plot_robinson = False
plot_B_obs = False
plot_vel = True

tol = 1e-6

data_dir = [
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_20k_200l/',
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_0.10mTrnse_20k_200l/',
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_0.30mTrnse_20k_200l/',
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_0.50mTrnse_20k_200l/',
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_0.70mTrnse_20k_200l/',
    '../data/m0_140km_1.00N_abs_dipole_0.62mTBd_0.10mTrnse_20k_200l/',
    '../data/m0_140km_1.00N_constant_0.62mTBr_20k_200l/',
]

notify_me_by_text = True
verbose = False
num_threads = None