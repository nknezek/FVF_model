#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

# T_list = [25/365.25, 50/365.25, 100/365.25, 200/365.25, -25/365.25, -50/365.25, -100/365.25, -200/365.25]
# T_list = [3000, 6000]
nev = 500
num_to_keep = 50

eq_split = 0.5
target_r_order = 1
target_th_order = 1
target_Q = 1e6
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 1.
wt_r_order = 1.
wt_th_order = 1.
wt_region = 10.
wt_sym = 1.
wt_r_sm = 10 # smooth variation in r-direciton
wt_th_sm = 10 # smooth variation in th-direction

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

tol = 1e-8

data_dir = [
    '../data/m6_140km_1.00N_constant_1.00e-04mTBr_30k_210l/',
    '../data/m6_140km_1.00N_constant_1.00e-03mTBr_30k_210l/',
    '../data/m6_140km_1.00N_constant_0.01mTBr_30k_210l/',
    '../data/m6_140km_1.00N_constant_0.10mTBr_30k_210l/',
]

notify_me_by_text = True
verbose = False
num_threads = None