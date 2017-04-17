#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

# T_list = [25/365.25, 50/365.25, 100/365.25, 200/365.25, -25/365.25, -50/365.25, -100/365.25, -200/365.25]
T_list = [100/365.25, 200/365.25, 2., 4., 8., 16., 32., 64.]
nev = 200
num_to_keep = 20

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
    # '../data/m0_80km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m0_80km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m4_80km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m4_80km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m6_80km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m6_80km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m0_140km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m0_140km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m4_140km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m4_140km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m6_140km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m6_140km_1.00N_set_1e-02m2s_20k_200l/',
    '../data/m0_140km_1.00N_set_30k_210l/',
    # '../data/m1_140km_1.00N_set_30k_210l/',
    # '../data/m2_140km_1.00N_set_30k_210l/',
    # '../data/m3_140km_1.00N_set_30k_210l/',
    # '../data/m4_140km_1.00N_set_30k_210l/',
    # '../data/m5_140km_1.00N_set_30k_210l/',
    # '../data/m6_140km_1.00N_set_30k_210l/',
]

notify_me_by_text = True
verbose = False
num_threads = None