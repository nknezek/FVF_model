#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [1/365.25, 10/365.25, 50/365.25, 1., -1/365.25, -10/365.25, -50/365.25, -1.]
# T_list = [3125, 10000, 30000, 90000, -3125, -10000, -30000, -90000]

nev = 250
num_to_keep = 75

eq_split = 0.5
target_r_order = 1
target_th_order = 1
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 1.
wt_r_order = 1e1
wt_th_order = 1e1
wt_region = 1e2
wt_sym = 1e1
wt_r_sm = 1. # smooth variation in r-direciton
wt_th_sm = 1. # smooth variation in th-direction

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
    '../data/m4_140km_1.00N_set_1e-02m2s_20k_200l/',
    # '../data/m6_140km_1.00N_set_1e+02m2s_20k_200l/',
    # '../data/m6_140km_1.00N_set_1e-02m2s_20k_200l/',
]

notify_me_by_text = True
verbose = False
num_threads = None