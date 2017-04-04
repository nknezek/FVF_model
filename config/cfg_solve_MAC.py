#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [64., 151., 181.]

nev = 10
num_to_keep = 10

eq_split = 0.5
target_r_order = 0
target_th_order = 4
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 1.
wt_r_order = 1.
wt_th_order = 1.
wt_region = 1.
wt_sym = 1.
wt_r_sm = 1. # smooth variation in r-direciton
wt_th_sm = 1e2 # smooth variation in th-direction

eq_var = 'uph'
r_ord_var = 'uph'
th_ord_var = 'uph'
real_var = 'uph'

filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

oscillate = False
plot_robinson = False
plot_B_obs = False
plot_vel = True

tol = 1e-8


dCyr_list = [65., 150., 300., 600., 1200.]
data_dir = [
    '../data/k20_l200_m0_nu1e-02_135km_constantN80_constantBrB62_MAC/',
]


notify_me_by_text = True
verbose = False