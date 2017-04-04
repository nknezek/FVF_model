#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [30., 100., 300., 900.]

nev = 30
num_to_keep = 10

eq_split = 0.5
target_r_order = 0
target_th_order = 4
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 0.
wt_r_order = 1.
wt_th_order = 1.
wt_region = 1.
wt_sym = 0.
wt_r_sm = 0.
wt_th_sm = 1e4

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

dCyr_list = [30., 100., 300., 900.]
data_dir = [
    '../data/k60_l120_m3_nu1e-01_60km_constantN100_constantBrB62_EMRthick/',
    '../data/k60_l120_m3_nu1e-01_100km_constantN100_constantBrB62_EMRthick/',
    '../data/k60_l120_m6_nu1e-01_60km_constantN100_constantBrB62_EMRthick/',
    '../data/k60_l120_m6_nu1e-01_100km_constantN100_constantBrB62_EMRthick/'
]


notify_me_by_text = True
verbose = False