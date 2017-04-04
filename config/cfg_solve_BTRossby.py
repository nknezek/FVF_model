#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [5.5/365., 10.0/365., 15./365., 50./365.]

nev = 50
num_to_keep = 10

eq_split = 0.5
target_r_order = 0
target_th_order = 4
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 0.
wt_Q = 0.
wt_r_order = 0.
wt_th_order = 0.
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

tol = 1e-12

dCyr_list = [1.0]
data_dir = [
    # '../data/k1_l180_m1_nu1e-01_100km_constantN0_constantBrB0_rossby/'
    # '../data/k1_l400_m2_nu1e-01_100km_constantN0_constantBrB0_rossby/',
    # '../data/k1_l180_m3_nu1e-01_100km_constantN0_constantBrB0_rossby/'
    # '../data/k1_l180_m4_nu1e-01_100km_constantN0_constantBrB0_rossby/'
    # '../data/k1_l200_m2_nu1e-01_100km_constantN0_constantBrB0_rossby/',
    # '../data/k1_l200_m2_nu1e+02_100km_constantN0_constantBrB0_rossby/'
    '../data/k1_l800_m2_nu1e+02_100km_constantN0_constantBrB0_rossby/'
]




notify_me_by_text = True
verbose = False
