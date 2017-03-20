#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [20., 60., 80.1, 100.,]

nev = 10
num_to_keep = 10

eq_split = 0.5
target_r_order = 0
target_th_order = 0
target_Q = 5.
target_symmetric = True

wt_T = 1.
wt_Q = 1.
wt_r_order = 4.
wt_th_order = 2.
wt_region = 1.
wt_sym = 1.

eq_var = 'p'
r_ord_var = 'p'
th_ord_var = 'p'
real_var = 'ur'

filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

oscillate = False
plot_robinson = False
plot_B_obs = False
plot_vel = True

tol = 1e-8

dCyr_list = [68.44627]
data_dir = [
    '../data/k10_l20_m0_nu8e-01_139km_constantN100_absDipoleBrB31_general/',
]