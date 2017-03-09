#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [1.4/365.25]
delta_T = 20./365.25

nev = 70
d2_filter = 10e8
dth_filter = 10e8
eq_split = 1.0
eq_var = 'p'
real_var = 'uth'
filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'
use_initial_guess = False
oscillate = False
plot_robinson = False
plot_B_obs = False
plot_vel = True
zeros_wanted = list(range(10))
min_Q = 0.0
target_Q = 2.
tol = 1e-10


dCyr_list = [1.0]
data_dir = [
    '../data/k1_l80_m2_nu1e-01_100km_constantN0_constantBrB0_rossby/',
]



