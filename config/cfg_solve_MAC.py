#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [65.]
delta_T = 80.

nev = 10
d2_filter = 10e8
dth_filter = 10e8
eq_split = 1.0
eq_var = 'p'
real_var = 'ur'
filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'
use_initial_guess = False
oscillate = False
plot_robinson = False
plot_B_obs = False
plot_vel = True
zeros_wanted = range(10)
min_Q = 1.
target_Q = 5.
tol = 1e-8


dCyr_list = [65.]
data_dir = [
    '../data/k20_l200_m0_nu1e-02_135km_constantN80_constantBrB62_MAC/',
]



