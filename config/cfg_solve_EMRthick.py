#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [30., 100., 300., 900.]
delta_T = 900.

nev = 200

eq_split = 0.5
target_zeros = 0
target_Q = 1.
target_symmetric = True

eq_var = 'p'
real_var = 'ur'

filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

use_initial_guess = False

oscillate = True
plot_robinson = False
plot_B_obs = False
plot_vel = True

tol = 1e-8




# dCyr_list = [120.]
# dCyr_list = [42., 110., 980.]
dCyr_list = [975.]
data_dir = [
    '../data/k60_l120_m6_nu1e-01_140km_constantN84_constantBrB62_EMRThick/',
]



