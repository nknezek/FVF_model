#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [-13., 13, -60., 60., -331, 331]

nev = 200
num_to_keep = 50

eq_split = 0.5
target_r_order = 2
target_th_order = 2
target_Q = 5.
target_symmetric = True
target_region='equator'

wt_T = 1.
wt_Q = 1.
wt_r_order = 1e1
wt_th_order = 1.
wt_region = 1.
wt_sym = 1.
wt_r_sm = 1e2 # smooth variation in r-direciton
wt_th_sm = 1e2 # smooth variation in th-direction

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


# dCyr_list = [7., 15., 30., 65., 150., 300., 600.]
dCyr_list = [35., 75., 150., 600.]
data_dir = [
    '../data/k20_l200_m6_nu1e-02_50km_constantN100_constantBrB62_EMR/',
    '../data/k20_l200_m6_nu1e-02_75km_constantN100_constantBrB62_EMR/',
    '../data/k20_l200_m6_nu1e-02_100km_constantN100_constantBrB62_EMR/',
    '../data/k20_l200_m6_nu1e-02_140km_constantN100_constantBrB62_EMR/',
]

notify_me_by_text = True
verbose = False
