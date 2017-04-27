#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [ -8, -10, -50, -500]

num_solutions_to_calculate = 300
num_solutions_to_plot = 40

# Parameters for deciding which solutions to keep or discard
filter_dict = {}
filter_dict['Q'] = {'minimum':0.01}
filter_dict['order_r'] = {'maximum':10, 'var':'uth'}
filter_dict['order_r'] = {'maximum':10, 'var':'uph'}
filter_dict['order_th'] = {'maximum':30, 'var':'uth'}
filter_dict['order_th'] = {'maximum':30, 'var':'uph'}

# Parameters for sorting solutions for plotting
sort_dict = {}
sort_dict['T'] = {'target': 1., 'weight':1.}
sort_dict['Q'] = {'target': 5., 'weight':1.}
sort_dict['order_r'] = {'target': 0., 'weight':1., 'var':'uth'}
sort_dict['order_th'] = {'target': 0., 'weight':1., 'var':'uth'}
sort_dict['region'] = {'target': 'equator', 'weight':1., 'var':'uth', 'split':0.5}
sort_dict['symmetry'] = {'target': 'symmetric', 'weight':1., 'var':'uth'}
sort_dict['smoothness_r'] = {'weight':10., 'var':'uth'}
sort_dict['smoothness_th'] = {'weight':10., 'var':'uth'}
sort_dict['power_in_layers'] = {'layer_wanted':1, 'weight':1., 'var':'ur', 'split_index':10}

# runtime parameters
notify_me_by_text = True
verbose = False
num_threads = None

# location of pre-computed data matrices
data_dir = [
    # '../data/m0_140km_1.00to4.00N_constant_0.60mTBr_70k_200l/',
    # '../data/m6_140km_1.00to4.00N_constant_0.60mTBr_70k_200l/',
    '../data/m0_200km_1.00to4.00N_constant_0.60mTBr_100k_200l/',
    '../data/m6_200km_1.00to4.00N_constant_0.60mTBr_100k_200l/',
]

# file names
filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

# numerical tolerance of solution
tol = 1e-6
