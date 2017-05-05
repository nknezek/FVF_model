#! /usr/bin/env python3
import numpy as np
""" Configuration File for SLEPc Run solving MAC model """

T_list = [ -8, 8., -10, 10, -50, 50., -500, 500.,]
target_Q = 5.

num_solutions_to_calculate = 200
num_solutions_to_plot = 50

# Parameters for deciding which solutions to keep or discard
filter_dict = {}
filter_dict['Q'] = {'minimum':0.01}
filter_dict['order_r'] = {'maximum':10, 'var':'vth'}
filter_dict['order_r'] = {'maximum':10, 'var':'vph'}
filter_dict['order_th'] = {'maximum':30, 'var':'vth'}
filter_dict['order_th'] = {'maximum':30, 'var':'vph'}

# Parameters for sorting solutions for plotting
sort_dict = {}
sort_dict['T'] = {'target': 1., 'weight':1.}
sort_dict['Q'] = {'target': 5., 'weight':1.}
sort_dict['order_r'] = {'target': 0., 'weight':1., 'var':'vth'}
sort_dict['order_th'] = {'target': 0., 'weight':1., 'var':'vth'}
sort_dict['region'] = {'target': 'equator', 'weight':1., 'var':'vth', 'split':0.5}
sort_dict['symmetry'] = {'target': 'symmetric', 'weight':1., 'var':'vth'}
sort_dict['smoothness_r'] = {'weight':10., 'var':'vth'}
sort_dict['smoothness_th'] = {'weight':10., 'var':'vth'}
sort_dict['power_in_layers'] = {'layer_wanted':1, 'weight':1., 'var':'ur', 'split_index':10}

# runtime parameters
notify_me_by_text = True
verbose = False
num_threads = None

# location of pre-computed data matrices
data_dir = [
    '../data/m0_140km_constant_1.00N_constant_0.62mTBr_51k_150l/',
    '../data/m6_140km_constant_1.00N_constant_0.62mTBr_51k_150l/',
]

# file names
filemodel = 'model.p'
fileA = 'A'
fileB = 'B'
savefile = 'data.p'

# numerical tolerance of solution
tol = 1e-6
