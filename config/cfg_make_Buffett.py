#! /usr/bin/env python3

import numpy as np
import macmodel_bdiv as mac
import macloglib as mlog
import pickle as pkl

# mode to simulate (longitudinal)
m = [0]

# Size of grid
Nk = 20 # Radial cells
Nl = 200 # Latitudinal cells

# Define Physical Constants
R = 3480e3  # Outer core radius in (m)
h = [139.2e3]  # layer thickness in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1.e4   # density in (kg/m^3)
nu = [0.8]   # momentum diffusivity in (m^2/s)
nu_th = 0.0
eta = 0.8  # magnetic diffusivity in (m^2/s)
eta_th = 0.0
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2
dCyr = [65.]

# background magnetic field (Tesla)
# choices: dipole, abs_dipole, constant, set
#   constant: must specify [Br, Bth, Bph] as floats
#   abs_dipole: must specify [Bd, Brnoise, Brconst, use_Bth, Bthnoise, Bthconst]
#   dipole: must specify [Bd, use_bth]
#   set: must specify [Br, Bth, Bph] as (Nk,Nl) arrays
B_type = 'constant'
Bd = 0.
Br = 0.62e-3
Brconst = 0.
Brnoise = 0.
Bth = 0.
Bthconst = 0.
Bthnoise = 0.
Bph = 0.
use_Bth = False

# background velocity field in (m/s)
Uphi = 0.0

# Buoyancy Frequency
# choices: constant, linear
buoy_type = 'constant'
buoy_ratio =  [1.0]

# model parameters
model_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p', 'r_disp')
dir_suf = '_Buffett'
ep = 1e-3

notify_me_by_text = True
verbose = False
num_threads = None