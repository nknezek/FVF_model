#! /usr/bin/env python3

import numpy as np

# model type
model_type = "MAC_model"

# mode to simulate (longitudinal)
m = [0]

# Size of grid
Nk = 10 # Radial cells
Nl = 20 # Latitudinal cells

# Define Physical Constants
R = 3480e3  # Outer core radius in (m)
h = [100.0e3, 139.2e3]  # layer thickness in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1.e4   # density in (kg/m^3)
nu = [0.8]   # momentum diffusivity in (m^2/s)
nu_th = 0.0
eta = 0.8  # magnetic diffusivity in (m^2/s)
eta_th = 0.0
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2
# dCyr = [2.5,5.,10.,20.,40.,60.,80.,160.,240.,320.,640.,960.,]
dCyr = [50.]

# background magnetic field (Tesla)
# choices: dipole, abs_dipole, constant, set
#   constant: must specify [Br, Bth, Bph] as floats
#   abs_dipole: must specify [Bd, Brnoise, Brconst, use_Bth, Bthnoise, Bthconst]
#   dipole: must specify [Bd, use_bth]
#   set: must specify [Br, Bth, Bph] as (Nk,Nl) arrays
B_type = 'abs_dipole'
Bd = 0.62e-3
Br = 0.
Brconst = 0.
Brnoise = 0.3e-3
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
boundary_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p')
dir_suf = 'general'
ep = 1e-3

notify_me_by_text = True
verbose = False
num_threads = None