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

# background magnetic field in (Tesla)
# chocies: dipole, dipoleBr, absDipole, absDipoleBr, constantBr, set, sinfuncBr
B_type = 'constantBr'

B_mag = [0.6e-3]
Bd = 0.0
Br = B_mag
Bth = 0.0
Bph = 0.0
const = 0.0
Bmax = 0.0
Bmin = 0.0
sin_exp = 0.0
Bnoise = 0.0

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
