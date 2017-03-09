#! /usr/bin/env python
import numpy as np
import macmodel_rlorentz as mac
import macloglib as mlog
import cPickle as pkl

# longitudinal mode to simulate
m = [6]

# Size of grid
Nk = 30 # Radial cells
Nl = 120 # Latitudinal cells

# Parameters that are usually constant
R = 3480e3  # Outer core radius in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1.e4   # density in (kg/m^3)
eta = 0.8  # magnetic diffusivity in (m^2/s)
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2

# Parameters that are usually varied
h = [140e3]  # layer thickness in (m)
dCyr = [9.]
nu = [1e-2]   # momentum diffusivity in (m^2/s)

# background magnetic field in (Tesla)
# chocies: dipole, dipoleBr, absDipole, absDipoleBr, constantBr, set, sinfuncBr
B_type = 'constantBr'
B_mag = [0.8e-3]
Bd = 0.0
Br = B_mag
Bth = 0.0
const = 0.0
Bmax = 0.0
Bmin = 0.0
sin_exp = 0.0
Bnoise = 0.0

# background velocity field in (m/s)
Uphi = np.ones((Nk+2, Nl+2))*0.0

# Buoyancy Frequency
# choices: constant, linear
buoy_type = 'constant'
buoy_ratio =  [4., 6. ]

# model parameters
model_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p', 'r_disp')
boundary_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p')
dir_suf = ''

# epsilon for numerical stability
ep = 1e-3
