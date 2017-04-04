#! /usr/bin/env python3
import numpy as np

# model type
model_type = "BTRossby_model"

# mode to simulate (longitudinal)
m = [3,6]

# Size of grid
Nk = 60 # Radial cells
Nl = 120 # Latitudinal cells

# Define Physical Constants
R = 3480e3  # Outer core radius in (m)
h = [60e3, 100e3]  # layer thickness in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1.e4   # density in (kg/m^3)
nu = [1e-1]   # momentum diffusivity in (m^2/s)
eta = 0.8  # magnetic diffusivity in (m^2/s)
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2
dCyr = [30., 100., 300., 900.]

# background magnetic field in (Tesla)
# chocies: dipole, dipoleBr, absDipole, absDipoleBr, constantBr, set, sinfuncBr
B_type = 'constantBr'

B_mag = [0.62e-3]
Bd = B_mag
Br = B_mag
Bth = None
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
boundary_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p')
dir_suf = '_EMRthick'
ep = 1e-6

notify_me_by_text = True