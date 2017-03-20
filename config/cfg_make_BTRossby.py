#! /usr/bin/env python3
import numpy as np

# model type
model_type = "BTRossby_model"

# mode to simulate (longitudinal)
m = [1,2,3,4]

# Size of grid
Nk = 1 # Radial cells
Nl = 800 # Latitudinal cells

# Define Physical Constants
R = 3480e3  # Outer core radius in (m)
h = 100e3  # layer thickness in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1e4   # density in (kg/m^3)
nu = 1e2  # momentum diffusivity in (m^2/s)
eta = 1.0  # magnetic diffusivity in (m^2/s)
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2
dCyr = 1.0 # Period of wave for magnetic boundary condition (years)

# background magnetic field in (Tesla)
# chocies: dipole, dipoleBr, absDipole, absDipoleBr, constantBr, set, sinfuncBr
B_type = 'constantBr'
B_mag = 0.0
Bd = 0.0
Br = 0.0
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
buoy_ratio = 0.0

# model parameters
model_variables = ('uth', 'uph', 'p')

# Directory suffix name to save model
dir_suf = '_rossby'
ep = 1e-8

