#! /usr/bin/env python3
import numpy as np

# model type
model_type = "MAC_model"

# mode to simulate (longitudinal)
m = [0]

# layer thickness in (m)
h = [140e3]

# Size of grid
Nk = 70 # Radial cells
Nl = 200 # Latitudinal cells

# Buoyancy Frequency
# choices: constant, linear, set
buoy_type = 'set'
h_upper = 40e3
dr = np.array(h)/Nk
N_upper = 4
N_lower = 1
buoy_type = 'set'
buoy_ratio = np.ones((Nk,Nl))*N_lower
buoy_ratio[-int(h_upper/dr):] *= N_upper/N_lower

# background magnetic field (Tesla)
# choices: dipole, abs_dipole, constant, set
#   constant: must specify [Br, Bth, Bph] as floats
#   abs_dipole: must specify [Bd, Brnoise, Brconst, use_Bth, Bthnoise, Bthconst]
#   dipole: must specify [Bd, use_bth]
#   set: must specify [Br, Bth, Bph] as (Nk,Nl) arrays
B_type = 'constant'
Bd = 0.
Br = 0.6e-3
Brconst = 0.
Brnoise = 0.
Brmult = 1.
Bth = 0.
Bthconst = 0.
Bthnoise = 0.
Bthmult = 0.
Bph = 0.
use_Bth = False

dCyr = [60., 180., 500.]  # Period of wave for magnetic boundary condition (years)


# Define Physical Constants
R = 3480e3  # Outer core radius in (m)
Omega = 2*np.pi/(23.9345*3600.0)  # rotation rate in (rad/s)
rho = 1.e4   # density in (kg/m^3)
nu = [1e-2]   # momentum diffusivity in (m^2/s)
eta = [0.8]  # magnetic diffusivity in (m^2/s)
mu_0 = 4.*np.pi*10.**-7  # vacuum permeability in (kg*m/(A^2s^2))
g = 10.  # Gravity in m/s^2

# background velocity field in (m/s)
Uphi = 0.0

# model parameters
model_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p', 'r_disp')
boundary_variables = ('ur', 'uth', 'uph', 'br', 'bth', 'bph', 'p')
dir_suf = '2layer'
ep = 1e-4

notify_me_by_text = True
verbose = False
num_threads = None




