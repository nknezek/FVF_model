import logging
import os
import FVF_utilities as futil

def setup_custom_logger(dir_name='./',filename='MAC.log', verbose=True):
    futil.ensure_dir(dir_name)
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler = logging.StreamHandler()
    fileHandler = logging.FileHandler(dir_name+filename)

    handler.setFormatter(formatter)
    fileHandler.setFormatter(formatter)

    logger = logging.getLogger(dir_name)
    logger.setLevel(logging.DEBUG)
    if verbose:
        logger.addHandler(handler)
    logger.addHandler(fileHandler)
    return logger

def log_model(logger,model):
    m = model.m_values[0]
    Nk = model.Nk
    Nl = model.Nl
    R = model.R
    Omega = model.Omega
    h = model.h
    rho = model.rho
    nu = model.nu
    eta = model.eta
    Bd = model.Bd
    B0 = model.B0
    dr = model.dr
    dth = model.dth
    t_star = model.t_star
    r_star = model.r_star
    u_star = model.u_star
    P_star = model.P_star
    B_star = model.B_star
    omega_g = model.omega_g
    G = model.G
    E = model.E
    Prm = model.Prm
    l1 = model.l1

    logger.info(
    '\nMAC model, m={0}, Nk={1}, Nl={2}\n'.format(m,Nk,Nl)
    +'\nPhysical Parameters\n'
    +'Omega = {0:.2e} rad/s\n'.format(Omega)
    +'R = {0} km\n'.format(R*1e-3)
    +'h = {0} km\n'.format(h*1e-3)
    +'rho = {0:.2e} kg/m^3\n'.format(rho)
    +'nu = {0:.2e} m/s^2\n'.format(nu)
    +'eta = {0:.2e} m/s^2\n'.format(eta)
    +'omega_g = {0:.2e} rad/s = {1:.1f}*Omega\n'.format(omega_g,omega_g/Omega)
    +'Bd = {0:.2e} T, {1:.2e} G\n'.format(Bd, Bd*1e4)
    +'\nNon-Dimensional Parameters\n'
    +'t_star = {0:.2e} s\n'.format(t_star)
    +'r_star = {0:.2e} m\n'.format(r_star)
    +'u_star = {0:.2e} m/s\n'.format(u_star)
    +'B_star = {0:.2e} T = {1:.2f} G\n'.format(B_star,B_star*1e4)
    +'P_star = {0:.2e} Pa = {1:.2f} GPa\n'.format(P_star,P_star*1e-9)
    +'\nG = {0:.2e}\n'.format(G[0,0])
    +'E = {0:.2e}\n'.format(E)
    +'E/Prm = {0:.2e}\n'.format(E/Prm)
    +'Bdnon-dim) = {0:.1e}\n'.format(Bd/B_star)
    +'B0 at equator = {0:.1e}\n'.format(B0[Nk/2,Nl/2+1])
    +'B0 at pole = {0:.1e}\n'.format(B0[1,1])
    +'\nGrid Spacing Evaluation\n'
    +'dr^2/dth^2 = {0:.3e}\n'.format(dr**2/dth**2)
    +'E/dr^2 = {0:.1e}\n'.format(E/dr**2)
    +'E/dth^2 = {0:.1e}\n'.format(E/dth**2)
    +'\nMomentum Equation Terms\n'
    +'lambda*u ~ l1 = {0:.2e}\n'.format(l1)
    +'drP ~ 1/dr = {0:.2e}\n'.format(1/dr)
    +'E/Prm*drB0b ~ E/Prm*Bd/dr = {0:.2e}\n'.format(E/Prm*Bd/dr)
    +'E dr^2 u ~ E/dr^2 = {0:.2e}\n'.format(E/dr**2)
    +'2cos(th)*u ~ 2 = {0:.2e}\n'.format(2)
    +'Gdelr ~ G*l1 = {0:.2e}\n'.format(G[0,0]*l1)
    +'\nLorentz Equation Terms\n'
    +'lambda*b ~ l1 = {0:.2e}\n'.format(l1)
    +'dr uxB ~ Bd/dr = {0:.2e}\n'.format(u_star*Bd/dr)
    +'E/Prm/dr^2 ~ {0:.2e}\n'.format(E/Prm/dr**2)
    )

