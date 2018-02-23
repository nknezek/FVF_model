import numpy as np
import FVF_analysis as fana
import FVF_utilities as futil
import scipy.optimize as op
import scipy.misc as ms
import sys

pmagpath = '/Users/nknezek/code/pymagmods/pymagmods/'
sys.path.append(pmagpath)
import pymagmods as pmag
mflow = pmag.flows()
chaos6 = pmag.chaos6()


def fit_polynomial(x, y, deg=20):
    pfit = np.polyfit(x, y, deg)
    return lambda x: np.polyval(pfit, x)

def hermite_fun(x, l):
    c = np.zeros(40)
    c[l] = 1.
    return (2 ** l * ms.factorial(l) * np.pi ** 0.5) ** -0.5 * np.exp(-x ** 2 / 2) * np.polynomial.hermite.hermval(x, c)

def hermite_th(th, l, thd):
    x = 3 / thd * th
    if l >= 0:
        return hermite_fun(x, l)
    else:
        return np.zeros_like(th)

def hermite_sum(x, coeffs, delta_x):
    out = np.zeros_like(x)
    for l in range(len(coeffs)):
        out += coeffs[l] * hermite_fun(x / delta_x, l)
    return out

def hermite_fit_fun(x, fit_c):
    return hermite_sum(x, fit_c[:-1], fit_c[-1])

def fit_with_hermite(lat, data, deg, return_coeffs=False):
    fitfun_data = lambda c: np.sum((hermite_fit_fun(lat, c) - data) ** 2)
    c0 = np.ones((deg + 1))
    c0[-1] = 10.
    res = op.fmin_bfgs(fitfun_data, c0)
    outfun = lambda x: hermite_fit_fun(x, res)
    if return_coeffs:
        return outfun, res
    else:
        return outfun

def fit_horiz_flows_hermite(vec, Nk, Nl, deg, return_coefficients=False):
    Nk = 40
    Nl = 200
    dth = 180 / Nl
    lat = np.linspace(-90 + dth / 2, 90 - dth / 2, Nl)
    vec = fana.shift_vec_real_nomodel(vec, Nk, Nl, var='vph')
    vph = futil.get_variable_from_vec(vec, 'vph', Nk, Nl)
    vec = vec / np.max(np.abs(vph.real))
    vph = futil.get_variable_from_vec(vec, 'vph', Nk, Nl)
    vphr = vph[-1, :].real
    vphi = vph[-1, :].imag
    vth = futil.get_variable_from_vec(vec, 'vth', Nk, Nl)
    vthr = vth[-1, :].real
    vthi = vth[-1, :].imag
    ffphr, cphr = fit_with_hermite(lat, vphr, deg, return_coeffs=True)
    ffphi, cphi = fit_with_hermite(lat, vphi, deg, return_coeffs=True)
    ffthr, cthr = fit_with_hermite(lat, vthr, deg, return_coeffs=True)
    ffthi, cthi = fit_with_hermite(lat, vthi, deg, return_coeffs=True)
    if return_coefficients:
        return (ffthr, ffthi, ffphr, ffphi), (cthr, cthi, cphr, cphi)
    else:
        return (ffthr, ffthi, ffphr, ffphi)

def horiz_flows_given_coeffs(lat, coeffs, delta_th_override=None):
    if delta_th_override is None:
        vthr = hermite_fit_fun(lat, coeffs[0])
        vthi = hermite_fit_fun(lat, coeffs[1])
        vphr = hermite_fit_fun(lat, coeffs[2])
        vphi = hermite_fit_fun(lat, coeffs[3])
    elif delta_th_override is not None:
        vthr = hermite_sum(lat, coeffs[0][:-1], delta_th_override)
        vthi = hermite_sum(lat, coeffs[1][:-1], delta_th_override)
        vphr = hermite_sum(lat, coeffs[2][:-1], delta_th_override)
        vphi = hermite_sum(lat, coeffs[3][:-1], delta_th_override)
    return vthr + vthi * 1j, vphr + vphi * 1j

def horiz_vel_accel(l, m, delta_th, c012, lat, ph, t=2010, period=7.5, peak_flow=2., phase=0.):
    ''' compute the horizontal velocity and acceleration of a wave, given its hermite fit constants

    :param l: latitudinal wavenumber of wave
    :param m: longitudinal wavenumber of wave
    :param delta_th: latitudinal width of wave
    :param c012: hermite fit constants
    :param lat: latitude in [degrees]
    :param ph: longitude in [degrees]
    :param t: time in [years]
    :param period: period of wave in [years]
    :param peak_flow: maximum flow velocity in [km/yr]
    :return: vth, vph, ath, aph
    '''
    lon_grid, lat_grid = np.meshgrid(ph,lat)
    vth1, vph1 = horiz_flows_given_coeffs(lat, c012[l], delta_th_override=delta_th)
    v_magnitude = peak_flow/np.max(np.abs((vth1**2+vph1**2)**0.5))
    w = 2*np.pi/period*np.sign(m)
    m = np.abs(m)
    deg2rad = np.pi / 180
    phase_offset = np.mod(phase*deg2rad-w*2000,2*np.pi)
    vphi = (vph1[None,:]*np.exp(1j*(m*lon_grid.T*deg2rad + w*t + phase_offset))).T*v_magnitude
    vthi = (vth1[None,:]*np.exp(1j*(m*lon_grid.T*deg2rad + w*t + phase_offset))).T*v_magnitude
    vph = np.real(vphi)
    vth = np.real(vthi)
    aph = np.real(vphi*1j*w)
    ath = np.real(vthi*1j*w)
    return vth, vph, ath, aph

def u_v_divv(l, m, delta_th, c012, lat, ph, t=2010, period=7.5, peak_flow=2., phase=0.):
    ''' computes the longitudinal (u) and latitudinal (v) flows and divergence at a grid of points specifed by lat and ph

    :param l:
    :param m:
    :param delta_th:
    :param c012:
    :param lat:
    :param ph:
    :param t:
    :param period:
    :param peak_flow: maximum flow velocity in [km/yr]
    :return:
    '''
    lon_grid, lat_grid = np.meshgrid(ph,lat)
    vth, vph = horiz_flows_given_coeffs(lat, c012[l], delta_th_override=delta_th)
    v_magnitude = peak_flow/np.max(np.abs((vth**2+vph**2)**0.5))
    w = 2*np.pi/period
    u = np.real(vph[:,None]*np.exp(1j*(m*lon_grid*np.pi/180 + w*t + phase))).T*v_magnitude
    v = np.real(vth[:,None]*np.exp(1j*(m*lon_grid*np.pi/180 + w*t + phase))).T*v_magnitude
    absv = (u**2 + v**2)**.05
    dudph = np.gradient(u,axis=0)
    dvdth = np.gradient(v,axis=1)
    divv = dudph+dvdth
    return u,v,divv

def SV_from_hermite_flows(l, B_lmax, c012, t=2010, period=7.5, delta_th=17, peak_flow=2, Nth=200, v_lmax=14, m=6):
    dth = 180/Nth
    lat = np.linspace(-90+dth/2,90-dth/2,Nth)
    lon = np.linspace(-180+dth/2, 180-dth/2, Nth*2)
    sh = chaos6.get_shtcoeffs_at_t(t, l_max = 14)
    u,v,divv = u_v_divv(l,m, delta_th, c012, lat, lon, t=t, period=period, peak_flow=peak_flow)
    uSH = mflow.v2vSH(u.T)
    vSH = mflow.v2vSH(v.T)
    SV = mflow.SV_from_flow(vSH, uSH, sh, B_lmax=B_lmax, v_lmax=v_lmax, Nth=Nth)
    return SV

def vel_accel_wave_allT(l_wave, m, T, delta_th_wave, c012, Nth, phase, period, vmax=1):
    ''' computed velocity and acceleration of the wave in units of km/yr and km/yr^2

    :param l_wave:
    :param m:
    :param T:
    :param delta_th_wave:
    :param c012:
    :param Nth:
    :param phase:
    :param period:
    :param vmax:
    :return:
    '''
    dlat = 180/Nth
    lat = np.linspace(-90+dlat/2, 90-dlat/2, Nth)
    ph = np.linspace(dlat/2, 360-dlat/2, Nth*2)
    vth_t = np.empty((len(T), len(lat), len(ph)))
    vph_t = np.empty((len(T), len(lat), len(ph)))
    ath_t = np.empty((len(T), len(lat), len(ph)))
    aph_t = np.empty((len(T), len(lat), len(ph)))
    for i,t in enumerate(T):
        vth_t[i,:,:],vph_t[i,:,:],ath_t[i,:,:],aph_t[i,:,:] = horiz_vel_accel(l_wave, m, delta_th_wave, c012, lat, ph, t=t, period=period, phase=phase)
    vmag = np.max(np.sqrt(vth_t**2+vph_t**2))
    vth_t = vth_t*vmax/vmag
    vph_t = vph_t * vmax / vmag
    ath_t = ath_t*vmax/vmag
    aph_t = aph_t * vmax / vmag
    return vth_t, vph_t, ath_t, aph_t

def SV_wave_allT(Br_t, dthB_t, dphB_t, vth_t, vph_t, divv_t):
    ''' computes the secular variaion from wave motion given the data

    :param Br_t:
    :param dthB_t:
    :param dphB_t:
    :param vth_t:
    :param vph_t:
    :param divv_t:
    :return:
    '''
    return -Br_t * divv_t - dthB_t * vth_t - dphB_t * vph_t

def SA_wave_fluidaccel_allT(Br_t, dthB_t, dphB_t, ath_t, aph_t, diva_t):
    ''' computes secular acceleration from fluid acceleration given data

    :param Br_t:
    :param dthB_t:
    :param dphB_t:
    :param ath_t:
    :param aph_t:
    :param diva_t:
    :return:
    '''
    return -Br_t*diva_t - dthB_t*ath_t - dphB_t*aph_t

def SA_wave_magSV_allT(SV_t, dthSV_t, dphSV_t, vth_t, vph_t, divv_t):
    ''' computes seccular acceleration from SV and fluid velocity given data

    :param SV_t:
    :param dthSV_t:
    :param dphSV_t:
    :param vth_t:
    :param vph_t:
    :param divv_t:
    :return:
    '''
    return -SV_t*divv_t - dthSV_t*vth_t - dphSV_t*vph_t

def div_allT(v_th_t, v_ph_t, Nth=None, l_max=None):
    '''computes the divergence given the velocity (or acceleration) of a vector field at a given set of points

    :param v_th_t:
    :param v_ph_t:
    :param Nth:
    :param l_max:
    :return:
    '''
    dth_vth_t,_ = mflow.gradients_v_allT(v_th_t, Nth=Nth, l_max=l_max)
    _,dph_vph_t = mflow.gradients_v_allT(v_ph_t, Nth=Nth, l_max=l_max)
    return dth_vth_t+dph_vph_t

def compute_SASVwave_allT(phase, period, l_wave, m, T, delta_th_wave, c012, Nth,
                          B=None, dthB=None, dphB=None, SV=None, dthSV=None, dphSV=None):
    ''' computes the secular acceleration and secular variation over a series of times given all the data required

    :param phase:
    :param period:
    :param l_wave:
    :param m:
    :param T:
    :param delta_th_wave:
    :param c012:
    :param Nth:
    :param v_lmax:
    :param B:
    :param dthB:
    :param dphB:
    :param SV:
    :param dthSV:
    :param dphSV:
    :return:
    '''
    if (B is None) or (dthB is None) or (dthB is None):
        Bsh = chaos6.get_sht_allT(T)
        if B is None:
            B = chaos6.B_sht_allT(Bsh)
        if dthB is None or dphB is None:
            _, dthB, dphB = chaos6.gradB_sht_allT(Bsh)
    if (SV is None) or (dthSV is None) or (dphSV is None):
        SVsh = chaos6.get_SVsht_allT(T)
        if SV is None:
            SV = chaos6.B_sht_allT(SVsh)
        if dthSV is None or dphSV is None:
            _, dthSV, dphSV = chaos6.gradB_sht_allT(SVsh)
    vth_t, vph_t, ath_t, aph_t = vel_accel_wave_allT(l_wave, m, T, delta_th_wave, c012, Nth, phase, period)
    divv_t = div_allT(vth_t, vph_t, Nth)
    diva_t = div_allT(ath_t, aph_t, Nth)

    SAwave_t = SA_wave_fluidaccel_allT(B, dthB, dphB, ath_t, aph_t, diva_t)
    SAwave_t += SA_wave_magSV_allT(SV, dthSV, dphSV, vth_t, vph_t, divv_t)
    SVwave_t = SV_wave_allT(B, dthB, dphB, vth_t, vph_t, divv_t)
    return SAwave_t, SVwave_t

def make_SASV_from_phaseperiod_wave_function(l_wave, m, T, delta_th_wave, c012, Nth,
                                             B=None, dthB=None, dphB=None, SV=None, dthSV=None, dphSV=None):
    return lambda phase, period: compute_SASVwave_allT(phase, period, l_wave, m, T, delta_th_wave, c012, Nth,
                                                       B=B, dthB=dthB, dphB=dphB, SV=SV, dthSV=dthSV, dphSV=dphSV)




