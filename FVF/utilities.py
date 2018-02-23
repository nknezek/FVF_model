import os
import shutil
import sys
import importlib
import dill
from . import analyze as fana
import numpy as np
import pandas as pd
import sys, os, re


def convert_model_freq_to_period_yrs(omega):
    return 2*np.pi/1j*omega

def get_directory_name(param_dict, dir_suf=None):
    c = param_dict
    folder_name = '../data/m{0:.0f}_{1:.0f}km_{2}'.format(c['m'], c['h'] * 1e-3, c['buoyancy_type'])
    if type(c['N']) is np.ndarray:
        folder_name += '{:.2f}to{:.2f}N'.format(np.min(c['N']), np.max(c['N']))
    else:
        folder_name += '{:.2f}N'.format(c['N'])
    folder_name += '_{}'.format(c['B_type'])
    if c['Br'] > 0:
        if (c['Br'] * 1e3) >= 0.01:
            folder_name += '{0:.2f}mTBr'.format(c['Br'] * 1e3)
        else:
            folder_name += '{0:.2e}mTBr'.format(c['Br'] * 1e3)
    if c['Bth'] > 0:
        folder_name += '{0:.2f}mTBth'.format(c['Br'] * 1e3)
    if c['Bd'] > 0.:
        folder_name += '{0:.2f}mTBd'.format(c['Bd'] * 1e3)
    if c['Brnoise'] > 0.0:
        folder_name += '{0:.2f}mTBrnse'.format(c['Brnoise'] * 1e3)
    if c['Brconst'] > 0.0:
        folder_name += '{0:.2f}mTBrcst'.format(c['Brconst'] * 1e3)
    if (c['Brmult'] != 1.) and (c['Brmult'] != 0.):
        folder_name += '{0:.2f}mTBrmult'.format(c['Brmult'])
    if c['Bthnoise'] > 0.0:
        folder_name += '{0:.2f}mTBthnse'.format(c['Bthnoise'] * 1e3)
    if c['Bthconst'] > 0.0:
        folder_name += '{0:.2f}mTBthcst'.format(c['Bthconst'] * 1e3)
    if (c['Bthmult'] != 1.) and (c['Bthmult'] != 0.):
        folder_name += '{0:.2f}mTBthmult'.format(c['Bthmult'])
    if (type(c['nu']) is list):
        if len(c['nu']) > 1:
            folder_name += '_{0:.0e}m2snu'.format(c['nu'])
    elif c['nu'] !=1e-2:
        folder_name += '_{0:.2e}m2snu'.format(c['nu'])
    if (type(c['eta']) is list):
        if len(c['eta']) > 1:
            folder_name += '_{0:.2e}m2seta'.format(c['eta'])
    elif c['eta'] != 0.8:
        folder_name += '_{0:.2e}m2seta'.format(c['eta'])
    folder_name += '_{0}k_{1}l'.format(c['Nk'], c['Nl'])
    if dir_suf:
        folder_name += '_'
        folder_name += dir_suf
    folder_name += '/'
    return folder_name

def get_out_dir(out_dir_base, data_dir, num_data_dirs, T, num_T):
    out_dir = out_dir_base
    if num_data_dirs > 1:
        subfolders = ((data_dir.split('/'))[2]).split('_')
        for subfolder in subfolders[:-3]:
            out_dir += subfolder+'/'
    if num_T > 1:
        if T>0:
            direction = 'w'
        else:
            direction = 'e'
        out_dir +='{2}{0:03.0f}yrs{1:03.0f}days/'.format(np.abs(T)//1, (np.abs(T)%1)*365.25, direction)
    return out_dir

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        try:
            os.makedirs(d)
        except:
            pass

def store_config_file(config_file, out_dir_base):
    shutil.copyfile('../config/' + config_file + '.py', out_dir_base + config_file + '.py')

def find_available_skin_depths(directory):
    dC_list = []
    for file in os.listdir(directory):
        if file.startswith("A"):
            if file.endswith(".dat"):
                dC_list.append(float(file[1:-4]))
    return dC_list

# function to find nearest skin-depth value for a wave period
def find_closest_CC(Target, dCyr_list):
    dist = [abs(x-abs(Target)) for x in dCyr_list]
    return dCyr_list[dist.index(min(dist))]

def load_stored_data(filename):
    with open(filename, 'rb') as f:
        d1 = dill.load(f)
    return d1['model'], d1['vals'], d1['vecs']

def make_lk_dict(model, vals, vecs, var='vth', cutoff=0.075):
    vdict = {}
    for val,vec in zip(vals,vecs):
        vec = fana.shift_vec_real(model, vec, var=var)
        l = fana.get_order_th(model, vec, var=var, cutoff=cutoff)
        k = fana.get_order_r(model, vec, var=var, cutoff=cutoff)
        if k not in vdict.keys():
            vdict[k] = {}
        if l not in vdict[k].keys():
            vdict[k][l] = []
        vdict[k][l].append((val,vec))
    return vdict

def get_information_from_directory(filepath):
    '''
    gets information on data run from directory name

    :param filepath:
    :return:
    '''
    brstr = filepath
    hre = re.match(".*/(\d+)km.*", brstr)
    mre = re.match("/m(\d).*", brstr)
    Nre = re.match(".*[^\d\.]([\d\.]+)N.*", brstr)
    Bdre = re.match(".*[^\d\.]([\d\.]+)mTBd.*", brstr)
    Bnre = re.match(".*[^\d\.]([\d\.]+)mTBrnse.*", brstr)
    if hre:
        h = float(hre.group(1))
    if mre:
        m = int(mre.group(1))
    if Nre:
        N = float(Nre.group(1))
    if Bdre:
        Bd = float(Bdre.group(1))
    if Bnre:
        Br = float(Bnre.group(1))
    return m, h, N, Bd, Br

def load_lkdict_into_datad(lkdict, datad, model, filepath, kmax=4, lmax=10, save_vec=False, sm_max=0.6, verbose=False):
    '''
    loads wave information from an lkdictionary into a data dictionary for wave

    :param lkdict:
    :param datad:
    :param model:
    :param filepath:
    :param kmax:
    :param lmax:
    :param save_vec:
    :param sm_max:
    :param verbose:
    :return:
    '''
    N_saved = 0
    for k in lkdict.keys():
        if k <= kmax:
            for l in lkdict[k].keys():
                if l <= lmax:
                    for item in lkdict[k][l]:
                        val = item[0]
                        vec = item[1]
                        vecth = fana.shift_vec_real(model, vec, var='vth')
                        vecph = fana.shift_vec_real(model, vec, var='vph')
                        smth = fana.misfit_smoothness_th(model, vecth, {'var': 'vth'})
                        smph = fana.misfit_smoothness_th(model, vecph, {'var': 'vph'})
                        if smth < sm_max and smph < sm_max:
                            N_saved += 1
                            period = fana.get_period(model, val)
                            q = fana.get_Q(val)
                            if period < 0:
                                direction = 'east'
                            if period > 0:
                                direction = 'west'
                            # print('k={},l={},T {:.2f},Q {:.2f},{},{:.4f}'.format(k,l,np.abs(period),q,direction,smth))
                            Bnre = re.match(".*[^\d\.]([\d\.]+)mTBrnse.*", filepath)
                            Brn = float(Bnre.group(1))
                            datad['h'].append(model.h / 1e3)
                            datad['N'].append(np.mean(model.N))
                            datad['m'].append(model.m)
                            datad['Bd'].append(model.Bd * 1e3)
                            datad['Brnse'].append(Brn)
                            datad['l'].append(l)
                            datad['k'].append(k)
                            datad['period'].append(np.abs(period))
                            datad['Q'].append(q)
                            datad['direction'].append(direction)
                            if save_vec:
                                datad['vec'].append(vec)
    print('{} solutions saved'.format(N_saved))

def make_datad_empty(save_vec=False):
    '''
    makes an empty dictionary of wave data

    :param save_vec:
    :return:
    '''
    datad = {}
    datad['h'] = []
    datad['N'] = []
    datad['m'] = []
    datad['Bd'] = []
    datad['Brnse'] = []
    datad['l'] = []
    datad['k'] = []
    datad['period'] = []
    datad['Q'] = []
    datad['direction'] = []
    if save_vec:
        datad['vec'] = []
    return datad

def get_data_from_run(root_dir, verbose=False, kmax=1, lmax=10, save_vec=False, return_model=False):
    '''search root directory for all data.p output files and stores outputs into giant data-dictionary
    '''
    datad = make_datad_empty(save_vec=save_vec)
    filepaths = []
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            hasdata = re.match("data.p", file)
            if hasdata:
                filepath = root + '/' + file
                filepaths.append(filepath)
    if verbose:
        print('found {} data files'.format(len(filepaths)))
    for i, filepath in enumerate(filepaths):
        N = len(filepaths)
        subfilepath = filepath[len(root_dir):]
        if verbose:
            print('getting data from {}/{}, {}'.format(i, N, subfilepath))
        data = dill.load(open(filepath, 'rb'))
        lkdict = make_lk_dict(data['model'], data['vals'], data['vecs'])
        model = data['model']
        load_lkdict_into_datad(lkdict, datad, model, filepath, kmax=kmax, lmax=lmax, save_vec=save_vec, verbose=verbose)
    if return_model:
        return datad, model
    else:
        return datad

def datadict_to_dataframe(datadict, model, save_vec=False):
    '''
    converts dictionary of data into a pandas dataframe

    :param datadict:
    :param save_vec:
    :return:
    '''
    df = pd.DataFrame(datadict)
    # filter duplicates by converting periods and quality factors into strings
    df['periodstr'] = df['period'].apply(lambda x: '{:.1e}'.format(x))
    df['Qstr'] = df['Q'].apply(lambda x: '{:.1e}'.format(x))
    df['max_loc'], df['mean_loc'], df['bot_loc'], df['top_loc'] = zip(*df['vec'].map(lambda x: find_flow_loc(model,x)))
    df['region'] = df['mean_loc'].apply(get_wave_region)
    dfcl = df.drop_duplicates(['h', 'N', 'Brnse', 'direction', 'k', 'l', 'periodstr', 'Qstr', 'mean_loc'])
    if save_vec:
        col_order = ['h', 'N', 'direction', 'k', 'l', 'm', 'Bd', 'Brnse', 'period', 'Q', 'periodstr', 'Qstr', 'max_loc',
                     'mean_loc', 'bot_loc', 'top_loc', 'region','vec']
    else:
        col_order = ['h', 'N', 'direction', 'k', 'l', 'm', 'Bd', 'Brnse', 'period', 'Q', 'periodstr', 'Qstr', 'max_loc',
                     'mean_loc', 'bot_loc', 'top_loc']
    dfcl = dfcl[col_order]
    return dfcl

def select_from_df(df, N=None, h=None, direction=None, l=None, periodmax=None, periodmin=None, k=None, Brnse=None,
                   Bd=None,
                   mean_loc=None, bot_loc=None, top_loc=None, max_loc=None,
                   long_period=None, region=None, l_max=None):
    '''
    selects a particular wave solution from the dataframe, given various parameter choices

    :param df:
    :param N:
    :param h:
    :param direction:
    :param l:
    :param periodmax:
    :param periodmin:
    :param k:
    :param Brnse:
    :param Bd:
    :param mean_loc:
    :param bot_loc:
    :param top_loc:
    :param max_loc:
    :param long_period:
    :return:
    '''
    td = df
    if N is not None:
        td = td[np.isclose(td['N'],N)]
    if h is not None:
        td = td[np.isclose(td['h'],h)]
    if direction is not None:
        td = td[td['direction'] == direction]
    if l is not None:
        td = td[np.isclose(td['l'],l)]
    if Bd is not None:
        td = td[np.isclose(td['Bd'],Bd)]
    if Brnse is not None:
        td = td[np.isclose(td['Brnse'], Brnse)]
    if periodmax is not None:
        td = td[td['period'] < periodmax]
    if periodmin is not None:
        td = td[td['period'] > periodmin]
    if k is not None:
        td = td[np.isclose(td['k'], k)]
    if mean_loc is not None:
        td = td[np.isclose(td['mean_loc'], mean_loc)]
    if bot_loc is not None:
        td = td[np.isclose(td['bot_loc'], bot_loc)]
    if top_loc is not None:
        td = td[np.isclose(td['top_loc'], top_loc)]
    if max_loc is not None:
        td = td[np.isclose(td['max_loc'], max_loc)]
    if long_period is not None:
        td = td[td['long_period'] == long_period]
    if region is not None:
        td = td[td['region'] == region]
    if l_max is not None:
        td = td[td['l'] <= l_max+0.5]
    return td

def select_from_df_keyvalue(df, k,v):
    '''
    selects a particular wave solution from the dataframe, given various parameter choices

    :param df:
    :param N:
    :param h:
    :param direction:
    :param l:
    :param periodmax:
    :param periodmin:
    :param k:
    :param Brnse:
    :param Bd:
    :param mean_loc:
    :param bot_loc:
    :param top_loc:
    :param max_loc:
    :param long_period:
    :return:
    '''
    td = df
    if k is 'N':
        td = select_from_df(df,N=v)
    elif k is 'h':
        td = select_from_df(df, h=v)
    elif k is 'direction':
        td = select_from_df(df, direction=v)
    elif k is 'l':
        td = select_from_df(df, l=v)
    elif k is 'periodmax':
        td = select_from_df(df, periodmax=v)
    elif k is 'periodmin':
        td = select_from_df(df, periodmin=v)
    elif k is 'k':
        td = select_from_df(df, k=v)
    elif k is 'Brnse':
        td = select_from_df(df, Brnse=v)
    elif k is 'Bd':
        td = select_from_df(df, Bd=v)
    elif k is 'mean_loc':
        td = select_from_df(df, mean_loc=v)
    elif k is 'bot_loc':
        td = select_from_df(df, bot_loc=v)
    elif k is 'top_loc':
        td = select_from_df(df, top_loc=v)
    elif k is 'max_loc':
        td = select_from_df(df, max_loc=v)
    elif k is 'region':
        td = select_from_df(df, region=v)
    return td

def drop_entry(df, N=None, h=None, direction=None, l=None, periodmax=None, periodmin=None, k=None, Brnse=None, Bd=None,
               mean_loc=None, bot_loc=None, top_loc=None, max_loc=None):
    '''
    drops entries from the dataframe matching particular parameter choices

    :param df:
    :param N:
    :param h:
    :param direction:
    :param l:
    :param periodmax:
    :param periodmin:
    :param k:
    :param Brnse:
    :param Bd:
    :param mean_loc:
    :param bot_loc:
    :param top_loc:
    :param max_loc:
    :return:
    '''
    td = select_from_df(df, N=N, h=h, direction=direction, l=l, periodmax=periodmax, periodmin=periodmin, k=k,
                        Brnse=Brnse, Bd=Bd, mean_loc=mean_loc, bot_loc=bot_loc, top_loc=top_loc, max_loc=max_loc)
    df.drop(td.index, inplace=True)

def find_flow_loc(model, vec, minflow=0.25):
    '''

    :param vec:
    :param minflow:
    :return:
        ma: maximum flow location at CMB in degrees co-latitude
        me: mean flow location between bottom and top of flow in degrees co-latitude
        b: bottom flow location at CMB in degrees co-latitude
        t: top flow location at CMB in degrees co-latitude
    '''
    vth = model.get_variable(vec, 'vth')
    vph = model.get_variable(vec, 'vph')
    vmag = np.sqrt(np.abs(vth) ** 2 + np.abs(vph) ** 2)
    vmth = np.mean(vmag, axis=0)
    vmax = np.max(vmth)
    th = np.linspace(-90, 90, len(vmth))
    maxind = np.where(vmax == vmth)[0][0]
    ma = th[maxind]
    try:
        botind = np.where(vmth[:maxind] < vmax * minflow)[0][-1]
        topind = maxind + np.where(vmth[maxind:] < vmax * minflow)[0][0]
        meanind = int(np.mean([botind, topind]))
        t = th[topind]
        b = th[botind]
        me = th[meanind]
    except:
        me = ma
        b = ma
        t = ma
    if me < 0:
        t *= -1
        b *= -1
        ma *= -1
        me *= -1
    if b > t:
        tmp = b
        b = t
        t = tmp
    return ma, me, b, t

def get_wave_region(mean_loc, eq_th=20):
    '''
    get flow region (equator or mid-latitude)

    :param df:
    :param eq_th:
    :return:
    '''
    if mean_loc < eq_th:
        return 'equator'
    else:
        return 'mid-latitude'

def get_branch(dfg, periodmin):
    if dfg['period'] > periodmin:
        return 'long'
    else:
        return 'short'

def get_variable_from_vec(vector, var, Nk, Nl, model_variables = ('vr', 'vth', 'vph', 'bth', 'bph', 'p', 'ur')):
        '''
        Takes a flat vector and a variable name, returns the variable in a
        np.matrix
        inputs:
            vector: flat vector array with len == SizeM
            var: str of variable name in model
        outputs:
            variable in np.array
        '''
        if (var not in model_variables):
            raise RuntimeError('variable not in model_variables')
        elif len(vector) != Nk*Nl*len(model_variables):
            raise RuntimeError('vector given is not correct length in this \
                               model')
        else:
            var_start = get_index(0, 0, var, Nk, Nl, model_variables=model_variables)
            var_end = get_index(Nk-1, Nl-1, var, Nk, Nl, model_variables=model_variables)+1
            variable = np.array(np.reshape(vector[var_start:var_end], (Nk, Nl), 'F'))
            return variable

def get_index(k, l, var, Nk, Nl, model_variables = ('vr', 'vth', 'vph', 'bth', 'bph', 'p', 'ur')):
    '''
    Takes coordinates for a point, gives back index in matrix.
    inputs:
        k: k grid value from 0 to K-1
        l: l grid value from 0 to L-1
        var: variable name in model_variables
    outputs:
        index of location in matrix
    '''
    Size_var = Nk*Nl

    if (var not in model_variables):
        raise RuntimeError('variable not in model_variables')
    elif not (l >= 0 and l <= Nl - 1):
        raise RuntimeError('l index out of bounds')
    elif not (k >= 0 and k <= Nk - 1):
        raise RuntimeError('k index out of bounds')
    return Size_var * model_variables.index(var) + k + l * Nk
