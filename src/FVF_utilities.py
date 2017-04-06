import os
import shutil


def get_directory_name(param_dict, include_nu=False):
    c = param_dict
    folder_name = '../data/m{0:.0f}_{1:.0f}km_{2:.2f}N_{3}'.format(c['m'], c['h'] * 1e-3, c['buoy_ratio'], c['B_type'])
    if c['Br'] > 0:
        folder_name += '_{0:.2f}mTBr'.format(c['Br'] * 1e3)
    if c['Bd'] > 0.:
        folder_name += '_{0:.2f}mTBd'.format(c['Bd'] * 1e3)
    if c['noise'] > 0.0:
        folder_name += '_{0:.2f}mTnse'.format(c['noise'] * 1e3)
    if c['const'] > 0.0:
        folder_name += '_{0:.2f}mTcst'.format(c['const'] * 1e3)
    if include_nu:
        folder_name += '_{0:.0e}m2s'.format(c['nu'])
    folder_name += '_{0}k_{1}l/'.format(c['Nk'], c['Nl'])
    return folder_name


def get_out_dir(out_dir_base, data_dir, num_data_dirs, T, num_T):
    out_dir = out_dir_base
    if num_data_dirs > 1:
        subfolders = data_dir.strip('/').split('_')
        for subfolder in subfolders[:-2]:
            out_dir += subfolder+'/'
    if num_T > 1:
        out_dir +='{0:.2f}yrs/'.format(T)
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
