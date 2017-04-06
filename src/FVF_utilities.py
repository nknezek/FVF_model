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

