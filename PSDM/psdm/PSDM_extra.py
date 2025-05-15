# -*- coding: utf-8 -*-
"""
Created on Mon May 10 08:02:30 2021

@author: JBurkhar
"""

import pandas as pd
import numpy as np

def calc_mesh_diameter(mesh='12x40'):
    large, small = mesh.split('x')

    sizes = {'4': 4.76, '5': 4.0, '6': 3.36, '7': 2.83, '8': 2.38, '10': 2.0,
             '12': 1.68, '14': 1.41, '16': 1.19, '18': 1.0,'20': 0.841,'25': 0.707,
             '30': 0.595,'35': 0.5,'40': 0.420, '45':0.354, '50':0.297,'60':0.25,
             '70':0.21, '80':0.177,'100':0.149,'120':0.125,'140':0.105,'170':0.088,
             '200':0.074, '230':0.063,'270':0.053,'325':0.044,'400':0.037}

    try:
        large_diam = sizes[large]
    except:
        large_diam = 0.
        print('Error: Upper mesh size not in available sizes')

    try:
        small_diam = sizes[small]
    except:
        small_diam = 0.
        print('Error: Lower mesh size not in available sizes')


    average_diam = (large_diam + small_diam)/2.0

    return average_diam


def _length_convert2cm(value, unit):
    unit = unit.lower()
    cm_val = np.nan
    if unit == 'cm':
        cm_val = value
    elif unit == 'm':
        cm_val = value * 100.
    elif unit == 'mm':
        cm_val = value / 10.
    elif unit == 'in':
        cm_val = value * 2.54
    elif unit == 'ft':
        cm_val = value * 12 * 2.54
    elif unit == 'yard' or unit == 'yd':
        cm_val = value * 3 * 12 * 2.54

    return cm_val


def _flrt_convert2cm3min(value, unit):
    unit = unit.lower()
    cm_val = np.nan
    if unit == 'cm3/min' or unit == 'ml/min':
        cm_val = value
    elif unit == 'm3/min':
        cm_val = value * (100.)**3
    elif unit == 'mm3/min':
        cm_val = value / (10.)**3
    elif unit == 'in3/min':
        cm_val = value * (2.54)**3
    elif unit == 'ft3/min':
        cm_val = value * (12 * 2.54)**3
    elif unit == 'mgd':
        cm_val = value * 3785.411784 * 1e6/1440.
    elif unit == 'lpm' or unit == 'l/min':
        cm_val = value * 1000.
    elif unit == 'gal/min' or unit == 'gpm':
        cm_val = value * 3785.411784
    elif unit == 'm3/day' or unit == 'm3pd':
        cm_val = value * 1e6/1440.

    return cm_val


def _mass_convert2g(value, unit):
    unit = unit.lower()
    g_val = np.nan

    if unit == 'gm' or unit == 'g' or unit == 'grams':
        g_val = value
    elif unit == 'kg':
        g_val = value * 1e3
    elif unit == 'lb' or unit == 'pounds':
        g_val = value * 453.59237
    elif unit == 'ton' or unit == 'tons':
        g_val = value * 2000. * 453.59237
    elif unit == 'mton' or unit == 'metric tons':
        g_val = value * 1000. * 1e3

    return g_val


def _time_convert2day(value, unit):
    unit = unit.lower()
    t_val = np.nan

    if unit == 'min':
        t_val = value / (24. * 60.)
    elif unit == 'sec' or unit == 's':
        t_val = value / (24 * 60 * 60.)
    elif unit == 'day' or unit == 'days' or unit == 'd':
        t_val = value
    elif unit == 'hr' or unit == 'hours' or unit == 'h':
        t_val = value / 24.

    return t_val


def build_large_scale_column(column_info=None, L=None, diam=None, flrt=None,
                             loading_rate=None, bulk_density=None,
                             adsorb_mass=None, operation_time=None,
                             particle_diameter=None, mesh=None):

    lc_obj = pd.DataFrame(columns=['value'], index=['flrt', 'L', 'diam','part_diam',
                                                    'operation time', 'mass',
                                                    'EBCT', 'density', 'loading rate'])

    if column_info == None:
        if L != None:
            if type(L) == list and len(L) == 2:
                value, unit = L
                lc_obj.loc['L','value'] = _length_convert2cm(value, unit)
            else:
                # assume provided in cm already
                lc_obj.loc['L','value'] = L

        if diam != None:
            if type(diam) == list and len(diam) == 2:
                value, unit = diam
                lc_obj.loc['diam','value'] = _length_convert2cm(value, unit)
            else:
                # assume provided in cm already
                lc_obj.loc['diam','value'] = diam

        if flrt != None:
            if type(flrt) == list and len(flrt) == 2:
                value, unit = flrt
                lc_obj.loc['flrt','value'] = _flrt_convert2cm3min(value, unit)
            else:
                # assume provided in cm**3/min already
                lc_obj.loc['flrt','value'] = flrt

        if flrt == None and loading_rate != None:
            if type(loading_rate) == list and len(loading_rate) == 2:
                value, unit = loading_rate

                unit = unit.lower()
                if unit == 'm/h' or unit == 'm/hr':
                    lc_obj.loc['loading rate', 'value'] = value
                elif unit == 'gpm/ft2':
                    lc_obj.loc['loading rate', 'value'] = value * 2.44475
                elif unit == 'cm/s':
                    print('here')
                    lc_obj.loc['loading rate', 'value'] = value * _time_convert2day(1,'hr') / _time_convert2day(1,'s') / _length_convert2cm(1,'m')

            else:
                # assumes this is provided in m/hr
                lc_obj.loc['loading rate', 'value'] = loading_rate


            lc_obj.loc['flrt','value'] = lc_obj['value']['loading rate'] * _length_convert2cm(1, 'm') * \
                                         (_time_convert2day(1,'min') /\
                                         _time_convert2day(1,'hr')) * np.pi/4.0 *\
                                         lc_obj['value']['diam']**2

        if particle_diameter != None:
            if type(particle_diameter) == list and len(particle_diameter) == 2:
                value, unit = particle_diameter
                lc_obj.loc['part_diam', 'value'] = _length_convert2cm(value, unit)
            else:
                lc_obj.loc['part_diam', 'value'] = particle_diameter
        elif mesh != None:
            if type(mesh) == str and len(mesh.split('x')) == 2:
                lc_obj.loc['part_diam', 'value'] = calc_mesh_diameter(mesh=mesh)


        if adsorb_mass != None:
            if type(adsorb_mass) == list and len(adsorb_mass) == 2:
                value, unit = adsorb_mass
                lc_obj.loc['mass','value'] = _mass_convert2g(value, unit)
            else:
                #assumes number already provided as grams
                lc_obj.loc['mass', 'value'] = adsorb_mass

        if operation_time != None:
            if type(operation_time) == list and len(operation_time) == 2:
                value, unit = operation_time
                lc_obj.loc['operation time', 'value'] = _time_convert2day(value, unit)
            else:
                #assumes number already provided in days
                lc_obj.loc['operation time', 'value'] = operation_time

        if bulk_density != None:
            if type(bulk_density) == list and len(bulk_density) == 2:
                value, unit = bulk_density
                unit = unit.lower()
                if unit == 'gm/ml' or unit == 'g/ml' or unit == 'gm/cm3' or unit == 'g/cm3':
                    lc_obj.loc['density', 'value'] = value
                if unit == 'lb/ft3':
                    lc_obj.loc['density', 'value'] = value / 62.4279605761

            else:
                #assumes g/cm**3
                lc_obj.loc['density', 'value'] = bulk_density
    else:
        print('functionality not yet included')


    lc_obj.loc['EBCT','value'] = lc_obj['value']['L'] * np.pi/4 * \
                                 (lc_obj['value']['diam']**2)/\
                                 lc_obj['value']['flrt']

    if np.isnan(lc_obj.loc['loading rate', 'value']):
        # assumes m/h
        lc_obj.loc['loading rate', 'value'] = lc_obj.loc['flrt','value'] / \
                                              (np.pi/4. * lc_obj.loc['diam','value']**2) /\
                                                   _length_convert2cm(1, 'm') /\
                                                   (_time_convert2day(1, 'min')/ _time_convert2day(1, 'hr'))

    return lc_obj


def calc_RSSCT_specs(lc_obj, column_diameter=1.,
                     particle_diameter=None, mesh=None, power_x=0):
    sc_obj = pd.DataFrame(columns=['value'], index=['diameter (cm)',
                                                    'length (cm)',
                                                    'particle diameter (mm)',
                                                    'EBCT (min)',
                                                    'duration (d)',
                                                    'loading rate (m/h)',
                                                    'Q (ml/min)',
                                                    'mass (g)',
                                                    'Volume Water Needed (L)'])

    if column_diameter != None:
        if type(column_diameter) == list and len(column_diameter) == 2:
            value, unit = column_diameter
            sc_obj.loc['diameter (cm)', 'value'] = _length_convert2cm(value, unit)
        else:
            sc_obj.loc['diameter (cm)', 'value'] = column_diameter

    if mesh != None and particle_diameter == None:
        part_diam = calc_mesh_diameter(mesh=mesh)
        sc_obj.loc['particle diameter (mm)'] = part_diam
        info_avail = True
    elif particle_diameter != None:
        if type(particle_diameter) == list and len(particle_diameter) == 2:
            value, unit = particle_diameter
            sc_obj.loc['particle diameter (mm)'] = _length_convert2cm(value, unit) *\
                                                   _length_convert2cm(1, 'mm')
        else:
            sc_obj.loc['particle diameter (mm)'] = particle_diameter
        info_avail = True
    else:
        info_avail = False
        print('Insufficient information. Must provide either particle_diameter or mesh')

    if info_avail:
        sc_obj.loc['EBCT (min)','value'] = lc_obj['value']['EBCT'] *\
                                   (sc_obj.loc['particle diameter (mm)', 'value']**(2-power_x))/\
                                   (lc_obj.loc['part_diam', 'value']**(2-power_x))

        sc_obj.loc['duration (d)', 'value'] = lc_obj['value']['operation time'] *\
                                   (sc_obj.loc['particle diameter (mm)', 'value']**(2-power_x))/\
                                   (lc_obj.loc['part_diam', 'value']**(2-power_x))

        sc_obj.loc['loading rate (m/h)', 'value'] = lc_obj['value']['loading rate'] *\
                                                   lc_obj['value']['part_diam']/\
                                                   sc_obj['value']['particle diameter (mm)']

        sc_obj.loc['length (cm)', 'value'] = sc_obj['value']['EBCT (min)'] *\
                                             sc_obj['value']['loading rate (m/h)'] *\
                                             _time_convert2day(1, 'min') / \
                                             _time_convert2day(1,'hr') *\
                                             _length_convert2cm(1, 'm')

        col_volume = sc_obj['value']['diameter (cm)']**2 * np.pi/4. * sc_obj['value']['length (cm)']
        sc_obj.loc['Q (ml/min)', 'value'] = col_volume / sc_obj['value']['EBCT (min)']

        sc_obj.loc['mass (g)', 'value'] = col_volume * lc_obj['value']['density']

        sc_obj.loc['Volume Water Needed (L)'] = sc_obj['value']['duration (d)'] *\
                                                sc_obj['value']['Q (ml/min)'] /\
                                                _time_convert2day(1, 'min') /\
                                                _flrt_convert2cm3min(1, 'lpm')

    return sc_obj


if __name__ == "__main__":

    # =============================================================================
    # Test
    # =============================================================================

    #replicating MWH - Crittenden's example 15-14
    c = calc_mesh_diameter()
    c = calc_mesh_diameter(mesh='140x200')

    lc = build_large_scale_column(L=[83.3, 'cm'],
                                  diam=[5.1, 'cm'],
                                  flrt=[170.1, 'ml/min'],
                                   # loading_rate=[.138889,'cm/s']
                                   # mesh='12x40',
                                  particle_diameter=1.0,
                                  adsorb_mass=[833.8, 'g'],
                                  operation_time=[100, 'day'],
                                  bulk_density=0.49,
                                  )
    print('Large Column Specs\n', 'x'*30)
    print(lc)

    sc = calc_RSSCT_specs(lc, particle_diameter=0.21, column_diameter=[1.10, 'cm'])
    # sc = calc_RSSCT_specs(lc, mesh='60x80', column_diameter=[.375, 'in'], power_x=.7)

    print('\nSmall Column Specs\n', 'x'*30)
    print(sc)

    expected_values = [1.1, 17.493, 0.21, 0.4411735945442268, 4.409999999999999,
       23.790635091937293, 37.68166089965399, 8.145835355873526,
       239.29361937716263]
    test_values = np.array(sc.value.values, dtype=np.float64)
    assert np.allclose(test_values, expected_values)


    #example 15-15
    print('\n Example #2 \n')
    lc = build_large_scale_column(L=[71.4, 'cm'],
                                  diam=[40.64, 'cm'],
                                  flrt=[18.9016, 'l/min'], #error in book. used volume (92618 cm3)/ebct (294 s)
                                  particle_diameter=1.62,
                                  adsorb_mass=[40751.92, 'g'], #volume * bulk density
                                  operation_time=[288, 'day'],
                                  bulk_density=0.44,
                                  )
    print('Large Column Specs\n', 'x'*30)
    print(lc)

    sc = calc_RSSCT_specs(lc, particle_diameter=0.212, column_diameter=[1.10, 'cm'])
    # sc = calc_RSSCT_specs(lc, mesh='60x80', column_diameter=[.375, 'in'], power_x=.7)

    print('\nSmall Column Specs\n', 'x'*30)
    print(sc)

    expected_values = [1.1, 9.343703703703705, 0.212, 0.0839148101156238,
        4.932126200274347, 66.8084955980663, 105.8170606465034,
        3.9070321625037034, 751.5404600409199]
    test_values = np.array(sc.value.values, dtype=np.float64)
    assert np.allclose(test_values, expected_values)

