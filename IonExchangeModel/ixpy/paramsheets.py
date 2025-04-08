# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:10:32 2020

@author: BDATSOV
"""

import pandas as pd
import numpy as np
import sys

pd.options.mode.chained_assignment = None 

def lowerEntries(item):    
    '''
        Retruns: a lowercase string iff the item is a string;
        
        Parameters: strings or numbers;
    '''        
            
    if isinstance(item, str):
        lower_item = map(lambda x: x.lower(), item)
        s = ''
        lower_item = s.join(list(lower_item))
    else:
        lower_item = item
    return lower_item

def lowercaseData(data):
    '''
        Returns: lowercase pandas dataframe;
        
        Parameters: pandas dataframe;
    '''
    
    lower_df = data.applymap(lowerEntries)
    return lower_df

def lowercaseIndex(data):
    '''
        Retruns: pandas dataframe with index in lowercase;
            
        Parameters: pandas dataframe;
    '''        
            
    lower_index = data.rename(str.lower, axis='index')
    return lower_index

   
def conv_units(u_in, u_out, u_list, u_coef, label, caller):
    '''
        Returns:
            a conversion factor based on "units in" and "units out";
            "units in" as string;
            "units out" as string;
          
        Parameters:
            u_in : string;
            u_out : string;
            u_list : list;
            u_coeff : list;
            label : string;
            caller : string;
         
        *** label and caller are optional parameters used for system messaging ***
        *** enter empty strings if label or caller are None ***
    '''    
    
    if u_in in u_list:
        u_in_pos = u_list.index(u_in)
    else:
 
        print('WARNING: ' + caller + ' is not in proper units!')
        print('Acceptable units for ' + label + ' are: ', *u_list, sep='\n')
        
        pos_assigned = True
        while pos_assigned:
            try:
                u_in = str(input('Please eneter the correct units here:'))
                if u_in in u_list:
                    u_in_pos = u_list.index(u_in)
                    pos_assigned = False
                    break
            except ValueError:
                print('Entered units are not accepted.')
                
    if u_out in u_list:
        u_out_pos = u_list.index(u_out)
    else:
        print(u_out + ' are not acceptable units.')
        print('Enter units in one of the following units list:', *u_list, sep='\n')
        u_out_pos = False
        sys.exit()
    
    base_vector = np.array(u_coef)
    base_data = np.array([base_vector]*len(base_vector))
    l_factor = base_data/base_vector[:,None]
    factor = l_factor[u_in_pos,u_out_pos]
    
    return factor, u_in, u_out



def conv_length(u_in, u_out, label, caller):
    units = ['m','cm','mm','ft','in', 'um']
    coefs = [1.,100.,1000.,3.2808,39.3701, 1.0e6]
    cv, u_in, u_out = conv_units(u_in,u_out, units, coefs, label, caller)
    
    return cv, u_in, u_out


def conv_volume(u_in, u_out, label, caller):
    units = ['m3','cm3','ml','L','mm3','gal','ft3','in3','cuf', 'l']
    coefs = [1., 1.0e6, 1.0e6, 1.0e3, 1.0e9, 264.172, 35.3147, 61023.801579099, 35.3147, 1.0e3]
    cv, u_in, u_out = conv_units(u_in,u_out, units, coefs, label, caller)
    
    return cv, u_in, u_out


def conv_time(u_in, u_out, label, caller):
    units = ['days','hours','min','sec','m','s', 'seconds', 'h','d', 'hr', 'day', 'hour']
    coefs = [1. , 24., 24.*60., 24*60*60, 24.*60., 24*60*60, 24*60*60, 24., 1., 24., 1., 24.]
    cv, u_in, u_out = conv_units(u_in,u_out, units, coefs, label, caller)
    
    return cv, u_in, u_out


def conv_vol_per_time(u_in, u_out, label, caller):
    vars_lst = [u_in, u_out]
    vol = []
    time = []
    for var in vars_lst:
        if '/' in var:
            vol.append(var.split('/')[0])
            time.append(var.split('/')[1])
        elif 'gpm' in var:
            vol.append('gal')
            time.append('min')
        else:
            print('else')
    v_f, u_in_v, u_out_v = conv_volume(vol[0],vol[1], label, caller)
    t_f, u_in_t, u_out_t = conv_time(time[0], time[1], label, caller)
    
    if u_in_v == 'gal':
        u_in_v = 'g'
        u_in = u_in_v + 'p' + u_in_t
    else:
        u_in = u_in_v + '/' + u_in_t        
    
    cv = v_f/t_f, u_in, u_out
    
    return cv


def conv_vel(u_in, u_out, label, caller):
    vars_lst = [u_in, u_out]
    lenght = []
    time = []
    for var in vars_lst:
        if '/' in var:
            lenght.append(var.split('/')[0])
            time.append(var.split('/')[1])
        elif 'gpm' in var:
            lenght.append('gal')
            time.append('min')
        else:
            print('else')
    l_f, u_in_l, u_out_l = conv_length(lenght[0],lenght[1], label, caller)
    t_f, u_in_t, u_out_t = conv_time(time[0], time[1], label, caller)
    
    u_in = u_in_l + '/' + u_in_t
    
    cv = l_f/t_f
    
    return cv, u_in, u_out

def conv_dens(u_in, u_out, label, caller):
    vars_lst = [u_in, u_out]
    wght = []
    vol = []
    for var in vars_lst:
        if '/' in var:
            wght.append(var.split('/')[0])
            vol.append(var.split('/')[1])
        elif 'gpm' in var:
            wght.append('gal')
            vol.append('min')
        else:
            print('else')
    w_f, u_in_w, u_out_w = conv_weight(wght[0],wght[1], label, caller)
    v_f, u_in_v, u_out_v = conv_volume(vol[0], vol[1], label, caller)
    
    if u_in_v == 'gal':
        u_in_v = 'g'
        u_in = u_in_w + 'p' + u_in_v
    else:
        u_in = u_in_w + '/' + u_in_v
    
    cv = w_f/v_f
    
    return cv, u_in, u_out

def conv_area(u_in, u_out, label, caller):
    units = ['m2', 'cm2', 'mm2', 'in2', 'ft2']
    coefs = [1., 10000., 1000000., 1550.0031, 10.7639]
    cv, u_in, u_out = conv_units(u_in, u_out, units, coefs, label, caller)
    
    return cv, u_in, u_out
    

def conv_area_per_time(u_in, u_out, label, caller):
    vars_lst = [u_in, u_out]
    area = []
    time = []
    for var in vars_lst:
        if '/' in var:
            area.append(var.split('/')[0])
            time.append(var.split('/')[1])
        elif 'gpm' in var:
            area.append('gal')
            time.append('min')
        else:
            print('else')
    l_f, u_in_l, u_out_l = conv_area(area[0],area[1], label, caller)
    t_f, u_in_t, u_out_t = conv_time(time[0], time[1], label, caller)
    
    u_in = u_in_l + '/' + u_in_t
    
    cv = l_f/t_f
    
    return cv, u_in, u_out
    
    
    
def conv_fr_to_vel(u_in, u_out, diam, diam_u, label, caller):

    if '/' in u_out:
        l_out, time_out = u_out.split('/')
    else:
        l_out = 'gal'
        time_out = 'min'
        
    tmp_vol_u = diam_u + '3' + '/' + time_out       # create diam_u per time_out
    vol_to_d, u_in_vd, u_out_vd = conv_vol_per_time(u_in, tmp_vol_u, label, caller)  #conv vol units to diam units
    tmp_cv, u_in_tmp, u_out_tmp = conv_length(diam_u, l_out, label, caller)           # conv linear units to cm
    cv = vol_to_d*tmp_cv/(np.pi*diam**2/4)
    u_in = u_in_vd + '/' + u_out_vd
    
    return cv, u_in, u_out


def conv_conc(u_in, u_out, label, caller, **kwargs):
    MW = kwargs.get('MW', 1.)
    val = kwargs.get('val', 1.)
    units = ['meq', 'mg', 'ug', 'ng', 'mgN', 'mgC', 'eq', 'g']
    coefs = [1, MW/val, 1000.*MW/val, 1.0e6*MW/val, 14.001/val, 12.011/val, 1/1000, 1/1000*MW/val]
    cv, u_in, u_out = conv_units(u_in, u_out, units, coefs, label, caller)
    return cv, u_in, u_out

def conv_weight(u_in, u_out, label, caller):
    units = ['kg', 'g', 'mg', 'ug', 'ng', 'lbs', 'lb', 'lbm', 'oz']
    coefs = [1., 1000., 1.0e6, 1.0e9, 1.0e12, 2.2046226218, 2.2046226218, 2.2046226218, 35.27396195]
    cv, u_in, u_out = conv_units(u_in, u_out, units, coefs, label, caller)
    return cv, u_in, u_out

def conv_capacity_mass(u_in, u_out, label, caller, **kwargs):
    vars_lst = [u_in, u_out]
    conc = []
    wght = []
    for var in vars_lst:
        if '/' in var:
            conc.append(var.split('/')[0])
            wght.append(var.split('/')[1])
        elif 'gpm' in var:
            conc.append('gal')
            wght.append('min')
        else:
            print('else')
    c_f, u_in_c, u_out_c = conv_conc(conc[0],conc[1], label, caller, **kwargs)
    w_f, u_in_w, u_out_w = conv_weight(wght[0], wght[1], label, caller)
    
    cv = c_f/w_f
    u_in = u_in_c + '/' + u_in_w
    
    return cv, u_in, u_out

def conv_capacity_vol(u_in, u_out, label, caller, **kwargs):
    vars_lst = [u_in, u_out]
    conc = []
    vol = []
    for var in vars_lst:
        if '/' in var:
            conc.append(var.split('/')[0])
            vol.append(var.split('/')[1])
        elif 'gpm' in var:
            conc.append('gal')
            vol.append('min')
        else:
            print('else')
    c_f, u_in_c, u_out_c = conv_conc(conc[0],conc[1], label, caller, **kwargs)
    v_f, u_in_v, u_out_v = conv_volume(vol[0], vol[1], label, caller)
    
    cv = c_f/v_f
    u_in = u_in_c + '/' + u_in_v
    
    return cv, u_in, u_out

cap_factor = conv_capacity_mass('meq/g','meq/kg', '','')

    
def conv_database(data_in, u_in, u_out, conv_fn, MW, val):
    '''
        Retruns: 
            pandas dataframe in requested units
        
        Parameters:
            data_in : pandas dataframe
            u_in : dictionary (units in);
            u_out : dictionary (unist out);
            conv_fn : dictionary (conversion functions);
            MW : dictionary (molecular weights);
            val : dictionary (valences);
            
            ***the keys must be the same for all dictionaries***
    '''
    for u in u_in.keys():
        args = {'MW':MW[u], 'val':val[u]}
        cf, u_in[u], u_out[u] = conv_fn(u_in[u], u_out[u], u, u, **args)
        data_in[u] = cf * data_in[u]
        
    return data_in, u_in, u_out

     

def conv_params(data_in, u_in, u_out, conv_fn):
    '''
        Retruns: 
            pandas dataframe in requested units
            
        Parameters:
            data_in : pandas dataframe ()
            u_in : dictionary (units in);
            u_out : dictionary (unist out);
            conv_fn : dictionary (conversion functions);
    
            ***the keys must be the same for all dictionaries***
    '''        
    for u in u_in.keys():
        tmp_conv_fn = conv_fn[u]
        cf, u_in[u], u_out[u] = tmp_conv_fn(u_in[u], u_out[u], u, u)
        try:
            data_in.loc[u, 'value'] *= cf
        except:
            # Compatability for Shiny model input files
            data_in.loc[u, 'value'] = cf # Calculate cf correctly and replace old value
        data_in.loc[u, 'units'] = u_out[u]        
        
    return data_in

def conv_the_same(u_in, u_out, label, caller):
    cv = 1.
    return cv, u_in, u_out

def Q_calc(**kwargs):
    '''
    Returns:
        Q (capacity by volume)
    Parameters:
        **kwargs : keyword argiments as a dictionary
    
    *** takes Qm and RHOP, or Qf and EBED to calculate Q ***
    '''
    
    if 'qm' in kwargs and 'rhop' in kwargs:
        Qm = kwargs.get('qm', None)
        RHOP = kwargs.get('rhop', None)
        Q = Qm * RHOP
    elif 'qf' in kwargs and 'ebed' in kwargs:
        Qf = kwargs.get('qf', None)
        EBED = kwargs.get('ebed', None)
        Q = Qf/(1-EBED)
    else:
        print('WARNING: Missing parameters! \n' \
              'Resin capacity by volume cannot been calculated.')
    
    return Q


def dict_from_data(data, val):
    '''
    Returns:
        dictionary
    Paremeters:
        data : pandas dataframe
        val : dataframe column name as a string
    '''
    data_dict = {}
    for index, row in data.iterrows():
        data_dict[index] = data.loc[index, val]
        
    return data_dict


def conv_params_data(data):
    '''    
        Returns:
            pandas dataframe in units acceptable by HSDMIX.solver();    
        
        Parameters:
            data : pandas dataframe
    
            *** takes parameter units from "units" column in the input file ***
            *** pass the input units to conv_data() ***
            *** takes the velocity if both velocity and flow rate are present;
                otherwise converts the flow rate to linear velocity ***         
    '''
    low_idx = lowercaseIndex(data)
    low_data = lowercaseData(low_idx)
    
    
    params_Dct = low_data.to_dict('index')

    u_in = {}
    
    u_out = {'time':'s', 'rhop':'g/ml', 'rb':'cm', 'kl':'cm/s', \
              'ds':'cm2/s', 'v':'cm/s', 'qm':'meq/kg', 'l':'cm', \
             'flrt':'cm3/s', 'diam':'cm', 'qf':'meq/L', 'ebed':None, \
             'nz':None, 'nr':None, 'epor':None, 'dp':'cm2/s', \
                'w':None, 'm':'g', 'vl':'L', 't_end':'s'}
    
    correct_idx = {'rhop':'RHOP', 'qm':'Qm', 'ebed':'EBED', 'l':'L', 'kl':'kL',\
                   'ds':'Ds', 'qf':'Qf', 'dp':'Dp', 'epor':'EPOR', 'vl':'VL', \
                    'vr':'VR','t_end':'t_end'}
    
    u_fn = {'time':conv_time, 'rhop':conv_dens, 'rb':conv_length, 'kl':conv_vel, \
             'ds':conv_area_per_time, 'v':conv_vel, 'qm':conv_capacity_mass, \
             'flrt':conv_vol_per_time, 'l':conv_length, 'diam':conv_length, \
             'qf':conv_capacity_vol, 'ebed':conv_the_same, 'nz':conv_the_same, \
             'nr':conv_the_same, 'epor':conv_the_same, 'dp':conv_area_per_time, \
                'w':conv_the_same, 'm':conv_weight, 'vl':conv_volume, \
                    't_end':conv_time}
    
    
    for p in params_Dct.keys():
        u_in[p] = params_Dct[p]['units']
        u_out[p] = u_out[p]
        u_fn[p] = u_fn[p]

    
    params_out = conv_params(low_data, u_in, u_out, u_fn)
    
    ####################################
    # check for velocity in parameters #
    ####################################
    if 'v' in params_out.index or 'flrt' in params_out.index:

        if 'v' not in params_out.index:
            flrt_f, u_in_fr, u_out_fr = conv_fr_to_vel(low_data.loc['flrt','units'], 'cm/s', \
                                    low_data.loc['diam','value'], \
                                    low_data.loc['diam','units'], '', '')    
        
            flrt_cv = params_out.loc['flrt','value']*flrt_f
            
            v_row = pd.DataFrame([[flrt_cv, 'cm/s']], columns = ['value', 'units'], \
                                index = ['v'])
                
            params_out = pd.concat([params_out, v_row], axis=0)
        else:
            params_out = params_out
            if 'flrt' in params_out.index:
                print('The linear velocity has been used instead of the flow rate.')
    else:
        '''
        If velocity of flow rate are not provided assumed is a batch model simulation
        Calculate VR and add it to the params dataframe
        '''
        vr = (1-params_out.loc['w','value'])*params_out.loc['m','value']/\
            params_out.loc['rhop','value']/1.0E3
        vr_row = pd.DataFrame([[vr, 'L']], columns = ['value', 'units'], \
                                index = ['VR'])
        params_out = pd.concat([params_out, vr_row], axis=0)
    
    
    Q_dict = dict_from_data(params_out, 'value')

    Q = Q_calc(**Q_dict)
    
    q_row = pd.DataFrame([[Q, 'meq/L']], columns = ['value', 'units'], \
                             index = ['Q'])

    params_out = pd.concat([params_out, q_row], axis=0)
    
       
    for cp in correct_idx.keys():
        params_out.rename(index={cp:correct_idx[cp]}, inplace=True)
 
    return params_out


