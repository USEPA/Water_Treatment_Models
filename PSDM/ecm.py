# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:00:15 2020

Equilibrium Theory Model
ECM

Originally coded by Thomas F. Speth, Mich. Tech. U.
Revised by Beth Erlanson, U. Texas at Austin

Recoded into Python
@author: UCChEJBB, Jonathan Burkhardt
"""

import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.optimize import minimize
from scipy.optimize import differential_evolution as de
from scipy.interpolate import interp1d

ml_pergal = 3785.411784
cm_perft = 12. * 2.54
sec_permin = 60.

class ecm():
    def __init__(self, filename, **kw):
        comp_sheet_name = kw.get('compound_sheet', 'compounds')
        bed_sheet_name = kw.get('column_sheet', 'bed')
        
        self.unittype = kw.get('units', 'molar')

        
        self.compdata, self.beddata = _process_input_file(filename, 
                                                          comp_sheet_name, 
                                                          bed_sheet_name,
                                                          values=self.unittype)
        
        
        self.n = len(self.compdata.columns)
    
    def run_ecm(self, print_results=False):
        zz = np.ones(self.n)
        bvf = np.ones(self.n)
        usage_rate = np.ones(self.n)
        self.vw = np.zeros(self.n)
        self.qstore = np.zeros((self.n, self.n))
        self.cstore = np.zeros((self.n, self.n))
        self.cstore[0] = self.compdata.loc['Co']
        dg = np.zeros((self.n, self.n))

        q = zz * self.compdata.loc['K'] * self.compdata.loc['Co'] ** self.compdata.loc['1/n']

        a = root(self.fcn, q.values)
        self.qstore[0] = a.x#*self.compdata.loc['MW'].values
        
        q11 = (a.x*self.compdata.loc['MW']).values[0]
        rhob = self.beddata['bed_density']
        C01 = self.compdata.loc['C0'].values[0]
        vf = self.beddata['vf']
        eps = self.beddata['void_i']
        self.vw[0] =  vf * C01 / (q11 * rhob + C01 * eps) # vf has eps in it, vf = vf * eps
        
        if 'L' in self.beddata.index:
            bed_length = self.beddata['L']
        else:
            bed_length = 2.765 # meters, assumed... may need to supply this as input
        
        bvf[0] = vf/self.vw[0] # vf has eps in it, dropped vf * eps
        
        dg[0] = rhob * a.x/(self.beddata['void_i'] * self.compdata.loc['Co'])
        self.sstc=(1000. * self.compdata.loc['Co']**(1-self.compdata.loc['1/n'])/\
               self.compdata.loc['K']).values        
        
        
        for i in np.arange(1, self.n):
            # print(i)
            q = zz * self.compdata.loc['K'] * \
                self.cstore[i-1] ** self.compdata.loc['1/n']
            q_vals = q.values[i:]
            b = root(self.fcn_extra, q_vals, method='lm')
            self.qstore[i][i:] = b.x
            self.cstore[i][i:] = self._fcn2_e(b.x)
            
            dg[i][i:] = rhob * b.x/(eps * self.cstore[i][i:])
            
              #wave velocity in parts
            vf_elem = self.beddata['vf'] * self.cstore[0][i]
            denom = (self.qstore[i][i] * rhob + self.cstore[i][i] * eps)
            mid1 = ((self.qstore.transpose()[i] * rhob + self.cstore.transpose()[i] * eps))[:i]
            vw_part = (self.vw - np.append([0],self.vw)[:self.n])[:i]
            mid_part = np.sum(mid1 * vw_part)
            third_part = denom * self.vw[i-1] #* eps
            self.vw[i] = (vf_elem - mid_part + third_part) / denom
            
            bvf[i] = vf/self.vw[i]
            
        usage_rate = rhob / bvf
        self.days_to_breakthrough = bed_length * 100./(self.vw * 3600 * 24)
        bvf[0] = vf/self.vw[0] # vf has eps in it, dropped vf * eps
        
        usage_rate = rhob / bvf
        self.dg = dg * 1.
        self.usage_rate = 1/usage_rate
        self.bvf = 1 * bvf        
        self.cratio = self.cstore/self.cstore[0]
        
        
        if print_results:
            print(self.compdata.columns.values)

            print('days to breakthrough')
            print(np.round(self.days_to_breakthrough,2))
 
            print('bvf')
            print(np.round(bvf,2))
        
            print('usage rate (m**3/kg)')
            print(np.round((1/usage_rate),3))
  
            print('single solute treatement capacity (mg Carbon/L water)')
            print(np.round(self.sstc,3)) # seems correct

            print('c store (ug/L)')
            print(np.round(self.cstore*self.compdata.loc['MW'].values,3))
            print('C/C0')
            print(np.round(self.cratio, 3))
        
        
        
    
    def _fcn2(self, x):
        qt = np.sum(x)
        qnq = np.sum((self.compdata.loc['xn'].values * x))
        
        ctmp = x/qt * (qnq/(self.compdata.loc['xn'].values *\
                            self.compdata.loc['K'].values))**self.compdata.loc['xn'].values
        return ctmp
        
    def fcn(self, x):
        ctmp = self._fcn2(x)
        
        return ctmp * self.compdata.loc['MW'].values - self.compdata.loc['C0'].values
    
    def _fcn2_e(self, x):
        start_idx = self.n - len(x)
        qt = np.sum(x)
        xn = self.compdata.loc['xn'].values[start_idx:]
        kval = self.compdata.loc['K'].values[start_idx:]
        qnq = np.sum((xn * x))
        
        ctmp = x/qt * (qnq/(xn * kval))**xn
        
        return ctmp
    
    def fcn_extra(self, x):
        start_idx = self.n - len(x)
        last_idx = start_idx - 1

        ctmp = self._fcn2_e(x)
        
        qI_last = self.qstore[last_idx][start_idx:]
        vf = self.beddata['vf']
        vw_last = self.vw[last_idx]
        eps = self.beddata['void_i']
        rhob = self.beddata['bed_density']
        c_last = self.cstore[last_idx][start_idx:]
        
        c_calc = (x - qI_last) * rhob * vw_last/((vf - vw_last*eps)) + c_last
                  
                 
        return c_calc - ctmp
    
    def fit_ecm(self, brk_days, peak_conc, xguess, xbounds, brk_weight=2):
        '''
        Attempts to backfit K & 1/n given days to breakthrough and peak C/C0

        Parameters
        ----------
        brk_days : TYPE
            DESCRIPTION.
        peak_conc : TYPE
            DESCRIPTION.

        Returns
        -------
        int
            DESCRIPTION.

        '''
        
        def _helper(array):
            ks = array[:self.n]
            xni = array[self.n:]
            self.compdata.loc['K'] = ks * 1.
            self.compdata.loc['1/n'] = xni * 1.
            self.compdata.loc['xn'] = 1./xni
            
            self.run_ecm()
            
            error1 = np.sum((self.days_to_breakthrough - brk_days)**2) * brk_weight 
            # adds weight to breakthrough day agreement
            error2 = np.sum((np.diag(self.cratio) - peak_conc)**2)

            
            tmp = error1 * error2 * 1e5

            return tmp

        fitter = minimize(_helper, xguess, method='L-BFGS-B', bounds=xbounds)
        
        # print(fitter)
        # print(_helper(xguess))
        
        return fitter
    
    def fit_ecm2(self, conc_dict, xguess, xbounds):
        '''
        Attempts to backfit K & 1/n given days to breakthrough and peak C/C0.
        Fits against a dictionary of interpolating functions.

        Parameters
        ----------
        brk_days : TYPE
            DESCRIPTION.
        peak_conc : TYPE
            DESCRIPTION.

        Returns
        -------
        int
            DESCRIPTION.

        '''
        
        def _helper(array):
            ks = array[:self.n]
            xni = array[self.n:]
            self.compdata.loc['K'] = ks * 1.
            self.compdata.loc['1/n'] = xni * 1.
            self.compdata.loc['xn'] = 1./xni
            
            self.run_ecm()
            self.plot_ecm()
            
            error = 0
            xdata = np.arange(0, np.max(conc_dict[list(conc_dict.keys())[0]].x))
            for comp in conc_dict.keys():
                ecmdata = self.plot_obj[comp](xdata)
                rawdata = conc_dict[comp](xdata)
                error += np.sum((ecmdata-rawdata)**2)
            
   
            return error

        
        fitter = minimize(_helper, xguess, method='L-BFGS-B', bounds=xbounds)

        return fitter    
    
    def plot_ecm(self):
        xdata = [0] + [i for i in np.flip(self.days_to_breakthrough) for j in range(2)] +\
                [self.days_to_breakthrough[0]+1]
        
        cmat = np.flip(np.transpose(self.cstore*self.compdata.loc['MW'].values))
        
        store_obj = {}
        cmps = np.flip(self.compdata.columns)
        for i in range(self.n):
            ydata = [0,0] + [j for j in cmat[i] for k in range(2)]
        
            fnct = interp1d(xdata, ydata, fill_value='extrapolate')
            store_obj[cmps[i]] = fnct
    
        self.plot_obj = store_obj #saves plottable dictionary back to self.
    
def _process_input_file(filename, comp_sn, bed_sn, values='molar'):
    comp_data = pd.read_excel(filename, sheet_name=comp_sn, 
                              index_col=[0], header=[0])
    
    # missing cok definition. not sure what c0_io is? or where it's defined
    comp_data.loc['xn'] = 1./comp_data.loc['1/n']
    comp_data.loc['Co'] = comp_data.loc['C0']/comp_data.loc['MW']
    comp_data.loc['Kmass'] = comp_data.loc['K']*1.
    comp_data.loc['Korig'] = comp_data.loc['K'] * 1.
    comp_data.loc['1/n_orig'] = comp_data.loc['1/n'] * 1.
    
    if values=='mass':
        comp_data.loc['K'] = comp_data.loc['K']/comp_data.loc['MW'] / \
                              (1. / comp_data.loc['MW'])**comp_data.loc['1/n']
    
    bed_data = pd.read_excel(filename, sheet_name=bed_sn, 
                             index_col=[0], header=[0])
    
    # may need to add unit conversion
    # for now, assuming the units are correct 
    bed_data = bed_data['values']
    #velocity of flow, cm/s
    bed_data['vf'] = bed_data['flrt'] * ml_pergal/(cm_perft**2)/sec_permin
    bed_data['bed_density'] *= 1e3 # double check
    
    return comp_data, bed_data



############# Some tests, commented out
# fn = 'ecm_data.xlsx'
# ecm_obj = ecm(fn)
# # print(ecm_obj.compdata)
# # print(ecm_obj.beddata)
# # print(ecm_obj.n)
# ecm_obj.run_ecm()

# print('\n\nProblem 2\n\n')

# fn = 'ecm_data2.xlsx'
# ecm_obj = ecm(fn, units='mass')
# ecm_obj.compdata.loc['C0'] *= 25.
# ecm_obj.compdata.loc['Co'] *= 25.
# ecm_obj.run_ecm()

