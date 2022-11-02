#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  equivalentWidth.py
#  
#  Copyright 2022  Joseph Chevalier joseph.chevalier@ijclab.in2p3.fr
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


"""
Author : Joseph Chevalier
Affiliation : IJClab-IN2P3
Creation date : October 20, 2022
"""

"""
Latest update : November 02, 2022
"""

"""
Purpose : to calculate the equivalent width of an absorption (or emission) band in a spectrum.
As used for the O2 band on Rubin AuxTel, for atmospheric calibration
"""
###########
# IMPORTS #
###########

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import (MinMaxInterval,SqrtStretch,ZScaleInterval,PercentileInterval,ImageNormalize)
from astropy.visualization.stretch import SinhStretch, LinearStretch,AsinhStretch,LogStretch
from itertools import repeat
from matplotlib.colors import LogNorm
import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.stats import norm
import pickle

from libradtranpy import libsimulateVisible

#####################################
## Main function definition, could ##
##    be done outside this block   ##
#####################################
def fit_gaussian(x, flux, abs_min, abs_max, sigma_flux=None, central_lambda=None):

    def fun_fit(t, a, b, mu, sigma, k):
        o2abs = norm(loc=mu, scale=sigma)
        return (a*t + b)-k*o2abs.pdf(t)
    
    mask_for_fit = (x >= abs_min) * (x <= abs_max)
    xmasked=x[mask_for_fit]
    if central_lambda is None:
        central_lambda = xmasked[ np.where( flux[mask_for_fit] == np.min(flux[mask_for_fit]) )[0] ][0] ## caution : the minimum could be in the continuum if the interval around the band is too wide
    
    #print("dbg - O2 line center : {} nm".format(central_lambda))
    if sigma_flux is not None:
		p1, cov = curve_fit(fun_fit, xmasked, flux[mask_for_fit], p0=[0., 1.0, central_lambda, 2.0, 1.0], sigma=sigma_flux[mask_for_fit])
	else:
		p1, cov = curve_fit(fun_fit, xmasked, flux[mask_for_fit], p0=[0., 1.0, central_lambda, 2.0, 1.0], sigma=None)
    min_cont, max_cont = p1[2] - 5*p1[3], p1[2] + 5*p1[3]
    #print("DEBUG: borne inf = {}nm, borne sup. = {}nm".format(min_cont, max_cont)) ## Add instruction for debug in fine mode
    xmask_continuum = (x >= min_cont) * (x <= max_cont)
    
    while len(x[xmask_continuum]) < len(p1):
	    min_cont -= 1.
	    max_cont += 1.
	    xmask_continuum = (x >= min_cont) * (x <= max_cont)
	if sigma_flux is not None:
		mod, cov = curve_fit(fun_fit, x[xmask_continuum], flux[xmask_continuum], p0=[p1[0], p1[1], p1[2], p1[3], p1[4]], sigma=sigma_flux[xmask_continuum])
	else:
		mod, cov = curve_fit(fun_fit, x[xmask_continuum], flux[xmask_continuum], p0=[p1[0], p1[1], p1[2], p1[3], p1[4]], sigma=None)
		
    min_lin, max_lin = mod[2] - 3*mod[3], mod[2] + 3*mod[3]
    min_cont, max_cont = mod[2] - 5*mod[3], mod[2] + 5*mod[3]
    limits_ = np.array([min_cont, min_lin, max_lin, max_cont])
    
    return mod, cov, limits_


def eqw_norm(x, flux, abs_min, abs_max, sigma_flux=None, central_lambda=None, gaussMod_band=None, limits=None, fit_band=True, return_fit=False, make_plot=True, plot_name='test.png'): 
    
    def continuum_error(l0, sigma_a, sigma_b):
        return np.sqrt( np.power(sigma_a*l0, 2.) + np.power(sigma_b, 2.) )

    def area_ul_error(flux_error, lambda_vec):
        delta_lambda = np.mean( np.diff(lambda_vec) )
        return np.sqrt( np.sum( np.power(flux_error, 2.) ) * np.power(delta_lambda, 2.) )

    def area_l_error(sigma_area_cont, sigma_area_ul):
        return np.sqrt( np.power(sigma_area_cont, 2.) + np.power(sigma_area_ul, 2.) )

    def eqw_error(cont_at_line, line_area, sigma_cont, sigma_area_line):
        sigma2_eqw = np.power(sigma_area_line / cont_at_line, 2.) + np.power(line_area*sigma_cont / np.power(cont_at_line, 2.), 2.)
        return np.sqrt(sigma2_eqw)
    
    if gaussMod_band is None:
        fit_band=True
    
    if fit_band:
        gaussMod_band, cov, lims = fit_gaussian(x, flux, abs_min, abs_max, sigma_flux, central_lambda)
        init=[gaussMod_band[0], gaussMod_band[1]]
        o2_band_centre = gaussMod_band[2]
        sigmaBand = gaussMod_band[3]
        kBand = gaussMod_band[4]
        min_cont, min_lin, max_lin, max_cont = lims[0], lims[1], lims[2], lims[3]
        #delta_line = np.sqrt(cov[2,2])
    else:
        o2_band_centre = gaussMod_band[2]
        sigmaBand = gaussMod_band[3]
        kBand = gaussMod_band[4]
        if limits is None:
            min_lin, max_lin = o2_band_centre - 3*sigmaBand, o2_band_centre + 3*sigmaBand
            min_cont, max_cont = o2_band_centre - 5*sigmaBand, o2_band_centre + 5*sigmaBand
        else:
            min_cont, min_lin, max_lin, max_cont = limits[0], limits[1], limits[2], limits[3]
        init = [ (flux[-1]-flux[0]) / (x[-1]-x[0]), flux[0] - ( (flux[-1]-flux[0]) / (x[-1]-x[0]) )*x[0] ]
    
    def lin_fun(x, a, b):
        return a*x+b
    
    limits_ = np.array([min_cont, min_lin, max_lin, max_cont])
    xmask_contForFit = (x >= min_cont)*(x<=min_lin) + (x<=max_cont)*(x>=max_lin)
    if sigma_flux is not None:
		contMod, cov = curve_fit(lin_fun, x[xmask_contForFit], flux[xmask_contForFit], p0=init, sigma=sigma_flux)
	else:
		contMod, cov = curve_fit(lin_fun, x[xmask_contForFit], flux[xmask_contForFit], p0=init, sigma=None)
    
    def continuum(x):
        return contMod[0]*x+contMod[1]

    xmask_continuum = (x>=min_cont)*(x<=max_cont)

    mask_for_fit = (x >= abs_min) * (x <= abs_max)
    xline0 = x[mask_for_fit]
    fline0 = flux[mask_for_fit]
    if sigma_flux is not None:
		fline_err0 = sigma_flux[mask_for_fit]
    
    xline = x[xmask_continuum]
    fline = flux[xmask_continuum]
    if sigma_flux is not None:
		fline_err = sigma_flux[xmask_continuum]
    xmask_line = (xline>=min_lin)*(xline<=max_lin)
    xline_model = xline[xmask_line]
    
    def spec_mod(x):
        o2abs = norm(loc=o2_band_centre, scale=sigmaBand)
        return continuum(x)-kBand*o2abs.pdf(x)
    
    fline_model = spec_mod(xline_model)
    
    def normed_cont(x): ## must return 1.
        return continuum(x) / continuum(x)
    def normed_spec(x):
        return spec_mod(x) / continuum(x)
    
    norm_array = continuum(xline)
    
    #area_c, area_c_err = quad(normed_cont, min_lin, max_lin)
    area_c, area_c_err = quad(normed_cont, min_cont, max_cont)
    #area_ul = quad(normed_spec, min_lin, max_lin)[0]
    #area_ul = np.trapz(fline[xmask_line]/norm_array[xmask_line], xline_model)
    area_ul = np.trapz(fline/norm_array, xline)
    #area_c = quad(normed_cont, min_cont, max_cont)[0]
    #area_ul = quad(normed_spec, min_cont, max_cont)[0]
    area_l = area_c-area_ul
    #print(area_c)
    #print(area_ul)
    #print(area_l)
    
    cont_min = continuum(o2_band_centre)
    eqw = area_l #/cont_min
    #print('EQW = ',eqw)
        
    sigma_a = np.sqrt(cov[0,0])
    sigma_b = np.sqrt(cov[1,1])
    
    #print(cont_min)
    #print(area_l)
    #print(area_c_err)
    
    cont_err_ = continuum_error(o2_band_centre, sigma_a, sigma_b)
    #print(cont_err_)
    if sigma_flux is not None:
		area_ul_err_ = area_ul_error(fline_err[xmask_line]/continuum(norm_array[xmask_line]), xline_model)
	else:
		area_ul_err_ = area_ul_error(np.ones_like(xline_model)*0.02, xline_model) ## default 2% error to ensure that the code runs
		
    #print(area_ul_err_)
    
    area_l_err_ = area_l_error(area_c_err, area_ul_err_)
    #print(area_l_err_)
    
    #eqw_err_ = eqw_error(cont_min, area_l_, cont_err_, area_l_err_)
    eqw_err_ = eqw_error(1., area_l, cont_err_, area_l_err_)
    #print(eqw_err_)
    
    #print('----------')
    
    eqwmask = (xline_model>=o2_band_centre-eqw/2.)*(xline_model<=o2_band_centre+eqw/2.)
    #print(continuum)
    #print(xmin_line-eqw/2.,xmin_line+eqw/2.)
    
    if make_plot:
        fig, axs=plt.subplots(1,2)#figsize=(10,6), constrained_layout=True)
        axs=axs.ravel()
        axs[0].scatter(xline0, fline0, marker="+", label="Observed fluxes")
        axs[0].plot(xline0, continuum(xline0), ls='--', color='orange', label="Continuum model")
        axs[0].plot(xline0, spec_mod(xline0), ls='-', color='purple', label="Flux model")
        axs[0].axvline(x=min_lin,color='b',ls='--', label="3 sigma : O2 absorption model boundaries")
        axs[0].axvline(x=max_lin,color='b',ls='--')
        axs[0].axvline(x=min_cont,color='r',ls='--', label="5 sigma : boundaries of the model")
        axs[0].axvline(x=max_cont,color='r',ls='--')
        axs[0].fill_between(xline, continuum(xline), fline, color='gray',alpha=0.2, label="Ref. surface")
        #axs[0].fill_between(xline_model[eqwmask], continuum(xline_model[eqwmask]), color='cyan', alpha=0.2, label="Rectangle of EQW")
        #axs[0].fill_between(xline_model, continuum(xline_model), spec_mod(xline_model), color='gray',alpha=0.2, label="Ref. surface")
        axs[0].fill_between(xline_model[eqwmask], continuum(xline_model[eqwmask]), color='cyan', alpha=0.2, label="Rectangle of EQW")
        axs[1].set_xlabel(r'$\lambda \, [nm]$',fontsize=12)
        axs[0].set_ylabel(r'$\gamma \, [s \cdot cm^2 \cdot nm]$',fontsize=12)
        axs[0].grid()
        axs[0].legend(loc="lower right")
        
        
        axs[1].scatter(xline, fline/norm_array, marker="+", label="Normalised observed fluxes")
        axs[1].plot(xline, normed_cont(xline), ls='--', color='orange', label="Continuum model")
        axs[1].plot(xline, normed_spec(xline), ls='-', color='purple', label="Flux model")
        axs[1].axvline(x=min_lin,color='b',ls='--', label="3 sigma : O2 absorption model boundaries")
        axs[1].axvline(x=max_lin,color='b',ls='--')
        axs[1].axvline(x=min_cont,color='r',ls='--', label="5 sigma : boundaries of the model")
        axs[1].axvline(x=max_cont,color='r',ls='--')
        axs[1].axvline(x=o2_band_centre,color='k',ls=':', label="O2 line center")
        axs[1].fill_between(xline, normed_cont(xline), fline/norm_array, color='gray',alpha=0.2, label="Ref. surface")
        #axs[1].fill_between(xline_model[eqwmask], normed_cont(xline_model[eqwmask]), color='cyan', alpha=0.2, label="Rectangle of EQW")
        #axs[1].fill_between(xline_model, normed_cont(xline_model), normed_spec(xline_model), color='gray',alpha=0.2, label="Ref. surface")
        axs[1].fill_between(xline_model[eqwmask], normed_cont(xline_model[eqwmask]), color='cyan', alpha=0.2, label="Rectangle of EQW")
        axs[1].set_xlabel(r'$\lambda \, [nm]$',fontsize=12)
        axs[1].set_ylabel(r'$\gamma \, / continuum [-]$',fontsize=12)
        axs[1].grid()
        axs[1].legend(loc="lower right")
        
        fig.suptitle("Figure : "+plot_name)
        '''
        plot_outdir = 'output_plots/eqw_line_fit_{0}/'.format(source_spec)
        if os.path.exists(plot_outdir)==False:
            os.mkdir(plot_outdir)
        plt.savefig(plot_outdir+plot_name) 
        '''
    if return_fit:
        return eqw, sigmaBand, eqw_err_, spec_mod, continuum, limits_
    else:
        return eqw, sigmaBand, eqw_err_, area_c_err, cont_err_, area_ul_err_, area_l_err_
        
def simAtmEqw(airmasses, pressures, wl_min=700., wl_max=800., wl_mid=761.0, lims=None, pwv=0.0, oz=300.0, atm='us', inter='sa', clouds=0.0, flagVerbose=False):
	'''
	airmasses and pressures must be 1D-lists with the same length.
	Returns the equivalent widths and associated errors for the specified atmosphere(s) (possibility to use "all")
	For each couple of values : airmass, pressure.
	The result is returned as a list of arrays : - dimension along axis 0 equals len(returned list)=len(airmasses)=len(pressures)
	                                             - dimension along axis 1 equals the number of atmospheres (1 or 6 if "all" is specified) -- array axis 0
	                                             - dimension along axis 2 equals 3 : atmosphere name, eq. width, eq. width error -- array axis 1
	                                             
    DEFAULT VALUES ARE FOR THE O2 ABSORPTION BAND OF THE ATMOSPHERE (test case of this program)
	'''
	results_am_press = []
	for am, press in zip(airmasses, pressures):
		outpaths, outfiles = libsimulateVisible.ProcessSimulation(am, pwv, oz, press, prof_str=atm, proc_str=inter, cloudext=clouds, FLAG_VERBOSE=flagVerbose)
		filepaths = [ os.path.join(path_, file_) for path_, file_ in zip(outpaths, outfiles) ]
		eqws = []
		atms = []
		eqwErrs = []
		for dataFile in filepaths:
			data = np.loadtxt(dataFile)
			atmType = dataFile.split('/')[-5]
			print(atmType)
			wl = data[:,0]
			transm = data[:,1]
			gaussMod, gaussCov, gaussLimits = fit_gaussian(wl, transm, wl_min, wl_max, central_lambda=wl_mid)
			eqw_, sigmaBand_, eqwErr_, areaC_err , cont_err_, areaUL_err_, areaL_err_= eqw_norm(wl, transm, wl_min, wl_max, central_lambda=gaussMod[2], gaussMod_band=gaussMod, limits=lims, fit_band=False, return_fit=False, make_plot=False)
			eqws.append(eqw_)
			eqwErrs.append(eqwErr_)
			atms.append(atmType)
		eqwArr = np.column_stack((atms, eqws, eqwErrs))
		results_am_press.append(eqwArr)
	return results_am_press

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
