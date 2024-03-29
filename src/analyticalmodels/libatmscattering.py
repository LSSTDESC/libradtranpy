#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 20:19:44 2018
Author: Sylvie Dagoret-Campagne 
Affliiation : LAL/IN2P3/CNRS
"""

import numpy as np
import pandas as pd


import matplotlib
import matplotlib.pyplot as plt



N_A = 6.0221409e23  # mol-1
R = 8.3144598       # J/(K.mol)
g0 = 9.80665        #  m/s2

M_air = 28.965338*1e-3  # u.g/u.mol  (kg/mol)
M_air_dry = 28.9644*1e-3 # *u.g/u.mol (kg/mol)
M_h2o = 18.016*1e-3      # *u.g/u.mol  (kg/mol)

P0 = 101325.         # *u.Pa;   /*!< Pa : pressure at see level */
T0 = 288.15          # *u.K;   #/*!< sea level temperature */  
L = 0.0065           #*u.K/u.m  # refroidissement en fonction de l'altitude

LSST_Altitude = 2750  # in meters from astropy package (Cerro Pachon)
CTIO_Altitude = 2200  # in meters from astropy package (Cerro Pachon)
OHP_Altitude = 650   # in km
PDM_Altitude= 2890.5 # in km Pic du Midi

altitude0=LSST_Altitude  # m altitude at LSST




#-----------------------------------------------------------------
def Temperature_adiabatic(h):
    """Temperature vs altitude
    
    :param h: altitude 
    :type h: float in meters

    :return: temperature
    :rtype: float in unit degree Kelvin 
    """
    return  T0 - L * h
#------------------------------------------------------------------------------------
def Pressure_isothermal(altitude):
    """
    Pressure at the altitude for an isothermal atmosphere.
    For dry air.
    
    :param altitude: input altitude
    :type altitude: float in meters
    :return: pressure  
    :rtype: float in unit Pa SI 
    """
    h = altitude
    P = P0*np.exp(-((g0*M_air_dry)/(R*T0))*h)
    return P  
#------------------------------------------------------------------------------------
def Pressure_adiabatic(altitude):
    """
    Pressure at the altitude for an adiabatic atmosphere.
    
    :param altitude: input altitude
    :type altitude: float in unit meters
    :return: pressure  
    :rtype: float in unit Pa SI 
    """
    h = altitude
    P=P0*np.exp(g0*M_air_dry/R/L*np.log(1-L*h/T0))
    return P
#---------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
def MassDensity_isothermal(altitude):
    """
    Mass density Density for Dry Air in kg per m^3
    Provide the density at the altitude.
    For dry air.

    :param altitude: input altitude in meters
    :type altitude: float in unit meters.
    :return: density 
    :rtype: float in unit kg/m^3
    """
    h = altitude
    T = T0
    rho = Pressure_isothermal(h) *  M_air_dry / (R * T)
    return rho


# ------------------------------------------------------------------------------------
def MassDensity_adiabatic(altitude,T=0):
    """
    Mass density Density for Dry Air in kg per m^3
    Provide the density at the altitude.
    For dry air.

    :param altitude: altitude
    :type altitude: flotat in unit meter
    :param T: Temperature. If T=0, use average altitude profile
    :type T: float in degree K.
    :return: density  
    :rtype: float in unit kg/m^3
    """
    h = altitude

    if T==0:
        T = T0 - L * h
    rho = Pressure_adiabatic(h) * M_air_dry / (R * T)
    return rho
# ------------------------------------------------------------------------------------
def MassDensity_adiabatic_humid(altitude,Hr,T=0):
    """
    Mass density Density for Dry Air in kg per m^3
    Provide the density at the corresponding altitude.
    
    :param altitude: altitude 
    :type altitude: float in unit meter
    :param Hr: relative humidity
    :type Hr: float
    :param T: Temperature. If T=0, use average altitude profile
    :type T: float in unit degree Kelvin.
    :return: density 
    :rtype: float in unit kg/m^3
    """

    h = altitude

    if T==0:
        T = T0 - L * h

    T_c=T-273.15

    p = Pressure_adiabatic(h)
    rho=1/(287.06*T)*(p-230.617*Hr*np.exp(17.5043*T_c/(241.2+T_c)))
    return rho

# ------------------------------------------------------------------------------------
def AtomeDensity_isothermal(altitude):
    """
    atome Density for Dry Air ( double altitude)

    Provide the atome density at the altitude.
    Attention, ici on considère de l'air sec.

    :param altitude: input altitude 
    :type: float in unit meters
    :return: density 
    :rtype: float in unit atom/ m^3
    """
    h = altitude
    T = T0
    rho = Pressure_isothermal(h) * N_A / (R * T)

    return rho


# ------------------------------------------------------------------------------------
def AtomeDensity_adiabatic(altitude,T=0):
    """
    atome Density for Dry Air ( double altitude)

    Provide the atome density at the altitude.
    Attention, ici on considère de l'air sec.

    :param altitude: altitude 
    :type altitude: float in unit meter
    :param T:  Temperature. If T=0, use average altitude profile
    :type T: float in unit degree Kelvin
    :return: density  
    :rtype:  float in atom / m^3
    """
    h = altitude

    if T==0:
        T = T0 - L * h
    rho = Pressure_adiabatic(h) * N_A / (R * T)
    return rho

# ---------------------------------------------------------------------------------



def XDepth_isothermal(altitude,costh=1):
    """
    Function : XDepth(altitude,costh)
    Provide the column depth in gr / cm^2 equivalent of airmass in physical units
    for an isothermal atmosphere.

    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return:  XDepth  column depth 
    :rtype: float in unit gr per cm squared
    """
    h=altitude
    XD=Pressure_isothermal(h)/g0/costh
    return XD
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
def XDepth_adiabatic(altitude,costh=1):
    """
    Function : XDepth(altitude,costh)
    Provide the column depth in gr / cm^2 equivalent of airmass in physical units
    for an adiabatic atmosphere.

    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return:  XDepth  column depth 
    :rtype: float in unit gr per cm squared
    """
    h=altitude
    XD=Pressure_adiabatic(h)/g0/costh
    return XD
#---------------------------------------------------------------------------------
def RayOptDepth_adiabatic(wavelength, altitude=altitude0, costh=1):
    """
    Function RayOptDepth(double wavelength, double altitude, double costh)
    Provide Rayleigh optical depth in an adiabatic atmosphere.

    :param wavelength: wavelength
    :type wavelength: float in unit nm
    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return: the optical depth no unit, for Rayleigh
    :rtype: float
    """

    h=altitude

    #A=(XDepth(h,costh)/(3102.*u.g/(u.cm*u.cm)))
    A = (XDepth_adiabatic(h,costh)/(3102.*1e-3/(1e-4)))
    B = np.exp(-4.*np.log(wavelength/(400.)))  
    C = 1-0.0722*np.exp(-2*np.log(wavelength/(400)))

    OD = A*B/C
        
    #double OD=XDepth(altitude,costh)/2970.*np.power((wavelength/400.),-4);

    return OD
#-----------------------------------------------------------------------------------
def RayOptDepth_isothermal(wavelength, altitude=altitude0, costh=1):
    """
    Function RayOptDepth(double wavelength, double altitude, double costh)
    Provide Rayleigh optical depth for an isothermal atmosphere

    :param wavelength: wavelength
    :type wavelength: float in unit nm
    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return: the optical depth no unit, for Rayleigh scattering
    :rtype: float
         
    """

    h = altitude

    #A=(XDepth_adiab(h,costh)/(3102.*u.g/(u.cm*u.cm)))
    #B=np.exp(-4.*np.log(wavelength/(400.*u.nm)))  
    #C= 1-0.0722*np.exp(-2*np.log(wavelength/(400.*u.nm)))
    
    A=(XDepth_isothermal(h,costh)/(31020.))
    B=np.exp(-4.*np.log(wavelength/(400.)))  
    C= 1-0.0722*np.exp(-2*np.log(wavelength/(400.)))

    OD=A*B/C
        
    #double OD=XDepth(altitude,costh)/2970.*np.power((wavelength/400.),-4);

    return OD   
#------------------------------------------------------------------------------------
def RayOptDepth2_adiabatic(wavelength, altitude=altitude0, costh=1):
    """
    Function RayOptDepth2(wavelength, altitude, costh)
    Provide Rayleigh optical depth for an adiabatic atmosphere
    
    :param wavelength: wavelength
    :type wavelength: float in unit nm
    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return: the optical depth no unit, for Rayleigh scattering
    :rtype: float

    WORSE !
    
    """
    h = altitude
    #A=XDepth(h,costh)/(2770.*u.g/(u.cm*u.cm))
    #B=np.exp(-4*np.log(wavelength/(400.*u.nm)))
    
    A = XDepth_adiabatic(h,costh)/(27700.)
    B = np.exp(-4*np.log(wavelength/(400)))
      
    OD = A*B         
  
    return OD
#-----------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
def RayOptDepth2_isothermal(wavelength, altitude=altitude0, costh=1):
    """
    Function RayOptDepth2(wavelength, altitude, costh)
    Provide Rayleigh optical depth for an isothermal atmosphere.

    :param wavelength: wavelength
    :type wavelength: float in unit nm
    :param  altitude: input altitude
    :type altitude: float in unit meter
    :param costh: cosimus of zenith angle
    :type costh: float 
    :return: the optical depth no unit, for Rayleigh scattering
    :rtype: float
    
    """
    h = altitude
    #A=XDepth(h,costh)/(2770.*u.g/(u.cm*u.cm))
    #B=np.exp(-4*np.log(wavelength/(400.*u.nm)))
    
    A = XDepth_isothermal(h,costh)/(27700.)
    B = np.exp(-4*np.log(wavelength/(400)))
      
    OD = A*B         
  
    return OD
#-----------------------------------------------------------------------------------
def AeroOptDepth(wavelength,tau_aerosols_500=0.05,alpha_ang=1) :
    """
    AeroOptDepth(wavelength, alpha)
    Provide Vertical Aerosols optical depth

    :param wavelength: wavelength
    :type wavelength: float in unit nm
    :param tau_aerosols_500: VAOD at 500 nm
    :type tau_aerosols_500: float
    :param alpha_ang: Angstrom exponent, must be positive.
    :type alpha_ang: float
    :return: aerosol optical depth
    :rtype: float
    
    """

    OD = tau_aerosols_500*np.exp(-alpha_ang*np.log(wavelength/(500)))
    return OD


#-----------------------------------------------------------------------------------


def RayOptDepthXD(wavelength, xdepth):
    """
    Function RayOptDepthXD(wavelength, xdepth)
    Provide Rayleigh optical depth

    :param wavelength: wavelength
    :type wavelength: float in nm
    :param xdepth: atmospheric depth
    :type xdepth: float in unit  g/cm2
    :return:  optical depth no unit, for Rayleigh
    :rtype: float 
    """

    #A = xdepth / (3102. * u.g / (u.cm * u.cm))
    #B = np.exp(-4. * np.log(wavelength / (400. * u.nm)))
    #C = 1 - 0.0722 * np.exp(-2 * np.log(wavelength / (400. * u.nm)))

    A = xdepth / 3102.
    B = np.exp(-4. * np.log(wavelength / 400.))
    C = 1 - 0.0722 * np.exp(-2 * np.log(wavelength / 400. ))

    OD = A * B / C

    return OD
#========================================================================================================

if __name__ == "__main__":
    
    
    #---------------------------------------
    h=np.linspace(0,20000.,20)
    p1=Pressure_isothermal(h)
    p2=Pressure_adiabatic(h)
    pmin=p2.min()
    pmax=p2.max()
    p_pdm=np.interp(altitude0, h, p2)
    print(p_pdm/100.,"hecto pascal")
    plt.figure()
    plt.plot(p1/P0,h/1000.,'b-',label='isothermal')
    plt.plot(p2/P0,h/1000.,'r-',label='adiabatic')
    plt.title('Thermodynamic model of Altitude versus Pressure')
    plt.plot([pmin/P0,pmax/P0],[altitude0/1000.,altitude0/1000.],"g:")
    plt.plot([p_pdm/P0, p_pdm/P0],[h.min()/1000., h.max()/1000.], "g:")
    plt.xlabel('Pressure (atm)')
    plt.ylabel('Altitude (km)')
    plt.legend()
    plt.grid(True)
    plt.show()
    #---------------------------------------
    
    
    XD_adiabatic=XDepth_adiabatic(h)
    XD_isothermal=XDepth_isothermal(h)

    xmin=XD_adiabatic.min()
    xmax = XD_adiabatic.max()
    x_pdm = np.interp(altitude0, h, XD_adiabatic)

    print(x_pdm / 100., "kg/m^2")

    plt.figure()
    plt.plot(XD_adiabatic,h/1000.,'r-',label='adiabatic')
    plt.plot(XD_isothermal,h/1000.,'b-',label='isothermal')

    plt.plot([xmin, xmax], [altitude0 / 1000., altitude0 / 1000.], "g:")
    plt.plot([x_pdm , x_pdm ], [h.min() / 1000., h.max() / 1000.], "g:")

    plt.title('Thermodynamic model of Altitude versus Atmospheric depth')
    plt.xlabel('X Depth $(kg/m^2)$')
    plt.ylabel('Altitude (km)')
    plt.legend()
    plt.grid(True)
    plt.show()
    #---------------------------------------
    
    
    wavelength=np.linspace(200.,1100.,100)  # in nm
    
    od_isothermal=RayOptDepth_isothermal(wavelength)
    od_adiabatic=RayOptDepth_adiabatic(wavelength)
    
    plt.figure()
    plt.plot(wavelength,np.exp(-od_isothermal),'b-',label='isothermal')
    plt.plot(wavelength,np.exp(-od_adiabatic),'r-',label='adiabatic')
    plt.title('Model 1 for Rayleigh Scattering')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Air transmittance')
    plt.grid(True)
    plt.legend()
    plt.show()
    #-------------------------------------------------------------------------
    
    wavelength=np.linspace(200.,1100.,100)  # in nm
    
    od=RayOptDepth_isothermal(wavelength)
    od1_adiab=RayOptDepth_adiabatic(wavelength)
    od2_adiab=RayOptDepth2_adiabatic(wavelength)
    
    plt.figure()
    plt.plot(wavelength,np.exp(-od_isothermal),'k:',label='formula 1 : isothermal')
    plt.plot(wavelength,np.exp(-od1_adiab),'r-',label='formula 1 : adiabatic')
    plt.plot(wavelength,np.exp(-od2_adiab),'b-',label='formula 2 : adiabatic')
    plt.title('Model for Rayleigh Scattering')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Air transmittance')
    plt.legend()
    plt.grid(True)
    #--------------------------------------------------------------------------------
    
    AOD=AeroOptDepth(wavelength)
    
    plt.figure()
    plt.plot(wavelength,np.exp(-od1_adiab),'r-',label='Rayleigh')
    plt.plot(wavelength,np.exp(-AOD),'b-',label='Aerosols')
    plt.title('Model for Rayleigh and Aerosols Scattering')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Air transmittance')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    #---------------------------------------------------------------------------------
    
    
    
    
    
