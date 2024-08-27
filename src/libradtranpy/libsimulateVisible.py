"""
Module libratranpy, an interface package to libradtran executable
library to simulate air transparency with LibRadTran

:author: sylvielsstfr
:creation date: November 2nd 2016
:last update: February 20th 2024
"""

import os
import re
import numpy as np
from libradtranpy import UVspec3

libradtranvers = "2.0.4"
FLAG_DEBUG = False

# Definitions and configuration
#-------------------------------------

# LibRadTran installation directory

var = 'HOME'
if var not in os.environ:
    raise EnvironmentError(f"Failed because {var} is not set.")
home = os.environ['HOME']+ '/'

if os.getenv("LIBRADTRAN_DIR"):
    libradtranpath = os.getenv('LIBRADTRANDIR')+ '/'
    libradtrandatapath = libradtranpath + "/share/libRadtran/data"
elif os.getenv("CONDA_PREFIX") != "" and os.path.isdir(os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/data")):
    libradtranpath = os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/")
    libradtrandatapath = os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/data")
else:
    raise EnvironmentError(f"\n\tYou should set a LIBRADTRAN_DIR environment variable or install rubin-libradtran package.")

#print("libradtranpath=",libradtranpath)

# Filename : RT_LS_pp_us_sa_rt_z15_wv030_oz30.txt
#          : Prog_Obs_Rte_Atm_proc_Mod_zXX_wv_XX_oz_XX
  
Prog='RT'  #definition the simulation programm is libRadTran
Obs='LS'   # definition of observatory site (LS,CT,OH,MK,...)
Rte='pp'   # pp for parallel plane of ps for pseudo-spherical
Atm=['us']   # short name of atmospheric sky here US standard and  Subarctic winter
Proc='sa'  # light interaction processes : sc for pure scattering,ab for pure absorption
           # sa for scattering and absorption, ae with aerosols default, as with aerosol special
Mod='rtvis'   # Models for absorption bands : rt for REPTRAN, lt for LOWTRAN, k2 for Kato2
ZXX='z'        # XX index for airmass z :   XX=int(10*z)
WVXX='wv'      # XX index for PWV       :   XX=int(pwv*10)
OZXX='oz'      # XX index for OZ        :   XX=int(oz/10)
AEXX='aer'
AEXX2='aer2'
CLD="cld"

# Rubin-LSST altitude from https://en.wikipedia.org/wiki/Vera_C._Rubin_Observatory
# altitude has been reduced from 2.75 km to  2.663 km     

# preselected sites 
Dict_Of_sitesAltitudes = {'LSST':2.663,
                          'CTIO':2.207,
                          'OHP':0.65,
                          'PDM':2.8905,
                          'OMK':4.205,
                          'OSL':0.000,
                           }
# pressure calculated by libradtran
Dict_Of_sitesPressures = {'LSST':731.50433,
                          'CTIO':774.6052,
                          'OHP':937.22595,
                          'PDM':710.90637,
                          'OMK':600.17224,
                          'OSL':1013.000,
                        }

Dict_Of_sitesTags = {'LSST':'LS',
                     'CTIO':'CT',
                     'OHP':'OH',
                     'PDM':'PM',
                     'OMK':'MK',
                     'OSL':'SL',
                    }

# july 2023 libradtran version
TOPTOPDIR=f"simulations/RT/{libradtranvers}/"

def CleanSimDir():
    """Remove simulation directory"""   
    os.system("rm -rf simulations")


############################################################################
def ensure_dir(f):
    """function to create a path"""
    d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)
#########################################################################



def ApplyAerosols(wl,tr,thelambda0,tau0,alpha0):
    """
     ApplyAerosols(wl,tr,thelambda0,tau0,alpha0)
     Function to provide the aerosol transmission from an analytical formula
     
     Parameters:

     :param wl: np array of wavelengths
     :type wl: float in nm unit
     :param tr: transparency array without aerosols, by example the one calculated by libRadtran
     :type tr: float
     :param thelambda0: the reference point where to have tau0 in nm
     :type thelambda0: float in nm unit
     :param tau0: is the extinction at thelambda0
     :type tau0: float
     :param alpha0: the Angstrom exponent
     :type alpha0: float
     :returns: array of aerosol transmission
     :rtype: float 
    """
    #extinc_aer=tau0*(thelambda0/wl)**alpha0
    extinc_aer=tau0*np.power(wl/thelambda0,-alpha0)
    tr_aer=np.exp(-extinc_aer)
    tr_tot=tr*tr_aer
    return tr_tot
    

def ProcessSimulation(airmass_num, pwv_num, oz_num, press_num, aer_num, angstrom_exponent_num=1.4,
                      prof_str='us', proc_str='as', cloudext=0.0, altitude="LSST", aer_lambda0=500., FLAG_VERBOSE=False):
    """
    ProcessSimulation(airmass_num,pwv_num,oz_num) 
    Function to simulate air transparency.
    No aerosol simulation is performed.
    
    Parameters:

    :param airmass_num: airmass
    :type airmass_num: float, unitless

    :param pwv_num: precipitable water vapor 
    :type pwv_num: float, in units mm

    :param oz_num: ozone colon depth 
    :type oz_num: float, in Dobson unit

    :param press_num: ground pressure 
    :type press_num: float in unit of in hPa or millibar

    :param aer_num: vertical aerosol depth
    :type aer_num: float

    :param angstrom_exponent_num: Angstrom exponent. If None or negative, default aerosol profile from libradtran is used (scaled by aer_num).
    :type angstrom_exponent_num: float

    :param prof_str: defines the type of atmosphere, such standard us,
    mid latitude summer, mid latitude winter, tropical,.., default standard us 
    :type prof_str: optional string among us,ms,mw,tp,ss,sw

    :param proc_str: activation of different processes light-air interaction, 
    like scattering and absorption (sa), absorption only (ab), scattering only (sc),..,
    default scattering and absorption and aersols (as)
    :type proc_str: optional string among sa,ab,sc,ae,as

    :param cloudext: cloud optical depth, default 0
    :type cloudext: float, optional

    :param altitude: observation site predefined (either the site name abbreviation (LSST/CTIO,OMP,OMK,MSL) or the altitude in km like 2.663 as a float
    :type altitude: string,float,int

    :param FLAG_VERBOSE: flag to activate libradtran verbose mode that print all the table generated
    :type FLAG_VERBOSE: bool

    :returns: wl, atm, arrays containing the simulated wavelengths and atmospheric transmission
    :rtype: two arrays
    """

    if FLAG_VERBOSE:
        print('--------------- ProcessSimulation -----------------------------')
        print(' 1) airmass = ', airmass_num)
        print(' 2) pwv = ', pwv_num)
        print(' 3) oz = ', oz_num)
        print(' 4) aer_num  = ', aer_num)
        print(' 5) angstrom_exp  = ', angstrom_exponent_num)
        print(' 6) wl0  = ', aer_lambda0)
        print(' 7) pressure  = ', press_num)
        print(' 8) atmospheric profile = ', prof_str)
        print(' 9) interaction processes = ', proc_str)
        print(' 10) cloud extinction = ', cloudext)
        print(' 11) site or altitude = ', altitude)
        print('--------------------------------------------')

    # altitude workaround
    if type(altitude) is str:
        if altitude in Dict_Of_sitesAltitudes.keys():
            altitude_num = Dict_Of_sitesAltitudes[altitude]
        else:
            raise Exception(f"Bad site string {altitude}")
    elif type(altitude) is float or type(altitude) is int:
        altitude_num = altitude
    else:
        raise TypeError("altitude argument must be a string (observatory name) or a numerical value in km.")
    
    # set the interaction process
    # Proc='sa'  # Pure absorption and Rayleigh scattering : Clear sky without aerosols
    if proc_str in ["sa","ab","sc","ae","as"]:
        Proc=proc_str
    else:
        raise ValueError(f'Unknown atmospheric profile {prof_str=}. Must be in ["sa","ab","sc","ae","as"].')

    # set the selected atmosphere
    if prof_str in ["us","ms","mw","tp","ss","sw"]:
        skyindex=prof_str
    else:
        raise ValueError(f'Unknown atmospheric profile {prof_str=}. Must be in ["us","ms","mw","tp","ss","sw"].')

    # Set up type of
    if Proc == 'sc':
        runtype='no_absorption'
    elif Proc == 'ab':
        runtype='no_scattering'
    elif Proc == 'sa':
        runtype='clearsky'
    elif Proc == 'ae':
        runtype='aerosol_default'
    elif Proc == 'as':
        runtype='aerosol_special'
    else:
        runtype='clearsky'

#   Selection of RTE equation solver        
    if Rte == 'pp':  # parallel plan
        rte_eq='twostr'  # 'disort' is slower for equivalent results
    elif Rte=='ps':   # pseudo spherical
        rte_eq='sdisort'
    else:
        raise ValueError(f'Unknown RTE equation solver {Rte=}.')
 
#   Selection of absorption model 
    molmodel='reptran'
    if Mod == 'rt':
        molmodel='reptran'
    if Mod == 'lt':
        molmodel='lowtran'
    if Mod == 'kt':
        molmodel='kato'
    if Mod == 'k2':
        molmodel='kato2'
    if Mod == 'fu':
        molmodel='fu'    
    if Mod == 'cr':
        molmodel='crs'     

    if re.search('us',skyindex):
        atmosphere = 'afglus'
    elif re.search('sw',skyindex):
        atmosphere = 'afglsw'
    elif re.search('ss',skyindex):
        atmosphere = 'afglss'
    elif re.search('mw',skyindex):
        atmosphere = 'afglmw'
    elif re.search('ms',skyindex):
        atmosphere = 'afglms'
    elif re.search('tp',skyindex):
        atmosphere = 'afglt'
    else:
        raise ValueError(f'Unknown atmospheric profile {skyindex=}.')

    # molecular molecular resolution: can be ['coarse','medium','fine']
    # select only COARSE Model
    molresol = 'coarse'

    # water vapor
    pwv_str='H2O '+str(pwv_num)+ ' MM'

    # Ozone
    oz_str='O3 '+str(oz_num)+ ' DU'

    # premare to create libradtran input file
    uvspec = UVspec3.UVspec()
    uvspec.inp["data_files_path"] = libradtrandatapath

    uvspec.inp["atmosphere_file"] = libradtrandatapath+'/atmmod/'+atmosphere+'.dat'
    # arbitrary earth albedo
    uvspec.inp["albedo"]           = '0.2'

    uvspec.inp["rte_solver"] = rte_eq

    if Mod == 'rtvis':
        uvspec.inp["mol_abs_param"] = molmodel + ' ' + molresol
    else:
        uvspec.inp["mol_abs_param"] = molmodel

    # Convert airmass into zenith angle (could be improved)
    sza=np.arccos(1. / airmass_num) * 180. / np.pi

    # Should be no_absorption
    if runtype=='aerosol_default':
        uvspec.inp["aerosol_default"] = ''
    elif runtype == 'aerosol_special':
        uvspec.inp["aerosol_default"] = ''
        if angstrom_exponent_num is None or angstrom_exponent_num < 0:
            uvspec.inp["aerosol_set_tau_at_wvl"] = f'500 {aer_num:.20f}'
        else:
            # below formula recover default aerosols models with angstrom_exponent_num=1.2
            tau = aer_num * (0.5 ** angstrom_exponent_num)
            uvspec.inp["aerosol_angstrom"] = f"{angstrom_exponent_num:.10f} {tau:.10f}"

    if runtype=='no_scattering':
        uvspec.inp["no_scattering"] = ''
    if runtype=='no_absorption':
        uvspec.inp["no_absorption"] = ''

    # set up the PWV and ozone value
    # Notice we have mol_modify key twice in libradtran input file, both for PWV and OZ 
    # But same key  in uvspec.inp dictionary would overwrite pwv    
    uvspec.inp["mol_modify"] = pwv_str
    uvspec.inp["mol_modify2"] = oz_str #### BUGUGUGUGUGUGUGUG ==> No because overwrite pwv. UVspec3 manage the same tag for PWV and OZ

    # rescale pressure if reasonable pressure values are provided
    if press_num > 200. and press_num < 1080.:
        uvspec.inp["pressure"] = press_num

    uvspec.inp["ic_file"] = "1D ./IC.DAT"
    uvspec.inp["ic_properties"] = "yang"
    uvspec.inp["ic_modify"] = "tau set "+str(cloudext)

    uvspec.inp["output_user"] = 'lambda edir'
    uvspec.inp["altitude"] = str(altitude_num)   # altitude observatory
    uvspec.inp["source"] = 'solar '+libradtrandatapath+'/solar_flux/kurudz_1.0nm.dat'
    uvspec.inp["sza"] = str(sza)
    uvspec.inp["phi0"] = '0'
    uvspec.inp["wavelength"] = '250.0 1200.0'
    uvspec.inp["output_quantity"] = 'reflectivity'  # 'transmittance'
    if FLAG_VERBOSE:
        uvspec.inp["verbose"] = ''
    else:
        uvspec.inp["quiet"] = ''

    # provide a usefull file
    fname = "IC.DAT"
    if not os.path.isfile(fname):
        with open(fname, 'w+') as f:
            f.write('#      z     LWC    R_eff\n')
            f.write('#     (km)  (g/m^3) (um) \n')
            f.write('     11.000   0      0   \n')
            f.write('     10.000   0.005  20  \n')

    # run libradtran
    wl, atm = uvspec.run(FLAG_VERBOSE, path=libradtranpath)
    # return wavelengths and transmissions
    return wl, atm

######################################################################################
# The main program for starting atmospheric simulation start here
#####################################################################################

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="Enter verbose (print more stuff).", default=False)
    parser.add_argument("-f", "--file", dest="filename", default="tatm.csv",
                        help="Write results in given output file name (default: tatm.csv).")
    parser.add_argument("-z", "--airmass", dest="airmass", default=1.,
                        help="Airmass (default: 1).")
    parser.add_argument("-w", "--pwv", dest="pwv", default=5.,
                        help="PWV in mm (default: 5).")
    parser.add_argument("-o", "--ozone", dest="ozone", default=300.,
                        help="Ozone in db (default: 300).")
    parser.add_argument("-c", "--cloud", dest="cloud", default=0.,
                        help="Cloud vertical optical depth (default: 0).")
    parser.add_argument("-p", "--pressure", dest="pressure", default=1000,
                        help="Pressure in hPa (default: 1000).")
    parser.add_argument("-s", "--altitudesite", dest="altitudesite", default=0.,
                        help="Site altitude in km or can be recognized observatory name as a string (default: 0).")
    parser.add_argument("-a", "--vaod", dest="vaod", default=0.,
                        help="Vertical aerosol optical depth (default: 0).")
    parser.add_argument("-e", "--exp", dest="exp", default=1.4,
                        help="Angstrom exponent for aerosols (default: 1.4).")
    parser.add_argument("-l", "--wl0", dest="wl0", default=500,
                        help="Reference wavelength for aerosol power law function in nm (default: 500).")
    parser.add_argument("-q", "--interactproc", dest="interactproc", default="sa",
                        help="Interaction processus (default: 'sa').")
    parser.add_argument("-m", "--atmmodel", dest="atmmodel", default="us",
                        help="Atmospheric model (default: 'us').")
    args = parser.parse_args()

    if args.verbose:
        FLAG_VERBOSE = True
    else:
        FLAG_VERBOSE = False

    airmass_nb = float(args.airmass)
    pwv_nb = float(args.pwv)
    oz_nb = float(args.ozone)
    press_nb = float(args.pressure)
    vaod_nb = float(args.vaod)
    exp_nb = float(args.exp)
    cld_nb = float(args.cloud)
    wl0_nb = float(args.wl0)

    # Check consistency of values
    if airmass_nb < 1 or airmass_nb > 3:
        raise ValueError("bad airmass value : z=", airmass_nb)

    if pwv_nb < 0 or pwv_nb > 50:
        raise ValueError("bad PWV value : pwv=", pwv_nb)

    if oz_nb < 0 or oz_nb > 600:
        raise ValueError("bad Ozone value : oz=", oz_nb)

    if press_nb < 0 or press_nb > 1500:
        raise ValueError("bad Pressure value : press=", press_nb)

    if cld_nb < 0 or cld_nb > 100:
        raise ValueError("bad cloud optical depth value : cld=", cld_nb)

    wl, atm = ProcessSimulation(airmass_num=airmass_nb, pwv_num=pwv_nb, oz_num=oz_nb, press_num=press_nb,
                                aer_num=vaod_nb, angstrom_exponent_num=exp_nb, aer_lambda0=wl0_nb,
                                prof_str=args.atmmodel, proc_str=args.interactproc,
                                cloudext=cld_nb, altitude=args.altitudesite, FLAG_VERBOSE=FLAG_VERBOSE)

    np.savetxt(args.filename, np.array([wl, atm]).T, delimiter=",")







   
