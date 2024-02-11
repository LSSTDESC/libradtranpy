"""
Module libratranpy, an interface package to libradtran executable
library to simulate air transparency with LibRadTran

:author: sylvielsstfr
:creation date: November 2nd 2016
:last update: October 22th 2023
"""

import logging
import os
import re
import math
import numpy as np
import sys,getopt
from libradtranpy import UVspec3

libradtranvers = "2.0.5"
FLAG_DEBUG = False

# Definitions and configuration
#-------------------------------------

# LibRadTran installation directory

var = 'HOME'
if var not in os.environ:
    raise EnvironmentError(f"Failed because {var} is not set.")
home = os.environ['HOME']+ '/'

if os.getenv("CONDA_PREFIX") != "" and os.path.isdir(os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/data")):
    libradtranpath = os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/")
    libradtrandatapath = os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/data")
elif os.getenv["LIBRADTRAN_DIR"]:
    libradtranpath = os.getenv('LIBRADTRANDIR')+ '/'
    libradtrandatapath = libradtranpath + "/share/libRadtran/data"
else:
    logging.warning(f"\n\tYou should set a LIBRADTRAN_DIR environment variable or install rubin-libradtran package.")

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



def usage():
    """Description of the usageof the script"""
    print("*******************************************************************")
    print(sys.argv[0],' [-v] -z <airmass> -w <pwv> -o <oz> -p <P> -c <cld> -m<mod> -q<proc> -s<site-string> -v<verbose-flag>')
    print(' \t - z   : airmass from 1.0 to 3.0, typical z=1 ')
    print(' \t - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm')
    print(' \t - oz  : ozone in Dobson units, typically from 100 DU to 600 DU')
    print(' \t - p   : Pressure in hPa, typical P=775.3 hPa, optional, may provide 0 to choose the standard pressure at the selected obs site')
    print(' \t - c   : Cloud vertical optical depth, optional, typical c=0')
    print(' \t - m   : Atmospheric model, optional ,typical m=\'us\' ')
    print(' \t - q   : Interaction processes, optional, typical q=\'sa\' for scattering and absorption')
    print(' \t - s   : provide site or altitude as a string : LSST/OHP/PDM/OMK/OSL or altitude in km like akm_2.663, optional, default LSST')
    print(' \t - v   : activate libradtran output verbose mode, optional, default no verbose mode')
    
    print('\t Examples : ')
    print('\t \t 1) python libsimulateVisible.py -z 1 -w 0 -o 0 -s LSST')
    print('\t \t 2) python libsimulateVisible.py -z 1 -w 4 -o 300 -c 0 -p 742 -m us -q sa -s LSST')
    
    print('\t To generate ascii printout of the used atmospheric model tables in a log file :')
    print('\t \t python libsimulateVisible.py -v -z 1 -w 0 -o 0 -s LSST >& output.log')
    
    print('\t Actually provided : ')
    print('\t \t Number of arguments:', len(sys.argv), 'arguments.')
    print('\t \t Argument List:', str(sys.argv))
    
    print("*******************************************************************")
    
def usageaer():
    """Description of the usage of the script"""
    print("*******************************************************************")
    print(sys.argv[0],' [-v] -z <airmass> -w <pwv> -o <oz> -a<aer> -p <P> -c <cld> -m<mod> -q<proc> -s<altitude> -v<verbose-flag>')
    print(' \t - z   : airmass from 1.0 to 3.0, typical z=1 ')
    print(' \t - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm')
    print(' \t - oz  : ozone in Dobson units, typically  from 100 DU to 600 DU')
    print(' \t - aer : Aerosols vertical optical depth, typical a=0.04')
    print(' \t - p   : Pressure in hPa, typical P=775.3 hPa, optional, can be set to 0 to get the standard pressure at that site  ')
    print(' \t - c   : Cloud vertical optical depth,optional ,typical c=0')
    print(' \t - m   : Atmospheric model, optional, typical m=\'us\' ')
    print(' \t - q   : Interaction processes, optional ,typical q=\'sa\' for scattering and absorption')
    print(' \t - s   : provide site or altitude : LSST/OHP/PDM/OMK/PDM or altitude in km like akm_2.663, optional, default is LSST')
    print(' \t - v   : activate libradtran output verbose to get atmospheric profile, optional, default no verbose mode')
   
    print('\t Examples : ')
    print('\t \t 1) python libsimulateVisible.py -z 1 -w 0 -o 0 -a 0 -s LSST')
    print('\t \t 2) python libsimulateVisible.py -z 1 -w 4 -o 300 -a 0.3 -c 0 -p 742 -m us -q sa -s LSST')
    
    print('\t To generate ascii printout of the used atmospheric model table in a log file :')
    print('\t python libsimulateVisible.py -v -z 1 -w 0 -o 0 -s LSST >& output.log')
    #python libsimulateVisible.py -v -z 1 -w 0 -o 0 -s LSST >& output.log

    print('\t Actually provided : ')
    print('\t \t Number of arguments:', len(sys.argv), 'arguments.')
    print('\t \t Argument List:', str(sys.argv))
    print("*******************************************************************")

#----------------------------------------------------------------------------
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
    
    
#-----------------------------------------------------------------------------

def ProcessSimulation(airmass_num,pwv_num,oz_num,press_num,vaod,angstrom_exponent,prof_str='us',proc_str='as',cloudext=0.0, altitude_str ="LSST",FLAG_VERBOSE=False):
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

    :param vaod: vertical aerosol depth
    :type vaod: float

    :param angstrom_exponent: Angstrom exponent
    :type angstrom_exponent: float

    :param prof_str: defines the type of atmosphere, such standard us,
    mid latitude summer, mid latitude winter, tropical,.., default standard us 
    :type prof_str: optional string among us,ms,mw,tp,ss,sw

    :param proc_str: activation of different processes light-air interaction, 
    like scattering and absorption (sa), absorption only (ab), scattering only (sc),..,
    default scattering and absorption and aersols (as)
    :type proc_str: optional string among sa,ab,sc,ae,as

    :param cloudext: cloud optical depth, default 0
    :type cloudext: float, optional

    :param altitude_str: observation site predefined (either the site name abrebiation (LSST/CTIO,OMP,OMK,MSL) or the string on altitude like akm_2.663)
    :type altitude_str: string

    :param FLAG_VERBOSE: flag to activate libradtran verbose mode that print all the table generated
    :type FLAG_VERBOSE: bool

    :returns: OUTPUTDIR,outputFilename, path and filename of datafile containing the simulated data     
    :rtype: two strings
    """


    if FLAG_DEBUG:
        print('--------------- ProcessSimulation -----------------------------')
        print(' 1) airmass = ', airmass_num)
        print(' 2) pwv = ', pwv_num)
        print(' 3) oz = ', oz_num)
        print(' 4) pressure  = ',press_num)
        print(' 5) atmospheric profile = ',prof_str)
        print(' 6) interaction processes = ',proc_str)
        print(' 7) cloud extinction = ',cloudext)
        print(' 8) site or altitude = ', altitude_str)
        print('--------------------------------------------')

   
    # altitude workaround
    #if Dict_Of_sitesAltitudes.get(altitude_str):
    if altitude_str in Dict_Of_sitesAltitudes.keys():    
        altitude_num = Dict_Of_sitesAltitudes[altitude_str]
        altitude_dir = altitude_str
        Obs = Dict_Of_sitesTags[altitude_str] # for the path of input/output
    elif altitude_str[:4] == "akm_":
        height_str = altitude_str[4:]
        altitude_num = float(height_str)
        altitude_dir = height_str.replace(".","_")
    else:
        raise Exception(f"Bad altitude/site string {altitude_str}")

    # keep the altutude to write in input file
    OBS_Altitude = altitude_num

    if FLAG_DEBUG:
        print(f"Observation site altitude for libradran sim : {OBS_Altitude} km")
    
    # set the interaction process
    
    #Proc='sa'  # Pure absorption and Rayleigh scattering : Clear sky without aerosols
    if proc_str in ["sa","ab","sc","ae","as"]:
        Proc=proc_str
    else:
        raise ValueError(f'Unknown atmospheric profile {prof_str=}. Must be in ["sa","ab","sc","ae","as"].')

    # set the selected atmosphere
    if prof_str in ["us","ms","mw","tp","ss","sw"]:
        skyindex=prof_str
    else:
        raise ValueError(f'Unknown atmospheric profile {prof_str=}. Must be in ["us","ms","mw","tp","ss","sw"].')

    # create output dir
    TOPDIR = os.path.join(TOPTOPDIR,altitude_dir)
    ensure_dir(TOPDIR)

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
    if Rte == 'pp': # parallel plan
        rte_eq='twostr'  # 'disort'
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

      
    # loop on molecular model resolution
    #molecularresolution = np.array(['COARSE','MEDIUM','FINE'])
    # select only COARSE Model
    molecularresolution = np.array(['COARSE'])
    for molres in molecularresolution:
        if molres=='COARSE':
            molresol ='coarse'
        elif molres=='MEDIUM':
            molresol ='medium'
        else:
            molresol ='fine'


    #water vapor
    pwv_val=pwv_num
    pwv_str='H2O '+str(pwv_val)+ ' MM'

    # airmass
    airmass=airmass_num

    # Ozone
    oz_str='O3 '+str(oz_num)+ ' DU'

    verbose = FLAG_DEBUG

    # premare to create libradtran input file

    uvspec = UVspec3.UVspec()
    uvspec.inp["data_files_path"]  =  libradtrandatapath

    uvspec.inp["atmosphere_file"] = libradtrandatapath+'/atmmod/'+atmosphere+'.dat'
    # arbitrary earth albedo
    uvspec.inp["albedo"]           = '0.2'

    uvspec.inp["rte_solver"] = rte_eq

    if Mod == 'rtvis':
        uvspec.inp["mol_abs_param"] = molmodel + ' ' + molresol
    else:
        uvspec.inp["mol_abs_param"] = molmodel

    # Convert airmass into zenith angle (could be inproved)
    sza=np.arccos(1. / airmass) * 180. / np.pi

    # Should be no_absorption
    if runtype=='aerosol_default':
        uvspec.inp["aerosol_default"] = ''
    elif runtype == 'aerosol_special':
        uvspec.inp["aerosol_default"] = ''
        if angstrom_exponent is None or angstrom_exponent < 0:
            uvspec.inp["aerosol_set_tau_at_wvl"] = f'500 {vaod:.20f}'
        else:
            # below formula recover default aerosols models with angstrom_exponent=1.2
            tau = vaod * (0.5 ** angstrom_exponent)
            uvspec.inp["aerosol_angstrom"] = f"{angstrom_exponent:.10f} {tau:.10f}"

    if runtype=='no_scattering':
        uvspec.inp["no_scattering"] = ''
    if runtype=='no_absorption':
        uvspec.inp["no_absorption"] = ''

    # set up the ozone value
    uvspec.inp["mol_modify"] = pwv_str
    uvspec.inp["mol_modify"] = oz_str

    # rescale pressure   if reasonable pressure values are provided
    if press_num>200. and press_num<1080.:
        uvspec.inp["pressure"] = press_num

    uvspec.inp["ic_file"] = "1D ./IC.DAT"
    uvspec.inp["ic_properties"] = "yang"
    uvspec.inp["ic_modify"] = "tau set "+str(cloudext)

    uvspec.inp["output_user"] = 'lambda edir'
    uvspec.inp["altitude"] = OBS_Altitude   # Altitude  observatory
    uvspec.inp["source"] = 'solar '+libradtrandatapath+'/solar_flux/kurudz_1.0nm.dat'
    uvspec.inp["sza"]        = str(sza)
    uvspec.inp["phi0"]       = '0'
    uvspec.inp["wavelength"]       = '250.0 1200.0'
    uvspec.inp["output_quantity"] = 'reflectivity' #'transmittance' #
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
    wl, atm = uvspec.run(verbose,path=libradtranpath)
    # return wavelengths and transmissions
    return wl, atm

#---------------------------------------------------------------------------


#------------------------------------------------------------------------------
def ProcessSimulationaer(airmass_num,pwv_num,oz_num,aer_num,press_num,prof_str='us',proc_str='sa',cloudext=0.0, altitude_str='LSST',FLAG_VERBOSE=False):
    """
    ProcessSimulationaer(airmass_num,pwv_num,oz_num,aer_num,press_num) 
    with aerosol simulation is performed
    
    Parameters:

    :param airmass_num: airmass
    :type airmass_num: float, unitless

    :param pwv_num: precipitable water vapor 
    :type pwv_num: float, in units mm

    :param oz_num: ozone colon depth 
    :type oz_num: float, in Dobson unit

    :param aer_num: vertical aerosol optical depth
    :type aer_num: fliat, unitless

    :param press_num: ground pressure 
    :type press_num: float in unit of in hPa or millibar

    :param prof_str: defines the type of atmosphere, such standard us, 
    mid latitude summer, mid latitude winter, tropical,.., default standard us 
    :type prof_str: optional string among us,ms,mw,tp,ss,sw

    :param proc_str: activation of different processes light-air interaction, 
    like scattering and absorption (sa), absorption only (ab), scattering only (sc),..,
    default scattering and absorption (sa)
    :type proc_str: optional string among sa,ab,sc,ae,as

    :param cloudext: cloud optical depth, default 0
    :type cloudext: float, optional

    :param altitude_str: observation site predefined (either the site name abrebiation (LSST/CTIO,OMP,OMK,MSL) or the string on altitude like akm_2.663)
    :type altitude_str: string

    :param FLAG_VERBOSE: flag to activate libradtran verbose mode that print all the table generated
    :type FLAG_VERBOSE: bool
    
    :returns: OUTPUTDIR,outputFilename, path and filename of datafile containing the simulated data     
    :rtype: two strings
    """

    
 
    if FLAG_DEBUG:
        print('------------- ProcessSimulationaer -------------------------------')
        print(' 1) airmass = ', airmass_num)
        print(' 2) pwv = ', pwv_num)
        print(' 3) oz = ', oz_num)
        print(' 4) aer = ',aer_num)
        print(' 5) pressure =',press_num)
        print(' 6) profile =',prof_str)
        print(' 7) interaction processes = ',proc_str)
        print(' 8) cloud extinction = ',cloudext)
        print(' 9) site or altitude = ', altitude_str)
        print('--------------------------------------------')
    

    # altitude workaround
    #if Dict_Of_sitesAltitudes.get(altitude_str):
    if altitude_str in Dict_Of_sitesAltitudes.keys():  
        altitude_num = Dict_Of_sitesAltitudes[altitude_str]
        altitude_dir = altitude_str
        Obs = Dict_Of_sitesTags[altitude_str] # for the path of input/output
    elif altitude_str[:4] == "akm_":
        height_str = altitude_str[4:]
        altitude_num = float(height_str)
        altitude_dir = height_str.replace(".","_")
    else:
        raise Exception(f"Bad altitude/site string {altitude_str}")
        altitude_num = Dict_Of_sitesAltitudes['LSST']
        altitude_dir = 'LSST0'

    # keep altitude to write it in input file
    OBS_Altitude = altitude_num

    if FLAG_DEBUG:
        print(f"Observation site altitude for libradran sim : {OBS_Altitude} km")
    
    #Proc='sa'  # Pure absorption and Rayleigh scattering : Clear sky without aerosols
    if proc_str in ["sa","ab","sc","as","ae"]:
        Proc=proc_str
        
        
    # set the selected atmosphere
    if prof_str in ["us","ms","mw","tp","ss","sw"]:
        Atm=[prof_str]
    
    # create output dir
    TOPDIR = os.path.join(TOPTOPDIR,altitude_dir)
    ensure_dir(TOPDIR)

    
    # build the part 1 of filename
    BaseFilename_part1=Prog+'_'+Obs+'_'+Rte+'_'
    
    aerosol_string = '500 '+str(aer_num)
    #aerosol_str=str(wl0_num)+ ' '+str(tau0_num)
    aer_index=int(aer_num*100.)

    # Set up type of run
    runtype='aerosol_special' #'no_scattering' #aerosol_special #aerosol_default# #'clearsky'#
    Proc='as'
    
    #Proc='as'  # Absoprtion + Rayleigh + aerosols special
    
    
    if Proc == 'sc':
        runtype='no_absorption'
        outtext='no_absorption'
    elif Proc == 'ab':
        runtype='no_scattering'
        outtext='no_scattering'
    elif Proc == 'sa':
        runtype='clearsky'
        outtext='clearsky'
    elif Proc == 'ae':   
        runtype='aerosol_default'
        outtext='aerosol_default'
    elif Proc == 'as':   
        runtype='aerosol_special'
        outtext='aerosol_special'
    else:
        runtype='clearsky'
        outtext='clearsky'

#   Selection of RTE equation solver        
    if Rte == 'pp': # parallel plan
        rte_eq='disort'
    elif Rte=='ps':   # pseudo spherical
        rte_eq='sdisort'
        
 
#   Selection of absorption model 
    molmodel='reptran'
    if Mod == 'rtvis':
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
               
    	  
    # for simulation select only two atmosphere   
    #theatmospheres = np.array(['afglus','afglms','afglmw','afglt','afglss','afglsw'])
    atmosphere_map=dict()  # map atmospheric names to short names 
    atmosphere_map['afglus']='us'
    atmosphere_map['afglms']='ms'
    atmosphere_map['afglmw']='mw'  
    atmosphere_map['afglt']='tp'  
    atmosphere_map['afglss']='ss'  
    atmosphere_map['afglsw']='sw'  
      
    theatmospheres= []

    for skyindex in Atm:
        if re.search('us',skyindex):
            theatmospheres.append('afglus')
        if re.search('sw',skyindex):
            theatmospheres.append('afglsw')
        if re.search('ss',skyindex):
            theatmospheres.append('afglss')
        if re.search('mw',skyindex):
            theatmospheres.append('afglmw')
        if re.search('ms',skyindex):
            theatmospheres.append('afglms')
        if re.search('tp',skyindex):
            theatmospheres.append('afglt')

    # 1) LOOP ON ATMOSPHERE
    for atmosphere in theatmospheres:
        #if atmosphere != 'afglus':  # just take us standard sky
        #    break
        atmkey=atmosphere_map[atmosphere]
       
        # manage input and output directories and vary the ozone
        TOPDIR2=TOPDIR+'/'+Rte+'/'+atmkey+'/'+Proc+'/'+Mod
        ensure_dir(TOPDIR2)
        INPUTDIR=TOPDIR2+'/'+'in'
        ensure_dir(INPUTDIR)
        OUTPUTDIR=TOPDIR2+'/'+'out'
        ensure_dir(OUTPUTDIR)
    
    
        # loop on molecular model resolution
        #molecularresolution = np.array(['COARSE','MEDIUM','FINE']) 
        # select only COARSE Model
        molecularresolution = np.array(['COARSE'])    
        for molres in molecularresolution:
            if molres=='COARSE':
                molresol ='coarse'
            elif molres=='MEDIUM':
                molresol ='medium'
            else:
                molresol ='fine'
           
        
        #water vapor   
        pwv_val=pwv_num
        pwv_str='H2O '+str(pwv_val)+ ' MM'
        wvfileindex=int(10*pwv_val)
        
       
        # airmass
        airmass=airmass_num
        amfileindex=int(airmass_num*10)
        
        # Ozone    
        oz_str='O3 '+str(oz_num)+ ' DU'
        ozfileindex=int(oz_num/10.)

        # Cloud
        cldindex = str(int(cloudext * 1000))     
        cld_str=cldindex.zfill(4)
        
        # root path of input/output files    
        BaseFilename = BaseFilename_part1+atmkey+'_'+Proc+'_'+Mod+'_z'+str(amfileindex)+'_'+WVXX+str(wvfileindex) +'_'+OZXX+str(ozfileindex)+'_'+AEXX+str(aer_index)+"_"+CLD+cld_str
                    
        verbose = FLAG_DEBUG

        # prepare libradtran input file

        uvspec = UVspec3.UVspec()
        
        #uvspec.inp["data_files_path"]  =  libradtranpath+'data'
        uvspec.inp["data_files_path"]  =  libradtrandatapath
                
        #uvspec.inp["atmosphere_file"] = libradtranpath+'data/atmmod/'+atmosphere+'.dat'
        uvspec.inp["atmosphere_file"] = libradtrandatapath+'/atmmod/'+atmosphere+'.dat'
        
        # choose arbitrary earth albedo
        uvspec.inp["albedo"]           = '0.2'
    
        uvspec.inp["rte_solver"] = rte_eq
              
                
        if Mod == 'rtvis':
            uvspec.inp["mol_abs_param"] = molmodel + ' ' + molresol
        else:
            uvspec.inp["mol_abs_param"] = molmodel

        # Convert airmass into zenith angle 
        am=airmass
        sza=math.acos(1./am)*180./math.pi

        # Should be no_absorption
        if runtype=='aerosol_default':
            uvspec.inp["aerosol_default"] = ''
        elif runtype=='aerosol_special':
            uvspec.inp["aerosol_default"] = ''
            uvspec.inp["aerosol_set_tau_at_wvl"] = aerosol_string
                   
        if runtype=='no_scattering':
            uvspec.inp["no_scattering"] = ''
        if runtype=='no_absorption':
            uvspec.inp["no_absorption"] = ''
     
        # set up the ozone value               
        uvspec.inp["mol_modify"] = pwv_str
        uvspec.inp["mol_modify"] = oz_str
        
        # rescale pressure   if reasonable pressure values are provided
        if press_num>200. and press_num<1080.:
            uvspec.inp["pressure"] = press_num
        else:
            if FLAG_VERBOSE:
                print("creazy pressure p=",press_num, ' hPa')

        uvspec.inp["ic_file"] = "1D ./IC.DAT"
        uvspec.inp["ic_properties"] = "yang"
        uvspec.inp["ic_modify"] = "tau set " + str(cloudext)

        uvspec.inp["output_user"] = 'lambda edir'
        uvspec.inp["altitude"] = OBS_Altitude   # Altitude LSST observatory
        uvspec.inp["source"] = 'solar '+libradtrandatapath+'/solar_flux/kurudz_1.0nm.dat'
        #uvspec.inp["source"] = 'solar '+libradtranpath+'data/solar_flux/kurudz_1.0nm.dat'
        #uvspec.inp["source"] = 'solar '+libradtranpath+'data/solar_flux/kurudz_0.1nm.dat'
        uvspec.inp["sza"]        = str(sza)
        uvspec.inp["phi0"]       = '0'
        uvspec.inp["wavelength"]       = '250.0 1200.0'
        uvspec.inp["output_quantity"] = 'reflectivity' #'transmittance' #

        if FLAG_VERBOSE:
            uvspec.inp["verbose"] = ''
        else:
            uvspec.inp["quiet"] = ''

  

        if "output_quantity" in uvspec.inp.keys():
            outtextfinal=outtext+'_'+uvspec.inp["output_quantity"]

           
            
        inputFilename=BaseFilename+'.INP'
        outputFilename=BaseFilename+'.OUT'
        inp=os.path.join(INPUTDIR,inputFilename)
        out=os.path.join(OUTPUTDIR,outputFilename)
        
         # provide a usefull file
        fname = "IC.DAT"
        if not os.path.isfile(fname):
            with open(fname, 'w+') as f:
                f.write('#      z     LWC    R_eff\n')
                f.write('#     (km)  (g/m^3) (um) \n')
                f.write('     11.000   0      0   \n')
                f.write('     10.000   0.005  20  \n')    

        # write libradtran input file    
        uvspec.write_input(inp)
        # execute libradtran
        uvspec.run_with_files(inp,out,verbose,path=libradtranpath)
        
    # return path for libratran output file    
    return OUTPUTDIR,outputFilename

#---------------------------------------------------------------------------




######################################################################################
# The main program for starting atmospheric simulation start here
#####################################################################################

if __name__ == "__main__":
    
    # Activate or not the mode with simulation of aerosols inside libradtran
    AerosolTest_Flag = False
    
    # init string variables
    # airmass
    airmass_str=""
    # precipitable water vapor
    pwv_str=""
    # ozone
    oz_str=""
    # ground pressure
    press_str=""
    # aerosol optical depth
    aer_str=""
    # aeorosol reference wavelength
    wl0_str=""
    # cloud optical depth
    tau0_str=""
    # cloud optical depth
    cld_str=""
    # atmospheric model to be used
    # 'us':'afglus','ms':'afglms','mw':'afglmw','tp':'afglt','ss':'afglss','sw':'afglsw'
    model_str=""
    # interaction process=""
    # 'sa': scattering and absorption, 'sc' : scattering only , 'ab' : absorption only
    proc_str=""

    # altitude or site string
    alt_str=""
    
    # Case No Aerosols
    # airmass_num,pwv_num,oz_num,press_num,prof_str='us',proc_str='sa',cloudext=0.0
    
    # No aerosol, just call function ProcessSimulation()
    if AerosolTest_Flag == False:
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hvz:w:o:p:c:m:q:s:",["z=","w=","o=","p=","c=","m=","q=","s="])
        except getopt.GetoptError:
            print(' Exception bad getopt with :: '+sys.argv[0]+ ' [-v] -z <airmass> -w <pwv> -o <oz> -p <press> -c <cldvod> -m <atmmodel> -q <interactproc> -s<altitudesite-string>')
            sys.exit(2)
        
    
        FLAG_VERBOSE=False
        print('opts = ',opts)
        print('args = ',args)
        
        
        for opt, arg in opts:
            if opt == '-h':
                usage()
                sys.exit()
            elif opt == '-v':
                FLAG_VERBOSE=True   
            elif opt in ("-z", "--airmass"):
                airmass_str = arg
            elif opt in ("-w", "--pwv"):
                pwv_str = arg
            elif opt in ("-o", "--oz"):
                oz_str = arg  
            elif opt in ("-p", "--pr"):
                press_str = arg
            elif opt in ("-c", "--cld"):
                cld_str = arg
            elif opt in ("-m", "--atm"):
                model_str = arg
            elif opt in ("-q", "--qp"):
                proc_str = arg
            elif opt in ("-s", "--site"):
                alt_str = arg

            else:
                print('Do not understand arguments : ',sys.argv)
            
         
        print('---Decoded the following args : ------')
        print('1) airmass-str = ', airmass_str)
        print('2) pwv-str = ', pwv_str)
        print("3) oz-str = ", oz_str)
        print("4) pr = ", press_str)
        print("5) cld = ", cld_str)
        print("6) mod = ",model_str)
        print("7) proc = ",proc_str)
        print("8) alt/site = ",alt_str)
        print("9) FLAG_VERBOSE = ",FLAG_VERBOSE)
        print("10) FLAG_DEBUG = ",FLAG_DEBUG)
        print('--------------------------------------------')

        # mandatory arguments
        if airmass_str=="":
            usage()
            sys.exit()

        if pwv_str=="":
            usage()
            sys.exit()

        if oz_str=="":
            usage()
            sys.exit()

        if alt_str=="":
            usage()
            sys.exit()
             
        
        airmass_nb=float(airmass_str)
        pwv_nb=float(pwv_str)
        oz_nb=float(oz_str)	
    
        
        # optional arguments
        # pressure
        if press_str=="":
            #this force to use pressure of the altitude
            press_str="0.0"
        press_nb=float(press_str)
          
        # cloud
        if cld_str=="":
            cld_str = "0"    
        cld_nb = float(cld_str)
        
        if model_str=="":
            model_str="us"
        if proc_str=="":
            proc_str="sa"
        
        print('----Selected parameters : -----------------')
        print('1) airmass  = ', airmass_nb)
        print('2) pwv = ', pwv_nb)
        print("3) oz = ", oz_nb)
        print("4) press = ", press_nb)
        print("5) cld = ", cld_nb)
        print("6) atm model = ", model_str)
        print("7) interaction model = ", proc_str)
        print('--------------------------------------------')
        
    
        # Check consistency of values
        
        if airmass_nb<1 or airmass_nb >3 :
            print("bad airmass value : z=",airmass_nb)
            sys.exit()
            
        if pwv_nb<0 or pwv_nb >50 :
            print("bad PWV value : pwv=",pwv_nb)
            sys.exit()
        
        if oz_nb<0 or oz_nb >600 :
            print("bad Ozone value : oz=",oz_nb)
            sys.exit()
              
        if press_nb<0 or press_nb >1500 :
            print("bad Pressure value : press=",press_nb)
            sys.exit()
        
        if cld_nb<0 or cld_nb >100 :
            print("bad cloud optical depth value : cld=",cld_nb)
            sys.exit()
        
        
        # do the simulation now 
        print("all arguments values are OK, start libradtran simulation")
        
        #ProcessSimulation(airmass_num,pwv_num,oz_num,press_num,prof_str='us',proc_str='sa',cloudext=0.0, FLAG_VERBOSE=False):
        path, outputfile=ProcessSimulation(airmass_nb,pwv_nb,oz_nb,press_nb,model_str,proc_str=proc_str,cloudext=cld_nb ,altitude_str=alt_str,FLAG_VERBOSE=FLAG_VERBOSE)
    
        print('*****************************************************')
        print(' path       = ', path)
        print(' outputfile =  ', outputfile)
        print('*****************************************************')
    
    # With aerosol, just call function ProcessSimulationaer()
    else:
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hvz:w:o:a:p:c:m:q:s:",["z=","w=","o=","a=","p=","c=","m=","q=","s="])
        except getopt.GetoptError:
            print(' Exception bad getopt with :: '+sys.argv[0]+ ' -z <airmass> -w <pwv> -o <oz> -a <aer> -p <press> -c <cldvod>-m <model> -q <interaction> -s <altitude-string>')
            sys.exit(2)
        
    
        
        print('opts = ',opts)
        print('args = ',args)
        
        FLAG_VERBOSE = False
        
        for opt, arg in opts:
            if opt == '-h':
                usageaer()
                sys.exit()
            elif opt == "-v":
                FLAG_VERBOSE = True
            elif opt in ("-z", "--airmass"):
                airmass_str = arg
            elif opt in ("-w", "--pwv"):
                pwv_str = arg
            elif opt in ("-o", "--oz"):
                oz_str = arg  
            elif opt in ("-a", "--aer"):
                aer_str = arg 
            elif opt in ("-p", "--pr"):
                press_str = arg
            elif opt in ("-c", "--cld"):
                cld_str = arg
            elif opt in ("-m", "--am"):
                model_str = arg
            elif opt in ("-q", "--qp"):
                proc_str = arg
            elif opt in ("-s", "--site"):
                alt_str = arg
            else:
                print('Do not understand arguments : ',sys.argv)
            
         
        print('---Decode the following arguments ----------')
        print('1) airmass-str = ', airmass_str)
        print('2) pwv-str = ', pwv_str)
        print("3) oz-str = ", oz_str)
        print("4) aer = ", aer_str)
        print("5) press = ", press_str)
        print("6) cld = ", cld_str)
        print("7) mod = ", model_str)
        print("8) proc = ", proc_str)
        print("9) alt/site = ",alt_str)
        print("10) FLAG_VERBOSE = ",FLAG_VERBOSE)
        print("11) FLAG_DEBUG = ",FLAG_DEBUG)
        print('--------------------------------------------')

        # mandatory arguments
        if airmass_str=="":
            usageaer()
            sys.exit()

        if pwv_str=="":
            usageaer()
            sys.exit()

        if oz_str=="":
            usageaer()
            sys.exit()

        if alt_str=="":
            usage()
            sys.exit()
            
        airmass_nb=float(airmass_str)
        pwv_nb=float(pwv_str)
        oz_nb=float(oz_str)	
        
        
        # optional arguments
        # pressure 
        if press_str== "":
            # set it to zero to force libradtran to use the pressure for the altitude
            press_str= "0.0"
        press_nb=float(press_str)  
        
        #aerosols
        if aer_str== "":
            aer_str = "0.0"    
        aer_nb = float(aer_str)
         
        # cloud
        if cld_str=="":
            cld_str = "0.0"    
        cld_nb = float(cld_str)
        
        if model_str == "":
            model_str = "us"
        if proc_str == "":
            proc_str = "sa"
	
        
        print('--------------------------------------------')
        print('1) airmass  = ', airmass_nb)
        print('2) pwv = ', pwv_nb)
        print("3) oz = ", oz_nb)
        print("4) aer = ", aer_nb)
        print("5) press = ", press_nb)
        print("6) cld = ", cld_nb)
        print("7) atm model = ", model_str)
        print("8) interaction process = ", proc_str)
        print("9) FLAG_VERBOSE = ",FLAG_VERBOSE)
        print("10) FLAG_DEBUG = ",FLAG_DEBUG)
        print('--------------------------------------------')
        
        
        # Check the consistency of values
    
        if airmass_nb<1 or airmass_nb >3 :
            print("bad airmass value z=",airmass_nb)
            sys.exit()
            
        if pwv_nb<0 or pwv_nb >50 :
            print("bad PWV value pwv=",pwv_nb)
            sys.exit()
        
        if oz_nb<0 or oz_nb >600 :
            print("bad Ozone value oz=",oz_nb)
            sys.exit()
            
        if aer_nb<0 or aer_nb >0.5 :
            print("bad Aerosol value aer=",aer_nb)
            sys.exit()
        
        
        if press_nb<0 or press_nb >1500 :
            print("bad Pressure value press=",press_nb)
            sys.exit()
            
        if cld_nb<0 or cld_nb >100 :
            print("bad cloud optical depth value : cld=",cld_nb)
            sys.exit()
        
        
        # do the simulation now 
        print("values are OK")
        #ProcessSimulationaer(airmass_num,pwv_num,oz_num,aer_num,press_num,prof_str='us',proc_str='sa',cloudext=0.0, FLAG_VERBOSE=False):
        path, outputfile=ProcessSimulationaer(airmass_nb,pwv_nb,oz_nb,aer_nb,press_nb,model_str,proc_str,cloudext=cld_nb,altitude_str=alt_str,FLAG_VERBOSE=FLAG_VERBOSE)
    
        print('*****************************************************')
        print(' path       = ', path)
        print(' outputfile =  ', outputfile)
        print('*****************************************************')





   
