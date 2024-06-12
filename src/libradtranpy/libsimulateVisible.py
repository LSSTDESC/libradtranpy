"""
Module libratranpy, an interface package to libradtran executable
library to simulate air transparency with LibRadTran

:author: sylvielsstfr
:creation date: November 2nd 2016
:last update: February 20th 2024
"""

import os
import numpy as np
from libradtranpy import UVspec3

libradtranvers = "2.0.5"
FLAG_DEBUG = False

# Definitions and configuration
# -------------------------------------

# LibRadTran installation directory

var = "HOME"
if var not in os.environ:
    raise EnvironmentError(f"Failed because {var} is not set.")
home = os.environ["HOME"] + "/"

if os.getenv("LIBRADTRANDIR"):
    libradtranpath = os.getenv("LIBRADTRANDIR")
    libradtrandatapath = os.path.join(libradtranpath, "share/libRadtran/data")
elif os.getenv("CONDA_PREFIX") != "" and os.path.isdir(
    os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/data")
):
    libradtranpath = os.path.join(os.getenv("CONDA_PREFIX"), "share/libRadtran/")
    libradtrandatapath = os.path.join(
        os.getenv("CONDA_PREFIX"), "share/libRadtran/data"
    )
else:
    raise EnvironmentError(
        "\n\tYou should set a LIBRADTRAN_DIR environment variable or install rubin-libradtran package."
    )

# Filename : RT_LS_pp_us_sa_rt_z15_wv030_oz30.txt
#          : Prog_Obs_Rte_Atm_proc_Mod_zXX_wv_XX_oz_XX

ZXX = "z"  # XX index for airmass z :   XX=int(10*z)
WVXX = "wv"  # XX index for PWV       :   XX=int(pwv*10)
OZXX = "oz"  # XX index for OZ        :   XX=int(oz/10)
AEXX = "aer"
AEXX2 = "aer2"
CLD = "cld"

# Rubin-LSST altitude from https://en.wikipedia.org/wiki/Vera_C._Rubin_Observatory
# altitude has been reduced from 2.75 km to  2.663 km

# preselected sites
Dict_Of_sitesAltitudes = {
    "LSST": 2.663,
    "CTIO": 2.207,
    "OHP": 0.65,
    "PDM": 2.8905,
    "OMK": 4.205,
    "OSL": 0.000,
}
# pressure calculated by libradtran
Dict_Of_sitesPressures = {
    "LSST": 731.50433,
    "CTIO": 774.6052,
    "OHP": 937.22595,
    "PDM": 710.90637,
    "OMK": 600.17224,
    "OSL": 1013.000,
}

Dict_Of_sitesTags = {
    "LSST": "LS",
    "CTIO": "CT",
    "OHP": "OH",
    "PDM": "PM",
    "OMK": "MK",
    "OSL": "SL",
}

# july 2023 libradtran version
TOPTOPDIR = f"simulations/RT/{libradtranvers}/"


def CleanSimDir():
    """Remove simulation directory"""
    os.system("rm -rf simulations")


############################################################################
def ensure_dir(f):
    """function to create a path"""
    if not os.path.exists(f):
        os.makedirs(f)
        
        
def zenith_angle_from_airmass(A):
    # Check if input is a float or an array
    if isinstance(A, (int, float)):
        A = np.array([A])  # Convert float to an array for consistent processing
    elif isinstance(A, list):
        A = np.array(A)  # Convert list to a NumPy array

    # Calculate zenith angle based on conditions
    result = np.empty_like(A, dtype=float)

    # For A < 1.5
    mask = A < 1.5
    result[mask] = np.rad2deg(np.arccos(1.0 / A[mask]))

    # For A >= 1.5
    mask = A >= 1.5
    term1 = 1
    term2 = 10.0 / 6371.0
    term3 = (A[mask]**2 / 2) * (10.0 / 6371.0)
    cos_z = (term1 + term2 - term3) / A[mask]
    result[mask] = np.rad2deg(np.arccos(cos_z))

    # If input was originally a single float, return a single float
    if result.size == 1:
        return result[0]
    return result


#########################################################################


def ApplyAerosols(wl, tr, thelambda0, tau0, alpha0):
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
    # extinc_aer=tau0*(thelambda0/wl)**alpha0
    extinc_aer = tau0 * np.power(wl / thelambda0, -alpha0)
    tr_aer = np.exp(-extinc_aer)
    tr_tot = tr * tr_aer
    return tr_tot


def ProcessSimulation(
    model="reptran",
    airmass=1.0,
    pwv=0.0,
    ozone=0.0,
    pressure=1013,
    atm_model="afglus",
    rte="twostr",
    proc="clearsky",
    cloudext=0.0,
    albedo=0,
    altitude="LSST",
    lambda_min=300,
    lambda_max=1200,
    wl_res="coarse",
    temperature=273.15,
    pseudospherical=True,
    output_quantity="reflectivity",
    FLAG_VERBOSE=False
):
    """
    ProcessSimulation(airmass,pwv,ozone)
    Function to simulate air transparency.
    No aerosol simulation is performed.

    Parameters:

    :param airmass: airmass
    :type airmass: float, unitless

    :param pwv: precipitable water vapor
    :type pwv: float, in units mm

    :param ozone: ozone colon depth
    :type ozone: float, in Dobson unit

    :param pressure: ground pressure
    :type pressure: float in unit of in hPa or millibar

    :param aod: vertical aerosol depth
    :type aod: float

    :param angstrom_exponent_num: Angstrom exponent. If None or negative, default aerosol profile from libradtran is used (scaled by aod).
    :type angstrom_exponent_num: float

    :param atm_model: defines the type of atmosphere, such standard us,
    mid latitude summer, mid latitude winter, tropical,.., default standard us
    :type atm_model: optional string among us,ms,mw,tp,ss,sw

    :param proc: activation of different processes light-air interaction,
    like scattering and absorption (sa), absorption only (ab), scattering only (sc),..,
    default scattering and absorption and aersols (as)
    :type proc: optional string among sa,ab,sc,ae,as

    :param cloudext: cloud optical depth, default 0
    :type cloudext: float, optional

    :param altitude: observation site predefined (either the site name abbreviation (LSST/CTIO,OMP,OMK,MSL) or the altitude in km like 2.663 as a float
    :type altitude: string,float,int

    :param FLAG_VERBOSE: flag to activate libradtran verbose mode that print all the table generated
    :type FLAG_VERBOSE: bool

    :returns: wl, atm, arrays containing the simulated wavelengths and atmospheric transmission
    :rtype: two arrays
    """

    # if FLAG_VERBOSE:
    #     print("--------------- ProcessSimulation -----------------------------")
    #     print(" 1) airmass = ", airmass)
    #     print(" 2) pwv = ", pwv)
    #     print(" 3) oz = ", ozone)
    #     print(" 4) aod  = ", aod)
    #     print(" 5) angstrom_exp  = ", angstrom_exponent_num)
    #     print(" 6) wl0  = ", aer_lambda0)
    #     print(" 7) pressure  = ", pressure)
    #     print(" 8) atmospheric profile = ", atm_model)
    #     print(" 9) interaction processes = ", proc)
    #     print(" 10) cloud extinction = ", cloudext)
    #     print(" 11) site or altitude = ", altitude)
    #     print("--------------------------------------------")
    
    
    # premare to create libradtran input file
    uvspec = UVspec3.UVspec()
    uvspec.inp["data_files_path"] = os.path.join(
        os.path.expanduser("~"), libradtrandatapath
    )
    uvspec.inp["atmosphere_file"] = os.path.join(
        libradtrandatapath, "atmmod/" + atm_model + ".dat"
    )

    # altitude workaround
    if isinstance(altitude, str):
        if altitude in Dict_Of_sitesAltitudes.keys():
            altitude_num = Dict_Of_sitesAltitudes[altitude]
        else:
            raise Exception(f"Bad site string {altitude}")
    elif isinstance(altitude, float) or isinstance(altitude, int):
        altitude_num = altitude
    else:
        raise TypeError(
            "altitude argument must be a string (observatory name) or a numerical value in km."
        )

    # set the interaction process
    if proc in ["no_absorption", "no_scattering", "clearsky", "aerosol_default", "aerosol_special"]:
        pass
    else:
        raise ValueError(
            f'Unknown atmospheric process {atm_model=}. Must be in ["no_absorption", "no_scattering", "clearsky", "aerosol_default", "aerosol_special"].'
        )

    # set the selected atmosphere
    if atm_model in [
        "afglus",
        "afglms",
        "afglmw",
        "afgltp",
        "afglss",
        "afglsw",
        "ohp_winter",
        "ohp_spring",
        "ohp_summer",
        "ohp_autumn",
    ]:
        pass
    else:
        raise ValueError(
            f'Unknown atmospheric profile {atm_model}. Must be in ["us","ms","mw","tp","ss","sw", "ohp_winter", "ohp_spring", "ohp_summer", "ohp_autumn"].'
        )

    # # Set up type of
    # if Proc == "sc":
    #     runtype = "no_absorption"
    # elif Proc == "ab":
    #     runtype = "no_scattering"
    # elif Proc == "sa":
    #     runtype = "clearsky"
    # elif Proc == "ae":
    #     runtype = "aerosol_default"
    # elif Proc == "as":
    #     runtype = "aerosol_special"
    # else:
    #     runtype = "clearsky"

    
    # Albedo
    uvspec.inp["albedo"] = str(albedo)
    # Radiative transfer solver
    uvspec.inp["rte_solver"] = rte
    # Molecular resolution
    uvspec.inp["mol_abs_param"] = model + " " + wl_res
    # Surface temperature
    uvspec.inp["sur_temperature"] = str(temperature)
    

    # Should be no_absorption
    if proc == "aerosol_default":
        uvspec.inp["aerosol_default"] = ""
    if proc == "no_scattering":
        uvspec.inp["no_scattering"] = ""
    if proc == "no_absorption":
        uvspec.inp["no_absorption"] = ""

    # PWV value
    uvspec.inp["mol_modify H2O"] = f"{pwv:.20f} MM"
    # Ozone value
    uvspec.inp["mol_modify O3"] = f"{ozone:.20f} DU"

    # rescale pressure if reasonable pressure values are provided
    if pressure > 200.0 and pressure < 1080.0:
        uvspec.inp["pressure"] = pressure
    else:
        print(f"ERROR, INPUT PRESSURE = {pressure} IS UNFEASIBLE")

    if cloudext > 0:
        uvspec.inp["ic_file"] = "1D ./IC.DAT"
        uvspec.inp["ic_properties"] = "yang"
        uvspec.inp["ic_modify"] = "tau set " + str(cloudext)

    # Configure desired output
    if pseudospherical:
        uvspec.inp["output_user"] = "lambda edir" + "\n" + "pseudospherical"
    else:
        uvspec.inp["output_user"] = "lambda edir"
    
    # Set up the altitude of the site
    uvspec.inp["altitude"] = str(altitude_num)  # altitude observatory
    
    

    # Only for than X < 2,a fter that need to implement approximate code
    # zenith_angle = np.rad2deg(np.arccos(1.0 / airmass))
    zenith_angle = zenith_angle_from_airmass(airmass)
    
    # print(zenith_angle)
    uvspec.inp["sza"] = str(zenith_angle)
    # looking downward = umu > 0
    # looking upward = umu < 0
    # The radiance is output at phi and umu.
    # uvspec.inp["umu"] = (
    #     str(-np.cos(np.deg2rad(zenith_angle)))
    #     .replace("[", "")
    #     .replace("]", "")
    #     .replace("\n", "")
    # )
    # uvspec.inp["sza"] = (
    #     str(zenith_angle)
    #     .replace("[", "")
    #     .replace("]", "")
    #     .replace("\n", "")
    # )
    
    
    # print(uvspec.inp["umu"])
    # uvspec.inp["source"] = "solar"
    uvspec.inp["phi"] = "0"
    uvspec.inp["phi0"] = "0"
    uvspec.inp["wavelength"] = f"{lambda_min} {lambda_max}"
    uvspec.inp["output_quantity"] = output_quantity
    
    
    if FLAG_VERBOSE:
        uvspec.inp["verbose"] = ""
    else:
        uvspec.inp["quiet"] = ""

    # provide a usefull file
    fname = "IC.DAT"
    if not os.path.isfile(fname):
        with open(fname, "w+") as f:
            f.write("#      z     LWC    R_eff\n")
            f.write("#     (km)  (g/m^3) (um) \n")
            f.write("     11.000   0      0   \n")
            f.write("     10.000   0.005  20  \n")

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
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Enter verbose (print more stuff).",
        default=False,
    )
    parser.add_argument(
        "-f",
        "--file",
        dest="filename",
        default="tatm.csv",
        help="Write results in given output file name (default: tatm.csv).",
    )
    parser.add_argument(
        "-z", "--airmass", dest="airmass", default=1.0, help="Airmass (default: 1)."
    )
    parser.add_argument(
        "-w", "--pwv", dest="pwv", default=5.0, help="PWV in mm (default: 5)."
    )
    parser.add_argument(
        "-o", "--ozone", dest="ozone", default=300.0, help="Ozone in db (default: 300)."
    )
    parser.add_argument(
        "-c",
        "--cloud",
        dest="cloud",
        default=0.0,
        help="Cloud vertical optical depth (default: 0).",
    )
    parser.add_argument(
        "-p",
        "--pressure",
        dest="pressure",
        default=1000,
        help="Pressure in hPa (default: 1000).",
    )
    parser.add_argument(
        "-s",
        "--altitudesite",
        dest="altitudesite",
        default=0.0,
        help="Site altitude in km or can be recognized observatory name as a string (default: 0).",
    )
    parser.add_argument(
        "-a",
        "--vaod",
        dest="vaod",
        default=0.0,
        help="Vertical aerosol optical depth (default: 0).",
    )
    parser.add_argument(
        "-e",
        "--exp",
        dest="exp",
        default=1.4,
        help="Angstrom exponent for aerosols (default: 1.4).",
    )
    parser.add_argument(
        "-l",
        "--wl0",
        dest="wl0",
        default=500,
        help="Reference wavelength for aerosol power law function in nm (default: 500).",
    )
    parser.add_argument(
        "-q",
        "--interactproc",
        dest="interactproc",
        default="sa",
        help="Interaction processus (default: 'sa').",
    )
    parser.add_argument(
        "-m",
        "--atmmodel",
        dest="atmmodel",
        default="us",
        help="Atmospheric model (default: 'us').",
    )
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

    wl, atm = ProcessSimulation(
        airmass=airmass_nb,
        pwv=pwv_nb,
        ozone=oz_nb,
        pressure=press_nb,
        aod=vaod_nb,
        angstrom_exponent_num=exp_nb,
        aer_lambda0=wl0_nb,
        atm_model=args.atmmodel,
        proc=args.interactproc,
        cloudext=cld_nb,
        altitude=args.altitudesite,
        FLAG_VERBOSE=FLAG_VERBOSE,
    )

    np.savetxt(args.filename, np.array([wl, atm]).T, delimiter=",")
