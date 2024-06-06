"""
Module libratranpy, an interface package to libradtran executable
library to simulate air transparency with LibRadTran

:author: sylvielsstfr
:creation date: November 2nd 2016
:last update: November 7th 2023
"""

import os
import re
import numpy as np
import sys
import argparse
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

Prog = "RT"  # definition the simulation programm is libRadTran
Obs = "LS"  # definition of observatory site (LS,CT,OH,MK,...)
Rte = "pp"  # pp for parallel plane of ps for pseudo-spherical
Atm = ["us"]  # short name of atmospheric sky here US standard and  Subarctic winter
Proc = (
    "sa"  # light interaction processes : sc for pure scattering,ab for pure absorption
)
# sa for scattering and absorption, ae with aerosols default, as with aerosol special
Mod = "rtthermal"  # Models for absorption bands : rt for REPTRAN, lt for LOWTRAN, k2 for Kato2
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
List_Of_Thermal_Outputs = [
    "BRIGHTNESS",
    "RADIANCE",
    "IRRADIANCE",
    "IRRADIANCE_INTEGRATED",
    "TRANSMITTANCE",
]

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


#########################################################################


def usage():
    """Description of the usageof the script"""
    print("*******************************************************************")
    print(
        sys.argv[0],
        " [-v] -z <airmass> -w <pwv> -o <oz> -p <P> -c <cld> -m<mod> -q<proc> -s<site-string> -v<verbose-flag>",
    )
    print(" \t - z   : airmass from 1.0 to 3.0, typical z=1 ")
    print(
        " \t - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm"
    )
    print(" \t - oz  : ozone in Dobson units from 200 DU to 400 DU")
    print(" \t - p   : Pressure in hPa, typical P=775.3 hPa  ")
    print(" \t - c   : Cloud vertical optical depth, typical c=0")
    print(" \t - m   : Atmospheric model, typical m='us' ")
    print(
        " \t - q   : Interaction processes, typical q='sa' for scattering and absorption"
    )
    print(
        " \t - s   : provide site or altitude as a string : LSST/OHP/PDM/OMK/OSL or altitude in km like akm_2.663"
    )
    print(" \t - v   : activate libradtran output verbose mode")
    print(
        " \t - t   : thermal output selected among (brightness,radiance,irradiance,irradiance_integrated,transmisstance). Optional. Default brightness"
    )

    print("\t Examples : ")
    print("\t \t 1) python libsimulateThermal.py -z 1 -w 0 -o 0 -s LSST")
    print(
        "\t \t 2) python libsimulateThermal.py -z 1 -w 4 -o 300 -c 0 -p 742 -m us -q sa -s LSST -t irradiance"
    )

    print(
        "\t To generate ascii printout of the used atmospheric model tables in a log file :"
    )
    print("\t \t python libsimulateThermal.py -v -z 1 -w 0 -o 0 -s LSST >& output.log")

    print("\t Actually provided : ")
    print("\t \t Number of arguments:", len(sys.argv), "arguments.")
    print("\t \t Argument List:", str(sys.argv))

    print("*******************************************************************")



def ProcessSimulation(
    airmass_num,
    pwv_num,
    oz_num,
    press_num,
    prof_str="us",
    proc_str="sa",
    cloudext=0.0,
    albedo=0.0,
    temperature=273.15,
    altitude="LSST",
    thermal_output="brightness",
    wl_min_nm=5000,
    wl_max_nm=20000,
    molresol="coarse",
    FLAG_VERBOSE=False,
):
    """
    ProcessSimulation(airmass_num,pwv_num,oz_num)
    Function to simulate air transparency.
    No aerosol simulation is performed.

    Parameters:

    :param airmass_num: airmass, useless in thermal mode
    :type airmass_num: float, unitless

    :param pwv_num: precipitable water vapor
    :type pwv_num: float, in units mm

    :param oz_num: ozone colon depth
    :type oz_num: float, in Dobson unit

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

    :param altitude: observation site predefined (either the site name abrebiation (LSST/CTIO,OMP,OMK,MSL) or the string on altitude like akm_2.663)
    :type altitude: string

    :param FLAG_VERBOSE: flag to activate libradtran verbose mode that print all the table generated
    :type FLAG_VERBOSE: bool

    :param thermal_output: select the kind of output for thermal mode, default set to brightness
    :type thermal_output: string among (brightness, irradiance, integrated_irradiance, radiance, transmittance)

    :returns: OUTPUTDIR,outputFilename, path and filename of datafile containing the simulated data
    :rtype: two strings



    """

    FLAG_BRIGHTNESS = False
    FLAG_IRRADIANCE = False
    FLAG_IRRADIANCE_INTEGRATED = False
    FLAG_RADIANCE = False
    FLAG_TRANSMITTANCE = False

    if FLAG_VERBOSE:
        print(
            "--------------- ProcessSimulation input args -----------------------------"
        )
        print(" 1) airmass = ", airmass_num)
        print(" 2) pwv = ", pwv_num)
        print(" 3) oz = ", oz_num)
        print(" 4) pressure  = ", press_num)
        print(" 5) atmospheric profile = ", prof_str)
        print(" 6) interaction processes = ", proc_str)
        print(" 7) cloud extinction = ", cloudext)
        print(" 8) site or altitude = ", altitude)
        print(" 9) thermal output required", thermal_output)
        print("--------------------------------------------")

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

    if FLAG_VERBOSE:
        print(f"Observation site altitude for libradran sim : {altitude_num} km")

    # set thermal mode

    thermal_output_uppercase = thermal_output.upper()

    if thermal_output_uppercase == List_Of_Thermal_Outputs[0]:
        FLAG_BRIGHTNESS = True
    elif thermal_output_uppercase == List_Of_Thermal_Outputs[1]:
        FLAG_RADIANCE = True
    elif thermal_output_uppercase == List_Of_Thermal_Outputs[2]:
        FLAG_IRRADIANCE = True
    elif thermal_output_uppercase == List_Of_Thermal_Outputs[3]:
        FLAG_IRRADIANCE_INTEGRATED = True
    elif thermal_output_uppercase == List_Of_Thermal_Outputs[4]:
        FLAG_TRANSMITTANCE = True
    else:
        msg = f"Bad thermal output selection {thermal_output}"
        raise Exception(msg)

    if not (
        FLAG_BRIGHTNESS
        or FLAG_IRRADIANCE
        or FLAG_IRRADIANCE_INTEGRATED
        or FLAG_RADIANCE
        or FLAG_TRANSMITTANCE
    ):
        msg = f"No thermal output selection provided {thermal_output}"
        raise Exception(msg)

    # set the interaction process

    # Proc='sa'  # Pure absorption and Rayleigh scattering : Clear sky without aerosols
    if proc_str in ["sa", "ab", "sc", "ae", "as"]:
        Proc = proc_str

    # set the selected atmosphere
    if prof_str in ["us", "ms", "mw", "tp", "ss", "sw", "ohp_winter",
        "ohp_spring",
        "ohp_summer",
        "ohp_autumn",]:
        Atm = [prof_str]

    # create output dir
    TOPDIR = os.path.join(TOPTOPDIR, str(altitude_num))
    ensure_dir(TOPDIR)

    # build the part 1 of filename
    BaseFilename_part1 = Prog + "_" + Obs + "_" + Rte + "_"

    # Set up type of

    runtype = (
        "clearsky"  #'no_scattering' #aerosol_special #aerosol_default# #'clearsky'#
    )

    if Proc == "sc":
        runtype = "no_absorption"
        outtext = "no_absorption"
    elif Proc == "ab":
        runtype = "no_scattering"
        outtext = "no_scattering"
    elif Proc == "sa":
        runtype == "clearsky"
        outtext = "clearsky"
    elif Proc == "ae":
        runtype = "aerosol_default"
        outtext = "aerosol_default"
    elif Proc == "as":
        runtype = "aerosol_special"
        outtext = "aerosol_special"
    else:
        runtype == "clearsky"
        outtext = "clearsky"

    #   Selection of RTE equation solver
    if Rte == "pp":  # parallel plan
        rte_eq = "disort"
    elif Rte == "ps":  # pseudo spherical
        rte_eq = "sdisort"
    else:
        raise ValueError(f"Unknown RTE equation solver {Rte=}.")

    #   Selection of absorption model
    molmodel = "reptran"
    if Mod == "rt":
        molmodel = "reptran"
    if Mod == "lt":
        molmodel = "lowtran"
    if Mod == "kt":
        molmodel = "kato"
    if Mod == "k2":
        molmodel = "kato2"
    if Mod == "fu":
        molmodel = "fu"
    if Mod == "cr":
        molmodel = "crs"

    # for simulation select only two atmosphere
    # theatmospheres = np.array(['afglus','afglms','afglmw','afglt','afglss','afglsw'])
    atmosphere_map = dict()  # map atmospheric names to short names
    atmosphere_map["afglus"] = "us"
    atmosphere_map["afglms"] = "ms"
    atmosphere_map["afglmw"] = "mw"
    atmosphere_map["afglt"] = "tp"
    atmosphere_map["afglss"] = "ss"
    atmosphere_map["afglsw"] = "sw"

    theatmospheres = []

    for skyindex in Atm:
        if re.search("us", skyindex):
            theatmospheres.append("afglus")
        if re.search("sw", skyindex):
            theatmospheres.append("afglsw")
        if re.search("ss", skyindex):
            theatmospheres.append("afglss")
        if re.search("mw", skyindex):
            theatmospheres.append("afglmw")
        if re.search("ms", skyindex):
            theatmospheres.append("afglms")
        if re.search("tp", skyindex):
            theatmospheres.append("afglt")
        if re.search("ohp_winter", skyindex):
            theatmospheres.append("ohp_winter")
        if re.search("ohp_spring", skyindex):
            theatmospheres.append("ohp_spring")
        if re.search("ohp_summer", skyindex):
            theatmospheres.append("ohp_summer")
        if re.search("ohp_autumn", skyindex):
            theatmospheres.append("ohp_autumn")

    # 1) LOOP ON ATMOSPHERE
    for atmosphere in theatmospheres:
        atmkey = atmosphere_map[atmosphere]

        # manage input and output directories and vary the ozone
        TOPDIR2 = TOPDIR + "/" + Rte + "/" + atmkey + "/" + Proc + "/" + Mod
        ensure_dir(TOPDIR2)
        INPUTDIR = TOPDIR2 + "/" + "in"
        ensure_dir(INPUTDIR)
        OUTPUTDIR = TOPDIR2 + "/" + "out"
        ensure_dir(OUTPUTDIR)

        # water vapor
        pwv_val = pwv_num
        wvfileindex = int(10 * pwv_val)

        # airmass
        airmass = airmass_num
        amfileindex = int(airmass_num * 10)

        # Ozone
        ozfileindex = int(oz_num / 10.0)

        # Cloud
        cldindex = str(int(cloudext * 1000))
        cld_str = cldindex.zfill(4)

        # root path for filename
        BaseFilename = (
            BaseFilename_part1
            + atmkey
            + "_"
            + Proc
            + "_"
            + Mod
            + "_z"
            + str(amfileindex)
            + "_"
            + WVXX
            + str(wvfileindex)
            + "_"
            + OZXX
            + str(ozfileindex)
            + "_"
            + CLD
            + cld_str
        )

        # verbose=True
        verbose = FLAG_DEBUG

        # premare to create libradtran input file

        uvspec = UVspec3.UVspec()
        uvspec.inp["data_files_path"] = os.path.join(
            os.path.expanduser("~"), libradtrandatapath
        )

        uvspec.inp["atmosphere_file"] = os.path.join(
            libradtrandatapath, "atmmod/" + atmosphere + ".dat"
        )

        # arbitrary earth albedo
        uvspec.inp["albedo"] = str(albedo)

        uvspec.inp["rte_solver"] = rte_eq

        if Mod == "rtthermal":
            uvspec.inp["mol_abs_param"] = molmodel + " " + molresol
        else:
            uvspec.inp["mol_abs_param"] = molmodel

        # Convert airmass into zenith angle (could be inproved)
        sza = np.arccos(1.0 / airmass_num) * 180.0 / np.pi

        # water vapor
        uvspec.inp["mol_modify H2O"] = f"{pwv_val:.20f} MM"
        # Ozone
        uvspec.inp["mol_modify O3"] = (
            f"{oz_num:.20f} DU"  # https://ozonewatch.gsfc.nasa.gov/NH.html
        )

        # rescale pressure   if reasonable pressure values are provided
        if press_num > 200.0 and press_num < 1080.0:
            uvspec.inp["pressure"] = press_num

        uvspec.inp["ic_file"] = "1D ./IC.DAT"
        uvspec.inp["ic_properties"] = "yang"
        uvspec.inp["ic_modify"] = "tau set " + str(cloudext)
        uvspec.inp["altitude"] = altitude_num  # Altitude  observatory

        uvspec.inp["sur_temperature"] = str(temperature)

        # in thermal mode
        uvspec.inp["source"] = "thermal "
        uvspec.inp["wavelength"] = f"{wl_min_nm} {wl_max_nm}"

        if FLAG_BRIGHTNESS:
            uvspec.inp["output_user"] = "lambda edn"
            uvspec.inp["output_quantity"] = "brightness"
        elif FLAG_IRRADIANCE:
            uvspec.inp["output_user"] = "lambda edn"
            uvspec.inp["output_process"] = "per_nm"
        elif FLAG_IRRADIANCE_INTEGRATED:
            uvspec.inp["output_user"] = "lambda edn"
            uvspec.inp["output_process"] = "sum"
        elif FLAG_RADIANCE:
            zenith_angle = np.rad2deg(np.arccos(1 / airmass))
            # looking downward = umu > 0
            # looking upward = umu < 0
            uvspec.inp["umu"] = (
                str(-np.cos(np.deg2rad(zenith_angle)))
                .replace("[", "")
                .replace("]", "")
                .replace("\n", "")
            )
            uvspec.inp["output_user"] = "lambda uu"
            uvspec.inp["output_process"] = "per_nm"
        elif FLAG_TRANSMITTANCE:
            uvspec.inp["output_user"] = "lambda edn"
            uvspec.inp["output_quantity"] = "reflectivity"

        if FLAG_VERBOSE:
            uvspec.inp["verbose"] = ""
        else:
            uvspec.inp["quiet"] = ""

        if "output_quantity" in uvspec.inp.keys():
            outtextfinal = outtext + "_" + uvspec.inp["output_quantity"]

        inputFilename = BaseFilename + ".INP"
        outputFilename = BaseFilename + ".OUT"
        inp = os.path.join(INPUTDIR, inputFilename)
        out = os.path.join(OUTPUTDIR, outputFilename)

        # provide a usefull file
        fname = "IC.DAT"
        if not os.path.isfile(fname):
            with open(fname, "w+") as f:
                f.write("#      z     LWC    R_eff\n")
                f.write("#     (km)  (g/m^3) (um) \n")
                f.write("     11.000   0      0   \n")
                f.write("     10.000   0.005  20  \n")

        # really write libradtran input file
        uvspec.write_input(inp)
        # run libradtran
        simu = uvspec.run(FLAG_VERBOSE, path=libradtranpath)
        wl, atm = simu[0], simu[1]

    # return path of output file
    return wl, atm


# ---------------------------------------------------------------------------


######################################################################################
# The main program for starting atmospheric simulation start here
#####################################################################################

if __name__ == "__main__":

    # Create the parser
    parser = argparse.ArgumentParser(description="Process some parameters.")

    # Define the arguments
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output",
        default=False,
    )
    parser.add_argument(
        "-f",
        "--file",
        dest="filename",
        default="latm.csv",
        help="Write results in given output file name (default: latm.csv).",
    )
    parser.add_argument(
        "-z", "--airmass", dest="airmass", default=1.0, help="Airmass value"
    )
    parser.add_argument("-w", "--pwv", dest="pwv", default=5.0, help="PWV value")
    parser.add_argument(
        "-o", "--ozone", dest="ozone", default=300.0, help="Ozone value in DU"
    )
    parser.add_argument(
        "-c", "--cloudd", dest="cloud", default=0.0, help="Cloud VOD value"
    )
    parser.add_argument(
        "-p", "--pressure", dest="pressure", default=1000, help="Pressure value in hPa"
    )
    parser.add_argument(
        "-s",
        "--site",
        dest="altitudesite",
        default=0.0,
        help="Site altitude in km or can be recognized observatory name as a string (default: 0).",
    )

    parser.add_argument(
        "-m", "--atmmodel", dest="atmmodel", default="us", help="Atmospheric model"
    )
    parser.add_argument(
        "-q", "--interactproc", dest="interactproc", default="clearsky", help="Interaction process"
    )

    parser.add_argument(
        "-t", "--thermal", dest="thermal", default="radiance", help="Thermal response output (radiance, brightness...)"
    )

    # Parse the arguments
    args = parser.parse_args()
    
    
    airmass_nb = float(args.airmass)
    pwv_nb = float(args.pwv)
    oz_nb = float(args.ozone)
    press_nb = float(args.pressure)
    cld_nb = float(args.cloud)
    model_str = str(args.atmmodel)
    altitude = args.altitudesite
    FLAG_VERBOSE = args.verbose
    interactproc = args.interactproc
    therm_str = args.thermal

    # Debug print statements (uncomment for debugging)
    # print("Parsed arguments:", args)

    
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

    if therm_str.upper() not in List_Of_Thermal_Outputs:
        print("bad thermal response selected : thermal response = ", therm_str)
        print("among :", List_Of_Thermal_Outputs)
        sys.exit()

    wl, atm = ProcessSimulation(
        airmass_num=airmass_nb,
        pwv_num=pwv_nb,
        oz_num=oz_nb,
        press_num=press_nb,
        prof_str=model_str,
        proc_str=interactproc,
        cloudext=cld_nb,
        altitude=altitude,
        FLAG_VERBOSE=FLAG_VERBOSE,
        thermal_output=therm_str,
    )
    
    np.savetxt(args.filename, np.array([wl, atm]).T, delimiter=",")