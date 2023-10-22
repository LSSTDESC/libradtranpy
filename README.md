# libradtranpy : libradtran python wrapper

- author : Sylvie Dagoret-Campagne
- affiliation : IJCLab/IN2P3/CNRS
- creation date : October 2022
- Last update : October 22th 2023


## Installation

### Installation of libradtran

The atmospheric simulation LibRadTran fir multiple site (different altitudes) must be installed according the instruction given in 
http://www.libradtran.org/doku.php

This documentation assumes libradtran version 2.0.5 is installed on your computer (version July 2023)


### To use libradtran inside librandtranpy wrapper


Environnement variable **LIBRADTRANDIR** must be set to libradtran installation path under which one have /bin /data /include /lib and /share of libradtran installation directory. 

ex:

	ls $LIBRADTRANDIR
	bin                     data                    include                 	lib                     libRadtran-2.0.5        share


### Installation of libradtranpy

Installation of libradtranpy from configuring setuptools defined in pyproject.toml file.
(see https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)
This package is maintained through the tool LINCC Frameworks Python Project Template (https://lincc-ppt.readthedocs.io/en/latest/index.html)

    cd libradtranpy
    # on Linux
	pip install -e .[dev]
    # on Mac M1 or M2 with zsh

    pip install -e '.[dev]'
       
## Use  libradtranpy


### Use in the shell

if **libradtranpy/src/libradtranpy/libsimulateVisible.py** is in the python path:


	libradtranpy/libsimulateVisible.py  [-v] -z <airmass> -w <pwv> -o <oz> -a<aer> -p <P> -c <cld> -m<mod> -q<proc> -s<site>
 	 - z   : airmass from 1.0 to 3.0, typical z=1 
 	 - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm
 	 - oz  : ozone in Dobson units from 200 DU to 400 DU
 	 - aer : Aerosols vertical optical depth, typical a=0.04
 	 - p   : Pressure in hPa, typical P=775.3 hPa  
 	 - c   : Cloud vertical optical depth, typical c=0
 	 - m   : Atmospheric model, typical m='us' 
 	 - q   : Interaction processes, typical q='sa' for scattering and absorption
     - s   : Observation site : LSST, CTIO, ....  
 	 - v   : activate verbose to get atmospheric profile
	 Examples : 
	 	 1) python libsimulateVisible.py -z 1 -w 0 -o 0 -a 0 -s LSST
	 	 2) python libsimulateVisible.py -z 1 -w 4 -o 300 -a 0.3 -c 0 -p 742 -m  us -q sa -s LSST
	 To generate ascii printout of the used atmospheric model table in a log file :
	 	 python libsimulateVisible.py -v -z 1 -w 0 -o 0 -a 0 -s LSST >& output.log
	 

### Outputs of libradtran

**Librandtran** generate output acii files consisting of rows of (wavelength, transmission).
 
**librandtranpy** manages the different simulations and their output files in a hierarchical directories. The top level directory is **simulations/**.

The output of libradtran can be found in subdirs of **simulations/RT/2.0.5/observationsite/pp/**.

	 	 
	 	 
### Use of libradtranpy as python package library

The call of libradtran through libradtranpy can be done as follow:

    from libradtranpy import libsimulateVisible
      
A call without aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,pressure,
                                                      prof_str='us',proc_str='sa',cloudext=cloudext,altitude_str="LSST")
A call with aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,aer,pressure,
                                                      prof_str='us',proc_str='sa',cloudext=cloudext,altitude_str="LSST")


**path,thefilename** are the path and filename of the output ascii file.

The result of the simulation can be obtaiend by:

    data = np.loadtxt(os.path.join(path,thefile))
    wl = data[:,0]
    atm = data[:,1]                                                 
                                                      
                                                      




## Documentation


In [docs/notebooks/intro_notebooks.ipynb](docs/notebooks/intro_notebooks.ipynb) a series of notebooks show the use of libradtranpy and a set of tools on atmospher to control its output. 



