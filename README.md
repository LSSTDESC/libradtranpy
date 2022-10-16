# libradtranpy : libradtran python wrapper

- author : Sylvie Dagoret-Campagne
- affiliation : IJCLab/IN2P3/CNRS
- creation date : October 2022


## Installation

### Installation of libradtran

The atmospheric simulation LibRadTran must be installed according the instruction given in 
http://www.libradtran.org/doku.php


### To use libradtran inside librandtranpy wrapper


Environnement variable **LIBRADTRANDIR** must be set to libradtran installation path under which one have /bin /data /include /lib and /share of libradtran installation directory. 

ex:

	ls $LIBRADTRANDIR
	bin                     data                    include                 	lib                     libRadtran-2.0.3        	libradtran-2.0.3.tar.gz share


### Installation of libradtranpy

    cd libradtranpy
    python setup.py install
       
## Use if libradtranpy


### Use in the shell

if **libradtranpy/libradtranpy/libsimulateVisible.py** is in the python path:


	libradtranpy/libsimulateVisible.py  [-v] -z <airmass> -w <pwv> -o <oz> -a<aer> -p <P> -c <cld> -m<mod> -q<proc>
 	 - z   : airmass from 1.0 to 3.0, typical z=1 
 	 - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm
 	 - oz  : ozone in Dobson units from 200 DU to 400 DU
 	 - aer : Aerosols vertical optical depth, typical a=0.04
 	 - p   : Pressure in hPa, typical P=775.3 hPa  
 	 - c   : Cloud vertical optical depth, typical c=0
 	 - m   : Atmospheric model, typical m='us' 
 	 - q   : Interaction processes, typical q='sa' for scattering and absorption
 	 - v   : activate verbose to get atmospheric profile
	 Examples : 
	 	 1) python libsimulateVisible.py -z 1 -w 0 -o 0 -a 0
	 	 2) python libsimulateVisible.py -z 1 -w 4 -o 300 -a 0.3 -c 0 -p 742 -m 'us' -q 'sa'
	 To generate ascii printout of the used atmospheric model table in a log file :
	 	 python libsimulateVisible.py -v -z 1 -w 0 -o 0 -a 0 >& output.log
	 

### Outputs of libradtran

**Librandtran** generate output acii files consisting of rows of (wavelength, transmission).
 
**librandtranpy** manages the different simulations and their output files in a hierarchical directories. The top level directory is **simulations/**.

The output of libradtran can be found in subdirs of **simulations/RT/2.0.3/LS/pp/**.



	 	 
	 	 
### Use of libradtranpy as python package library

The call of libradtran through libradtranpy can be done as follow:

    from libradtranpy import libsimulateVisible
      
A call without aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,pressure,
                                                      prof_str='us',proc_str='sa',cloudext=cloudext)
A call with aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,aer,pressure,
                                                      prof_str='us',proc_str='sa',cloudext=cloudext)


**path,thefilename** are the path and filename of the output ascii file.

The result of the simulation can be obtaiend by:

    data = np.loadtxt(os.path.join(path,thefile))
    wl = data[:,0]
    atm = data[:,1]                                                 
                                                      
                                                      

- See the notebooks for more examples in **libradtranpy/notebooks/**



