Quickstart
==========
      

Installation
------------

.. code-block:: bash

    git https://github.com/LSSTDESC/libradtranpy.git
    cd librandtranpy
    pip install -e .'[dev]'
   


Usage
-----



Usage in a shell script
```````````````````````

if **libradtranpy/src/libradtranpy/libsimulateVisible.py** is in the python path:


	libradtranpy/libsimulateVisible.py  
    [-v] -z <airmass> -w <pwv> -o <oz> -a<aer> -p <P> -c <cld> -m<mod> -q<proc> -s<site>

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
	 

Outputs of libradtran
~~~~~~~~~~~~~~~~~~~~~

* ``Librandtran`` generates output acii files consisting of rows of (wavelength, transmission).
 
* ``librandtranpy`` manages the different simulations and their output files in a hierarchical directories. The top level directory is **simulations/**.

The output of libradtran can be found in subdirs of 
``simulations/RT/2.0.5/observationsite/pp/``.

	 	 
	 	 
Use of libradtranpy as python package library
`````````````````````````````````````````````````

The call of libradtran through libradtranpy can be done as follow:

    from libradtranpy import libsimulateVisible
      
* A call without aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,pressure,
                prof_str='us',proc_str='sa',cloudext=cloudext,altitude_str="LSST")

* A call with aerosols:

    path,thefile=libsimulateVisible.ProcessSimulation(am[index],pwv,ozone,aer,pressure,
                prof_str='us',proc_str='sa',cloudext=cloudext,altitude_str="LSST")


* ``path,thefilename`` are the path and filename of the output ascii file.

The result of the simulation can be obtaiend by:

        data = np.loadtxt(os.path.join(path,thefile))
        wl = data[:,0]
        atm = data[:,1]                                                 
                                                      
                                                      

Remarks on the documentation on readthedocs
```````````````````````````````````````````


As ``libRadtran`` is not installed on ``readthedocs`` computer, the following example
below cannot appear.

The access is shown as follow:

These are detailed in :doc:`apidocs`.

.. code::
   >>> import os
   >>> import numpy as np
   >>> from libradtranpy import libsimulateVisible
   >>> # check libradtran is in your path
   >>> os.getenv('LIBRADTRANDIR')
   >>> am=1.2  # set the airmass
   >>> pwv =4.0  # set the precipitable water vapor in mm
   >>> oz=300. # set the ozone depth on DU
   >>> pressure = 0. # use default value
   >>> cloudext=0 # use default
   >>> path,thefile=libsimulateVisible.ProcessSimulation(am,pwv,ozone,pressure,
         prof_str='us',proc_str='sa',cloudext=cloudext,altitude_str='LSST',FLAG_VERBOSE=False)
   >>> data = np.loadtxt(os.path.join(path,thefile))
   >>> wl = data[:,0]   # wavelength array
   >>> transm = data[:,1] # transmission array

Moreover the ``libradtranpy.libsimulateVisible`` API cannot appear on readthedocs website.
This API may appear on user's computer if ``libRadtran`` is installed correctly. 