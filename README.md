# libradtranpy : a libradtran python wrapper


## Installation

### Installation of libradtran

The atmospheric simulation LibRadTran fir multiple site (different altitudes) must be installed according the instruction given in 
http://www.libradtran.org/doku.php

This documentation assumes libradtran version 2.0.5 is installed on your computer (version July 2023)


### To use libradtran inside librandtranpy wrapper


Environnement variable **LIBRADTRANDIR** must be set to libradtran installation path under which one have each of these directories:

- /bin 
- /share/libratran/data
- /include 
- /lib and 
- /share 

from libradtran installation directory. 

example:

/Users/dagoret/MacOSX/External>tree -L 1 libRadtran
libRadtran
├── bin
├── include
├── lib
├── libRadtran-2.0.5
└── share

and the inside the share directory the data dir must be available:

/Users/dagoret/MacOSX/External/libRadtran>tree -L 2 share/
share/
└── libRadtran
    ├── GUI
    ├── data
    ├── doc
    └── examples
	

### Installation of libradtranpy

Installation of libradtranpy from configuring setuptools defined in pyproject.toml file.
(see https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)
This package is maintained through the tool LINCC Frameworks Python Project Template (https://lincc-ppt.readthedocs.io/en/latest/index.html)

    cd libradtranpy
    # on Linux
	pip install -e .[dev]
    # on Mac M1 or M2 with zsh

    pip install -e '.[dev]'


### run tests

A simple test to check if libradtran is installed correctly runs:

    python -m unittest tests/libradtranpy/*.py

or in verbose mode:

    python -m unittest -v tests/libradtranpy/test_libsimulateVisible.py

it checks:
- the existence of the environnement variable `LIBRADTRANDIR` point to the libradtran installation top directory
- the existence of the `bin` directory and the executable `uvspec`
- the existence of the `data` directory which hold all the internal data that `libRadtran`requires for its execution.   
- the execution of a libradtran simulation works well returning some data. 
       
## Standard Use of libradtranpy

Two libRadtran running modes are available:
- visible mode from wavelength range : 250.0 nm -  1200.0 nm
- thermal mode from wavelength range : 2500 nm -  100000.0 nm

(these ranges are hardcoded, but it will be configurable in future).

Then whe have two interface modules for these modes : 
- **libradtranpy.libsimulateVisible.py** for the visible mode,
- **libradtranpy.libsimulateThermal.py** for the thermal mode.


### Use in the shell

if **libradtranpy/src/libradtranpy/libsimulateVisible.py** is in the python path:


	libradtranpy/libsimulateVisible.py  [-v] -z <airmass> -w <pwv> -o <oz> -a<aer> -p <P> -c <cld> -m<mod> -q<proc> -s<site>
 	 - z   : airmass from 1.0 to 3.0, typical z=1 
 	 - pwv : precipitable watr vapor in kg per m2 or mm, typical pwv = 5.18 mm
 	 - oz  : ozone in Dobson units from 200 DU to 400 DU
 	 - aer : Aerosols vertical optical depth, typical a=0.04
 	 - p   : Pressure in hPa, typical P=775.3 hPa, optional  
 	 - c   : Cloud vertical optical depth, optional ,typical c=0
 	 - m   : Atmospheric model, typical m='us' 
 	 - q   : Interaction processes, typical q='sa' for scattering and absorption
     - s   : Observation site : LSST, CTIO, ....  
 	 - v   : activate verbose to get atmospheric profile
	 Examples : 
	 	 1) python libsimulateVisible.py -z 1 -w 0 -o 0 -a 0 -s LSST
	 	 2) python libsimulateVisible.py -z 1 -w 4 -o 300 -a 0.3 -c 0 -p 742 -m  us -q sa -s LSST
	 To generate ascii printout of the used atmospheric model table in a log file :
	 	 python libsimulateVisible.py -v -z 1 -w 0 -o 0 -a 0 -s LSST >& output.log
	 
By example just run the following command in the shell:
     
    python libsimulateVisible.py -z 1 -w 0 -o 0 -a 0 -s LSST 


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


In [docs/notebooks/intro_notebook.ipynb](docs/notebooks/intro_notebook.ipynb) a series of notebooks show the use of libradtranpy and a set of tools on atmospher to control its output. 

- or at https://libradtranpy.readthedocs.io/en/latest/


## Dev installation Guide from LINCC-Frameworks - Getting Started with python project template


Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. LINCC-Frameworks engineers primarily use `conda` to manage virtual
environments. If you have conda installed locally, you can run the following to
create and activate a new environment.

.. code-block:: bash

   >> conda create env -n <env_name> python=3.10
   >> conda activate <env_name>


Once you have created a new environment, you can install this project for local
development using the following commands:

.. code-block:: bash

   >> pip install -e .'[dev]'
   >> pre-commit install
   >> conda install pandoc


Notes:

1) The single quotes around ``'[dev]'`` may not be required for your operating system.
2) ``pre-commit install`` will initialize pre-commit for this local repository, so
   that a set of tests will be run prior to completing a local commit. For more
   information, see the Python Project Template documentation on
   `pre-commit <https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html>`_.
3) Install ``pandoc`` allows you to verify that automatic rendering of Jupyter notebooks
   into documentation for ReadTheDocs works as expected. For more information, see
   the Python Project Template documentation on
   `Sphinx and Python Notebooks <https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks>`_.

