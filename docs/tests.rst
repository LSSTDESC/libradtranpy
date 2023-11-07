Tests of installation
=====================

Test of libRadtran installation
-------------------------------

``libRadtran`` must be installed according the instructions given at 
`libRadtran web site <http://www.libradtran.org/>`_ .


In addition ``libradtranpy`` requires that one
environnement variable **LIBRADTRANDIR** must be set to ``libRadtran`` installation path,
under which one have each of these directories:

* /bin 
* /share/libRadtran/data
* /include 
* /lib and 
* /share 

from ``libRadtran`` installation directory. 

example:

| >> tree -L 1 libRadtran
| libRadtran
| ├── bin
| ├── include
| ├── lib
| ├── libRadtran-2.0.5
| └── share


and the inside the share directory the ``data/`` directory must be available:


| >> cd libRadTran
| /> tree -L 2 share/
| share/
| └── libRadtran
|    ├── GUI
|    ├── data
|    ├── doc
|    └── examples


and inside the ``data/`` directory you must have *libRadtran data* installed as folow


| >> tree -L 1 data
| data
| ├── aerosol
| ├── albedo
| ├── altitude
| ├── atmmod
| ├── correlated_k
| ├── crs
| ├── filter
| ├── ic
| ├── nca_lookup
| ├── scripts
| ├── solar_flux
| └── wc





libradtranpy tests
------------------

After ``libRadtran and libradtranpy installation`` (using the command ``pip install -e '.[dev]'``), 
a simple test to check if libradtran is installed correctly runs:


| # call unit tests
| >> python -m unittest tests/libradtranpy/*.py


or in verbose mode:


| # call unit tests in verbose mode
| >> python -m unittest -v tests/libradtranpy/*.py


it checks:
* the existence of the environnement variable `LIBRADTRANDIR` point to the libradtran installation top directory
* the existence of the `bin` directory and the executable `uvspec`
* the existence of the `share/libRadtran/data` directory which hold all the internal data that `libRadtran`requires for its execution.   
* the execution of a libradtran simulation works well returning some data.
       

