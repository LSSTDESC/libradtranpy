Tests of installation
=====================

Test of libRadtran installation
-------------------------------



The atmospheric simulation LibRadTran for multiple site (different altitudes) must be installed according the instruction given in 
http://www.libradtran.org/doku.php

This documentation assumes libradtran version 2.0.5 is installed on your computer (version July 2023)


#### To use libradtran inside librandtranpy wrapper


Environnement variable **LIBRADTRANDIR** must be set to ``libRadtran`` installation path under which one have each of these directories:

- /bin 
- /share/libRadtran/data
- /include 
- /lib and 
- /share 

from ``libRadtran`` installation directory. 

example:

```bash
>> tree -L 1 libRadtran
libRadtran
├── bin
├── include
├── lib
├── libRadtran-2.0.5
└── share
```

and the inside the share directory the ``data/`` directory must be available:

```bash
>> cd libRadTran
/> tree -L 2 share/
share/
└── libRadtran
    ├── GUI
    ├── data
    ├── doc
    └── examples
```

and inside the ``data/`` directory you must have *libRadtran data* installed as folow

```bash
>> tree -L 1 data
data
├── aerosol
├── albedo
├── altitude
├── atmmod
├── correlated_k
├── crs
├── filter
├── ic
├── nca_lookup
├── scripts
├── solar_flux
└── wc
```




libradtranpy tests
------------------

After ``libRadtran and libradtranpy installation``, a simple test to check if libradtran is installed correctly runs:

```bash
# call unit tests
>> python -m unittest tests/libradtranpy/*.py
```

or in verbose mode:

```bash
# call unit tests in verbose mode
>> python -m unittest -v tests/libradtranpy/*.py
```

it checks:
- the existence of the environnement variable `LIBRADTRANDIR` point to the libradtran installation top directory
- the existence of the `bin` directory and the executable `uvspec`
- the existence of the `share/libRadtran/data` directory which hold all the internal data that `libRadtran`requires for its execution.   
- the execution of a libradtran simulation works well returning some data. 
       

