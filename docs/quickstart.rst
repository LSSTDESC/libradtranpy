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

The accessas follow.
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

