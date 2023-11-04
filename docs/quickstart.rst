Quickstart
==========
      

Installation
------------

.. code-block:: bash

   git https://github.com/LSSTDESC/libradtranpy.git
   cd librandtranpy
   python setup.py install
   


About libradtran
----------------

Libradtran is an atmospheric transmission full simulation package which can be downloaded
from the web site http://www.libradtran.org/.

.. figure:: images/libradtran.png
   :width: 200

This emulator provides interpolations from atmospheric transmissions for scattering and absorption
processes photon-air which are calculated by libradtran (version 2.0.5 for this current release).  


 

Usage
-----

The accessas follow.
These are detailed in :doc:`apidocs`.

.. code::

   >>> from atmemulator.atmemulator import AtmEmulator
   >>> emul =  AtmEmulator()
   >>> # or
   >>> emul =  AtmEmulator('CTIO')
   >>> # or 
   >>> emul =  AtmEmulator('LSST',743.0)
   
   >>> wl = [400.,800.,900.] # define the wavelength array
   >>> am=1.2  # set the airmass
   >>> pwv =4.0  # set the precipitable water vapor in mm
   >>> oz=300. # set the ozone depth on DU
   >>> transm = emul.GetAllTransparencies(wl,am,pwv,oz)

