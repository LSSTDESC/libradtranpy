Quickstart
==========
      

Installation
------------

.. code-block:: bash

   git clone git@github.com:sylvielsstfr/AtmEmulator.git
   cd AtmEmulator
   python setup.py install
   


About libradtran
----------------

Libradtran is an atmospheric transmission full simulation package which can be downloaded
from the web site http://www.libradtran.org/.

.. figure:: images/libradtran.png
   :width: 200

This emulator provides interpolations from atmospheric transmissions for scattering and absorption
processes photon-air which are calculated by libradtran (version 2.0.5 for this current release).  


Details
-------

The architecture of the ``AtmEmulator`` module relies on the implementation of 3 classes
with inheritage relationship.
The base class ``ObsAtmoGrid`` implements datafile IO and multi-dimentional interpolation
of transparencies.


.. figure:: images/inheritancediagram.png
    :width: 400

#.. inheritance-diagram:: getObsAtmo
#   :top-classes: getObsAtmo.getObsAtmo.ObsAtmoGrid, getObsAtmo.getObsAtmo.ObsAtmoPressure, getObsAtmo.getObsAtmo.ObsAtmo

.. tabs::

   .. tab::  ``Class AtmEmulator``

      Top-leve User-Endpoint class to call the emulator. This class inherit from the implemented class ``ObsAtmoPressure``.
      This interface class should remain as stable as possible given it does not implement behaviors than can
      be improved in later versions.

   .. tab:: ``Class AtmPressEmulator``

      This class inherit from the base class ``SimpleAtmEmulator``.
      This class implement corrections due to pressure effect on atmospheric transmissions.

   .. tab:: ``Class SimpleAtmEmulator``

      This is the base class which implement the interpolation of atmospheric transmissions from a series
      of dataset extracted from libradtran (http://www.libradtran.org/).
 

Usage
-----

The access to the emulator is done as follow.
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

