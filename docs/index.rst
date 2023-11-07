.. libradtranpy documentation main file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to libradtranpy's documentation!
========================================================================================

Welcome to the libradtranpy python interface to the libRadtran program.

``libradtranpy`` is a python package which has been developped to provide an easy access to the ``libRadtran`` radiative simulation package.
It is intended to provide an atmospheric transmission model that can  be used by the astronomical community for dealing with atmospheric transmission calibration issues.

This package can be downloaded from https://github.com/LSSTDESC/libradtranpy

   git clone https://github.com/LSSTDESC/libradtranpy.git


The interface can be used for the two modes of ``libRadtran``, 
* the visible mode which is usefull for ground observatories which wavelength range is covered by the u,g,r,i,z,y filters, thus in a
wavelength range from 300 nm to 1200 nm.
* the thermal mode which is usefull for observatories which monitor clouds in the Infra-Red range with an Infrared Camera.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home page <self>
   introduction
   installation
   quickstart
   quickstartther
   tests
   apidocs
   Notebooks <notebooks>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`