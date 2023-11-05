.. libradtranpy documentation main file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to libradtranpy's documentation!
========================================================================================

Welcome to the libradtran interface.
This package is designed to provide a user interface for atmospheric transmission models in the Libradtran program.

This package is intended to provide an easy configuration through python wrapper.
to libradtran program by giving a model of atmospheric transparency for astronomical observatories.
It is particularly suited to studies on atmospheric transmission.

It's important to stress that not all of libradtran's computational capabilities are 
available through this interface. Only those functions that are useful for atmospheric calibration in the wavelength range corresponding to the u,g,r,i,z,y filters are offered.

It restricts the wavelength to range covered by the u,g,r,i,z,y filters, thus in a
wavelength range from 300 nm to 1200 nm.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home page <self>
   introduction
   installation
   quickstart
   apidocs
   Notebooks <notebooks>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`