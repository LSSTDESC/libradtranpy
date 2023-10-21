#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
# Filename:          AtmosphericExtinction.py
# Copyright          2012, Clément Buton, Yannick Copin
# Author:            Clément Buton <buton@physik.uni-bonn.de>
# Author:            $Author: ycopin $
# Version:           $Revision: 1.11 $
# Modified at:       $Date: 2017/01/05 10:28:59 $
# $Id: AtmosphericExtinction.py,v 1.11 2017/01/05 10:28:59 ycopin Exp $
###############################################################################

"""
.. _module:

AtmosphericExtinction (module)
==============================
"""

import os
import numpy as N
import astropy.io.fits as F

__author__ = "Yannick Copin <y.copin@ipnl.in2p3.fr>, " \
             "Clément Buton <c.buton@ipnl.in2p3.fr>"
__version__ = '$Id: AtmosphericExtinction.py,v 1.11 2017/01/05 10:28:59 ycopin Exp $'

# Data ==============================

# Default ozone template
O3Template = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'o3data/ozoneTemplate.fits')

EXT2OPT = .92103403719761834  # LOG10/2.5 = Extinction to opt. thickness

# Classes ======================================================================


class ExtinctionModel(object):

    def __init__(self, lbda=None, ozoneTemplate=None, lrefAero=1e4):
        """
        Extinction model, from:

        :param lbda: wavelength vector [AA] (default to extended optical range)

        :param ozoneTemplate: name of the ozone template table (see
          :func:`readOzoneTemplate`). By default, use the provided
          ozone template.
        :param lrefAero: aerosol reference wavelength [AA]
        """

        if lbda is None:
            self.lbda = N.arange(3200, 10001, 10, dtype=float)
        else:
            self.lbda = N.asanyarray(lbda)  # Wavelength [AA]

        # Rayleigh extinction template [mag/airmass] for a pressure of 1 mbar
        self.rayleigh = self.rayleigh_HT74(self.lbda, 1.)

        # Ozone extinction [mag/airmass]
        self.ozoneName = ozoneTemplate if ozoneTemplate else O3Template

        # Read transmission ozone template, interpolate at input
        # wavelengthes, and convert to extinction [mag/airmass] for an
        # ozone column density of 1 DU
        self.ozone, self.ozoneNorm = readOzoneTemplate(self.ozoneName, self.lbda)
        self.ozone /= self.ozoneNorm

        # Aerosols
        self.lrefAero = lrefAero        # Aerosol reference wavelength [AA]
        self.lbdaN = self.lbda / self.lrefAero

    def __str__(self):

        s = """\
Wavelength domain: %(m).1f-%(M).1f A by step of %(s).1f A (%(n)d px)
Ozone template: %(template)s (%(column)s DU)
Aerosol reference wavelength: %(l).0f A
""" % dict(m=self.lbda[0], M=self.lbda[-1],
           s=self.lbda[1]-self.lbda[0], n=len(self.lbda),
           template=self.ozoneName, column=self.ozoneNorm,
           l=self.lrefAero)

        if hasattr(self, 'p'):
            p, o3, tau, ang = self.p
            dp, do3, dtau, dang = self.dp
            s += """\
Input extinction parameters:
  Pressure: %(p).0f +/- %(dp).0f mbar
  Ozone:    %(o).0f +/- %(do).0f DU
  Aerosols: optical depth @ refLbda: %(t).2g +/- %(dt).2g
            angstrom exponent:         %(a).2f +/- %(da).2f
""" % dict(p=p, dp=dp, o=o3, do=do3, t=tau, dt=dtau, a=ang, da=dang)

        else:
            s += """\
Input extinction parameters: not set yet
"""

        return s

    def setParams(self, pars, dpars=None):
        """
        Set physical extinction parameters: pressure, ozone column
        density [Dobson units], aerosol optical depth at reference
        wavelength and aerosol angstrom exponent.

        :param pars: extinction parameters (*pr*, *oi*, *ai*, *ap*)
          where:

          - *pr*: surface pressure [mbar]
          - *oi*: ozone intensity [Dobson units]
          - *ai*: aerosol optical depth at reference wavelength
          - *ap*: aerosol angstrom exponent

        :param dpars: associated standard errors

        The total atmospheric extinction will then be the sum of
        three components:

        - Rayleigh extinction: `pr[mbar] * HT74(1 mbar)`
        - Ozone extinction: `oi[DU] * OzoneTemplate(1 DU)`
        - Aerosol extinction: `ai/EXT2OPT * (lbda/lRef)**(-ap)`

        .. Note:: if `ndim(dpars)==2`, `dpars` is considered as the
           *covariance* matrix of input extinction
           parameters. Therefore, when `ndim(dpars)==1`,
           `self.extinctionErrors(pars, dpars)` is equivalent to
           `self.extinctionErrors(pars, N.diag(dpars)**2)` (note the
           square power).
        """

        self.p = N.asanyarray(pars)
        if dpars is None:
            self.dp = N.zeros(4)
        else:
            self.dp = N.asanyarray(dpars)

    def setDefaultParams(self, location='Mauna Kea'):
        """
        Set default physical extinction parameters from predefined location.

        :param location: predefined location.

        =================================  =================
        Parameter                          Value ± Error
        =================================  =================
        *Mauna Kea*
        ----------------------------------------------------
        Pressure                           616 ± 2 mbar
        Ozone column                       257 ± 23 DU
        Aerosols optical depth @ 1 micron  0.0076 ± 0.0014
        Aerosols angstrom exponent         1.26 ± 1.33
        =================================  =================
        """

        if location == "Mauna Kea":
            o3, do3 = 257., 23.         # Ozone column density [DU]
            ang, dang = 1.26, 1.33      # Ångström exponent
            tau, dtau = 7.6e-3, 1.4e-3  # Aerosol optical depth at 1 micron
            p, dp = 616., 2.            # Surface pressure [mbar]
        else:
            raise ValueError("Unknown location '%s'" % location)

        self.setParams([ p, o3, tau, ang], dpars=[dp, do3, dtau, dang])

    def extinctionComponents(self):
        """
        Compute extinction individual components from extinction
        parameters (see :meth:`setParams`)

        :return: extinction components 2D-array [rayleigh,ozone,aerosols]
        """

        return N.array([
            self.p[0] * self.rayleigh,                       # Rayleigh component
            self.p[1] * self.ozone,                          # Ozone component
            self.p[2] / EXT2OPT * self.lbdaN**(-self.p[3]),  # Aerosols component
        ])

    def extinctionErrors(self):
        """
        Compute total extinction (diagonal) standard error from extinction
        parameters and associated standard errors (see :meth:`setParams`)

        :return: total extinction standard error
        """

        jac = self.jac()
        if N.ndim(self.dp) == 1:    # dp is a vector of std (independant) errors
            vExt = N.dot(self.dp**2, jac**2)
        elif N.ndim(self.dp) == 2:  # dp is actually a covariance matrix
            vExt = N.dot(N.dot(jac.T, self.dp), jac).diagonal()

        return N.sqrt(vExt)

    def extinction(self, pars=None, dpars=None, components=False):
        """
        Compute total extinction (and associated standard error)
        from extinction parameters (and associated standard errors).

        :param pars: extinction parameters (see :meth:`setParams`)
        :param dpars: extinction parameter errors (see :meth:`setParams`)
        :param components: return individual extinction components if True
        :return: 2D-array [ext,dext,[components]]
        """

        if None not in (pars, dpars):
            self.setParams(pars, dpars)

        comp = self.extinctionComponents()    # (ncomp,nlbda)
        ext = comp.sum(axis=0)                # (nlbda,)
        dext = self.extinctionErrors()        # (nlbda,)

        if not components:      # Return [lbda,ext,dext]
            return N.vstack((ext, dext))
        else:                   # Return individual components as well
            return N.vstack((ext, dext, comp))

    def jac(self):
        """Jacobian of total extinction with respect to extinction
        parameters.

        :return: jacobian 2D-array (nparam=4,nlbda)
        """

        jac = N.empty((len(self.p), len(self.lbda)), 'd')
        jac[0] = self.rayleigh                            # dext/dP
        jac[1] = self.ozone                               # dext/do3
        jac[2] = self.lbdaN**(-self.p[3]) / EXT2OPT       # dext/dtau
        jac[3] = -self.p[2] * jac[2] * N.log(self.lbdaN)  # dext/dang

        return jac

    @staticmethod
    def rayleigh_HT74(lbda, pressure):
        """
        Rayleigh extinction from `Hansen & Travis (1974)
        <http://cdsads.u-strasbg.fr/abs/1974SSRv...16..527H>`_.

        :param lbda: wavelength vector [AA]
        :param pressure: effective surface pressure [mbar]
        :return: Rayleigh extinction [mag/airmass]
        """

        lm = lbda * 1e-4                # Wavelength from A to microns

        # Optical depth
        tau = 0.008569 / lm**4 * (1 + 0.0113 / lm**2 + 0.00013 / lm**4)
        tau *= pressure / 1013.25

        return tau / EXT2OPT    # Convert to attenuation [mag/airmass]

    def write(self, outname, ext=None, format='txt'):
        """
        Write extinction curve in output file.

        :param outname: output filename
        :param ext: explicit extinction curve(s) to be written out
        :param format: output file format ('txt' or 'fits')
        """

        if ext is None:
            if hasattr(self, 'p'):
                ext = self.extinction(components=True)
            else:
                raise RuntimeError("Cannot evaluate extinction "
                                   "without extinction parameters")

        if format == 'txt':             # ASCII table

            ext = N.absolute(ext.round(6))  # Avoid rounding imprecisions

            outFile = open(outname, 'w')
            # Header
            outFile.write('\n# '.join([''] + str(self).split('\n')) + '\n')
            outFile.write('# Reference: Buton et al. (2013A&A...549A...8B)\n')
            outFile.write('# Wavelength in AA\n')
            outFile.write('# Extinctions in mag/airmass\n')
            outFile.write('# lbda   Ext   dExt    Ray    O3     Aero \n')
            # Values
            for l, e, de, r, o, a in zip(self.lbda,
                                         ext[0], ext[1], ext[2], ext[3], ext[4]):
                outFile.write(' %5d  %.3f  %.3f  %.3f  %.3f  %.3f\n' %
                              (l, e, de, r, o, a))
            outFile.close()

        elif format == 'fits':       # FITS table

            p, o3, tau, ang = self.p
            dp, do3, dtau, dang = self.dp

            keywords = [
                # Generic keywords
                ('EXTMODEL', "Rayleigh+Ozone+Aerosols", "Extinction model"),
                ('EXTREF',
                 "Buton et al. (2013A&A...549A...8B)", "Bibliographical ref."),
                # Extinction parameters and errors
                ('RA_P', p, "Surface pressure [mbar]"),
                ('RA_DP', dp, "Pressure stddev [mbar]"),
                ('OZ_INT', o3, "Ozone intensity [DU]"),
                ('OZ_DINT', do3, "Ozone intensity stddev [DU]"),
                ('AE_TAU', tau, "Aerosol optical depth"),
                ('AE_DTAU', dtau, "Aerosol optical depth stddev"),
                ('AE_ANG', ang, "Aerosol Angstrom exponent"),
                ('AE_DANG', dang, "Aerosol Angstrom exponent stddev"),
                ('AE_LREF', self.lrefAero, "Aerosol ref. wavelength [AA]"),
            ]

            # Extinction table
            arrays = [self.lbda, ext[0], ext[1], ext[2], ext[3], ext[4]]
            names = ['LAMBDA', 'EXT', 'DEXT', 'RAYLEIGH', 'OZONE', 'AEROSOLS']
            units = ['Angstrom'] + ['mag/airmass']*5

            # Extinction BinTableHDU, with keywords
            table = createTable(arrays, names, units=units, keywords=keywords,
                                extname='EXTINCTION')
            table.writeto(outname, clobber=True)

        else:
            raise IOError("Unknown output format '%s'" % format)

    def plot(self, ext=None, ax=None, components=True, transmission=False):
        """
        Plot the atmospheric extinction/transmission and its physical components.

        :param ext: extinctions to be plotted (or None)
        :param ax: matplotlib Axes instance (or None)
        :param components: display individual components if True
        :param transmission: display transmission rather than extinction
        :return: matplotlib Axes instance
        """

        if ext is None:
            if hasattr(self, 'p'):
                ext = self.extinction(components=True)
            else:
                raise RuntimeError("Cannot evaluate extinction "
                                   "without extinction parameters")

        p, o3, tau, ang = self.p

        # Non-default colors
        blue, red, green, orange = ('#0066CC', '#CC0033', '#009966', '#FF9900')

        if transmission:  # Extinction [mag/airmass] -> Transmission
            title = "SNfactory atmospheric transmission"
            ylbl = "Transmission"
            ext[0] = 10**(-0.4 * ext[0])    # Total transmission
            ext[1] *= -EXT2OPT * ext[0]     # Error on total transmission
            ext[2:] = 10**(-0.4 * ext[2:])  # Component transmissions
        else:
            title = "SNfactory atmospheric extinction"
            ylbl = "Extinction [mag/airmass]"
        title += " (Buton et al., 2013A&A...549A...8B)"

        if ax is None:                      # Create a default axes
            import matplotlib.pyplot as plt

            fig = plt.figure(figsize=(8, 5))
            ax = fig.add_subplot(1, 1, 1,
                                 title=title,
                                 xlabel=u"Wavelength [Å]",
                                 xlim=(self.lbda[0], self.lbda[-1]),
                                 ylabel=ylbl)
            # ax.ticklabel_format(style='plain')

        # Total extinction and errorband
        ax.plot(self.lbda, ext[0], color=green, lw=2, label='Total')
        #xp, yp = plt.mlab.poly_between(self.lbda, ext[0] - ext[1], ext[0] + ext[1])
        #ax.fill(xp, yp, alpha=0.3, fc=green, ec=green, label='_')

        if components:                      # Physical components
            ax.plot(self.lbda, ext[2], color=red, ls='--',
                    label='Rayleigh [%.0f mbar]' % p)
            ax.plot(self.lbda, ext[3], color=blue, ls=':',
                    label='Ozone [%.0f DU]' % o3)
            ax.plot(self.lbda, ext[4], color=orange, ls='-.',
                    label=u'Aerosols [τ=%.4f, å=%.2f]' % (tau, ang))

        ax.legend(loc='best', frameon=False)

        return ax

# Functions ==================================================================


def readOzoneTemplate(ozoneName, lbda,
                      colLbda='LAMBDA', colTrans='OZONE', ext=1):
    """
    Read ozone transmission template, interpolate over
    wavelengthes, and convert to extinction [mag/airmass].

    :param ozoneName: input FITS table, with columns *colLbda*
      (wavelength in AA) and *colTrans* (fractional transmission), and
      key 'REFO3COL' specifing the reference ozone column density [DU]
    :param lbda: output wavelengthes [AA]
    :param colLbda: name of the wavelength (in AA) column
    :param colTrans: name of the ozone transmission column
    :param ext: extension in which to look for wavelength and
      transmission columns
    :return: ozone extinction [mag/airmass], refO3col
    """

    # Read wavelength and transmission columns
    ffile = F.open(ozoneName)
    x = ffile[ext].data.field(colLbda)   # Wavelength
    y = ffile[ext].data.field(colTrans)  # Transmission
    refO3col = ffile[ext].header["REFO3COL"]

    # Interpolate transmissions over lbda
    from scipy.interpolate import UnivariateSpline
    trans = UnivariateSpline(x, y, s=0)(lbda)

    # Convert to extinction [mag/airmass]
    return N.absolute(-2.5 * N.log10(trans)), refO3col


def createTable(arrays, names,
                units=None, formats=None, keywords=(), extname=None):
    """
    Create a FITS-table from a set of arrays.

    :param arrays: list of input arrays (ncols,)
    :param names: list of column names (ncols,)
    :param units: list of column units (ncols,) ('none' by default)
    :param formats: list of column formats (ncols,) ('1E' by default)
    :param keywords: list of keys (key,val[,comment]) to add to table header
    :param extname: name of the binary table extension
    :return: table HDU
    """

    assert len(arrays) == len(names)
    for arr in arrays[1:]:
        assert len(arr) == len(arrays[0])
    if units is None:
        units = ['none'] * len(arrays)
    else:
        assert len(units) == len(arrays)
    if formats is None:
        formats = ['1E'] * len(arrays)
    else:
        assert len(formats) == len(arrays)

    cols = [ F.Column(array=a, name=n, unit=u, format=f)
             for a, n, u, f in zip(arrays, names, units, formats) ]

    thdu = F.new_table(cols)
    if extname is not None:
        thdu.header['EXTNAME'] = extname
    thdu.header['FCLASS'] = (26, 'Table file class')
    for key, val, cmt in keywords:        # Add some keywords if any
        thdu.header[key] = (val, cmt)

    return thdu

# End of AtmosphericExtinction.py ==============================================
