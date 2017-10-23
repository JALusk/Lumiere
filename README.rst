SuperBoL: Supernova Bolometric Lightcurves
==========================================

Version 0.3.4

SuperBoL is a python package for calculating the bolometric lightcurves of Type II
supernovae using observed photometry. Three different methods for calculating
the bolometric luminosity are currently included. Those are:

* **Quasi-bolometric**: converts observed magnitudes to monochromatic fluxes at
  the effective wavelengths of the filters, then integrates using the
  trapezoidal rule to get an approximation of the total observed flux.

* **Direct**: Calculates the quasi-bolometric lightcurve, then makes UV and IR
  corrections by fitting a blackbody curve to the fluxes, and integrating that
  function beyond the wavelength limits of the observations.

* **Bolometric Correction**: Performs a bolometric correction based on B-V, V-I,
  or B-I color, using the method of Bersten & Hamuy (2009).

Typical usage often looks like this::

    from superbol import sn

    my_supernova = sn.SN('sn1998a')
    my_supernova.lqbol()                # quasi-bolometric lightcurve
    my_supernova.lbol_direct_bh09()     # direct lightcurve
    my_supernova.lbol_bc_bh09('B', 'V') # B-V bolometric correction lightcurve

SuperBoL propagates uncertainties in the input data through the calculations made
by the code, allowing for errorbars to be included in plots of the lightcurve.

Installation
============

Source code can be found at https://github.com/JALusk/SuperBoL

In order to install SuperBoL system-wide, use::

    python setup.py install

If you do not have root priviliges on your machine, then use::

    python setup.py install --user

Documentation
=============

Documentation is hosted at `ReadTheDocs <http://superbol.readthedocs.io>`_.

Documentation is automatically generated via Sphinx.
To generate the documentation locally, navigate to the ``docs/`` directory and use::

    make html

This will generate html files in the ``docs/build/html/`` directory.
Double-clicking the ``index.html`` file should open the documentation in your
web browser.

Development
===========

Bug reports, feature requests, and contributions are welcome.

Please issue bug reports and feature requests using https://github.com/JALusk/SuperBoL/issues - signing up for a GitHub account is free and easy.

To contribute code, please use the following procedure:

* Fork the SuperBoL repository on GitHub (follow the steps here:
  https://help.github.com/articles/fork-a-repo/)

* Create a branch off of ``develop`` with a descriptive name::

    git checkout -b my_feature develop

* Make your changes and additions

* Write documentation for classes and functions using the Google style guide:
  http://www.sphinx-doc.org/en/stable/ext/example_google.html#example-google

* Write new unit tests for your new code

* Run the test suite using::

    make test

* Issue a pull request (https://help.github.com/articles/using-pull-requests/)
  to the ``develop`` branch on GitHub

* Your pull request will be reviewed and, after any conflicts are resolved,
  merged into the ``develop`` branch and eventually into the next release!

