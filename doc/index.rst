.. IAST documentation master file, created by
   sphinx-quickstart on Thu May 28 12:45:03 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for IAST Package
================================

This Python package takes pure component gas adsorption isotherms in a nanoporous material and predicts a mixture isotherm, how much of each gas will adsorb when the material is in equilibrium with a multi-component mixture. Ideal Adsorbed Solution Theory (IAST) is the framework used to predict the mixture adsorption isotherm from the pure component adsorption isotherms.

This code has three options to apply IAST to the pure component adsorption isotherms:

1. Fit a Langmuir isotherm model
2. Fit a quadratic isotherm model
3. Use linear interpolation (numerical quadrature for spreading pressures)

============
Installation
============

===
Use
===

As an example for use, see the `/test` directory. We test the IAST code with a binary mixture of Xe and Kr in IRMOF-1.

Simulated pure component adsorption isotherms at 298 K are present in:
* IRMOF-1_clean_Xe_isotherm_298K.csv
* IRMOF-1_clean_Kr_isotherm_298K.csv

We ran dual component GCMC mixture isotherms of Xe/Kr in IRMOF-1 at 1 bar total pressure and different mixture conditions. This data is present in mixture_Xe_Kr_IRMOF-1_clean_298K.csv. We use IAST to reproduce this result.

-----------------------------------
Import the pure component isotherms
-----------------------------------

Use the `Pandas` package to load the pure component adsorption isotherms as DataFrames:

.. code-block:: python
    
    import pandas as pd
    df_Xe = pd.read_csv("IRMOF-1_clean_Xe_isotherm_298K.csv", skiprows=1)
    df_Kr = pd.read_csv("IRMOF-1_clean_Kr_isotherm_298K.csv", skiprows=1)

The units for pressure and loading in both DataFrames must be consistent; loading of gas must be in a molar quantity for IAST to apply (e.g. mmol/g or mmol/cm3). The `IAST` package will then work with these units throughout. 

--------------------------
Construct isotherm objects
--------------------------

First, import the `IAST` package:

.. code-block:: python

   import IAST

We separate the process of characterizing the pure component adsorption isotherms from performing IAST. Construct the models by passing the DataFrame with the pure component adsorption isotherm and the names of the columns that correspond to the loading and pressure.

* Langmuir isotherm model

.. code-block:: python

    xe_isotherm = IAST.LangmuirIsotherm(df_Xe, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)")
    kr_isotherm = IAST.LangmuirIsotherm(df_Kr, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)")

* Quadratic isotherm model

.. code-block:: python

    xe_isotherm = IAST.QuadraticIsotherm(df_Xe, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)")
    kr_isotherm = IAST.QuadraticIsotherm(df_Kr, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)")

* Linear interpolation model

The linear interpolation model has an additional, optional argument `fill_value` that tells us what loading to assume when we attempt to extrapolate beyond the data point with the highest pressure. In this example, we assume that the loading at the highest pressure point is equal to the saturation loading. By default, `fill_value=None` and an error is thrown when the IAST code tries to extrapolate loading beyond the highest pressure point in the pure component adsorption isotherm.

.. code-block:: python

    xe_isotherm = IAST.InterpolatorIsotherm(df_Xe, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)", fill_value=df_Xe["Loading(mol/m3)"].max())
    kr_isotherm = IAST.InterpolatorIsotherm(df_Kr, loading_key="Loading(mol/m3)", pressure_key="Pressure(bar)", fill_value=df_Kr["Loading(mol/m3)"].max())

Once an isotherm model is constructed, we can view the fit of the isotherm model to the data by:

.. code-block:: python

   IAST.plot_isotherm(xe_isotherm)
   
The plot can be rendered on the log-scale by passing options `xlogscale=True` and/or `ylogscale=False`.

Using the adsorption isotherm objects, we can calculate the predicted loading at a new pressure, for example 1.0, by:

.. code-block:: python

   xe_isotherm.loading(1.0)
   
or the spreading pressure via:

.. code-block:: python

   xe_isotherm.spreading_pressure(1.0)

-----------------------
Peform IAST calculation
-----------------------

Given the pure component isotherm models, we now illustrate how to use IAST to predict loading when the material is in equilibrium with a mixture of the said gases.

As an example, given `xe_isotherm` and `kr_isotherm` above, we seek the loading [at the same temperature as the pure component isotherms] at a 20/80 mol % Xe/Kr mixture at a pressure of 1.0. To do this, we call:

.. code-block:: python
    
    p = [.2, .8]  # list/array of partial pressures
    q = IAST.IAST(p, [xe_isotherm, kr_isotherm], verboseflag=True)

which will return `q`, an array of component loadings at these mixture conditions predicted by IAST. Entry 0 will correspond to Xe; entry 1 will correspond to Kr.

======
Theory
======

Ideal Adsorbed Solution Theory was developed by Myers and Prausnitz:

Myers, A. L., & Prausnitz, J. M. (1965). Thermodynamics of mixed‐gas adsorption. AIChE Journal, 11(1), 121-127.
    
We follow the method outlined in the reference:

Tarafder, A. and Mazzotti, M. A method for deriving explicit binary isotherms obeying ideal adsorbed solution theory. Chem. Eng. Technol. 2012, 35, No. 1, 102-108.

=====
Tests
=====

This code was tested using pure component Xe and Kr adsorption isotherms in IRMOF-1 to predict the uptake of Xe and Kr at 1 bar in a variety of Xe mole fractions. The test is displayed in [this IPython notebook](http://nbviewer.ipython.org/github/CorySimon/IAST.jl/blob/master/test/test.ipynb), and the files are in the `/test` directory.

TL;DR The following plot shows the simulated Xe/Kr adsorption (points) using binary Grand Canonical Monte Carlo simulations against the IAST prediction from pure-component adsorption isotherms (solid lines).

.. image:: validation.png
   :align: center

==============================
Some notes for troubleshooting
==============================

For using linear interpolation: this is not always a feasible option, depending on how high of a pressure to which your adsorption isotherms extend. It may be required to extrapolate beyond the last pressure point, which is where e.g. the Langmuir model fit may be better. At present, if the IAST calculation needs to extrapolate the isotherm beyond the highest pressure that you provide, it spits out an `NA`. You can change the interpolation strategy from `BCnan` to `BCnearest` if you want to assume the loading beyond the highest pressure provided is equal to the loading at the highest pressure provided. This is a dangerous option, so I did not make this default.

For fitting models (Langmuir, Quadratic) to the pure-component isotherms: the Optim package in Python may fail to fit your adsorption isotherm. A common cause is that your adsorption isotherm does not exhibit enough curvature in order to estimate the saturation loading; that is, you need to obtain data on the uptake at higher pressures. This can be an advantage of the approach using numerical quadrature and linear interpolation: you may not need to know the loading beyond a certain pressure (depends on the conditions in which you are interested).

===============================
Class documentation and details
===============================

.. automodule:: isotherms
   :members:

.. automodule:: IAST
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

