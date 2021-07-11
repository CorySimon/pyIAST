# coding: utf-8

# # Test

# We test the IAST code with a binary mixture of methane (CH$_4$) and ethane (C$_2$H$_6$) adsorbed in metal-organic framework IRMOF-1.
#
# Simulated pure-component adsorption isotherms at 298 K are present in the files:
# * `IRMOF-1_methane_isotherm_298K.csv`
# * `IRMOF-1_ethane_isotherm_298K.csv`
#
# We ran dual component GCMC mixture isotherms of methane/ethane in IRMOF-1 at 65 bar total pressure and 298 K at different mixture compositions. This data is present in `IRMOF-1_methane_ethane_mixture_isotherm_65bar_298K.csv`.
#
# The goal of this test is to use pyIAST to predict the mixture isotherms from the pure-component isotherms and compare to the binary GCMC mixture simulations.

# In[1]:

from __future__ import absolute_import
from __future__ import print_function
import pyiast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from six.moves import range

# Matplotlib settings
plt.style.use('bmh')

# colors to use for plots
color_key = {'methane': 'g', 'ethane': 'r'}

# ## Load pure-component isotherm data as Pandas dataframes

# ### Methane

# In[2]:

df_ch4 = pd.read_csv("../IRMOF-1_methane_isotherm_298K.csv")
df_ch4.head()

# ### Ethane

# In[3]:

df_ch3ch3 = pd.read_csv("../IRMOF-1_ethane_isotherm_298K.csv")
df_ch3ch3.head()

# ### Plot isotherm data

# In[4]:

fig, ax = plt.subplots()

plt.scatter(
    df_ch3ch3['Pressure(bar)'],
    df_ch3ch3['Loading(mmol/g)'],
    label='Ethane',
    color=color_key['ethane'],
    s=50)
plt.scatter(
    df_ch4['Pressure(bar)'],
    df_ch4['Loading(mmol/g)'],
    label='Methane',
    color=color_key['methane'],
    s=50,
    marker='s')

plt.xlabel('Pressure (bar)')
plt.ylabel('Gas uptake (mmol/g)')

plt.xlim([0, 100])
plt.ylim(ymin=0)

plt.legend(loc='lower right')

plt.tight_layout()
plt.savefig("pure_component_isotherms.pdf", format='pdf')
plt.savefig("pure_component_isotherms.png", format='png', dpi=250)
plt.show()

# ## Import binary-component GCMC mixture isotherms

# This data has the component uptakes in IRMOF-1 at 298 K immersed in a methane/ethane mixture at 65.0 bar total pressure. The column `y_ethane` in the data is the mole fraction of ethane in the mixture.

# In[5]:

df_mixture = pd.read_csv(
    "../IRMOF-1_methane_ethane_mixture_isotherm_65bar_298K.csv")
df_mixture.head()

# ### Plot mixture data

# In[6]:

fig = plt.figure()

plt.scatter(
    df_mixture['y_ethane'],
    df_mixture['EthaneLoading(mmol/g)'],
    color=color_key["ethane"],
    label='Ethane',
    s=50)
plt.scatter(
    df_mixture['y_ethane'],
    df_mixture['MethaneLoading(mmol/g)'],
    color=color_key["methane"],
    label='Methane',
    s=50)

plt.xlabel('Mole fraction ethane in gas phase, $y_{C_2H_6}$')
plt.ylabel('Gas uptake (mmol/g)')

plt.ylim(ymin=0)
plt.xlim([0, .12])

plt.legend(loc='center left')

plt.show()

# ## Use IAST to predict mixture data, compare to dual component GCMC

# Construct isotherm objects. Use the interpolator isotherm here, as Langmuir and Quadratic isotherms do not fit well.
#
# We use fill_value = largest loading so that, when the linear interpolation routine calls a pressure beyond our data, it will yield this value. Essentially, the assumption is that the saturation loading = the highest loading observed in the data.

# ### Interpolator isotherm for Methane

# In[7]:

ch4_isotherm = pyiast.InterpolatorIsotherm(
    df_ch4,
    loading_key="Loading(mmol/g)",
    pressure_key="Pressure(bar)",
    fill_value=df_ch4['Loading(mmol/g)'].max())
pyiast.plot_isotherm(ch4_isotherm)

# ### Interpolator isotherm for ethane

# In[8]:

ch3ch3_isotherm = pyiast.InterpolatorIsotherm(
    df_ch3ch3,
    loading_key="Loading(mmol/g)",
    pressure_key="Pressure(bar)",
    fill_value=df_ch3ch3["Loading(mmol/g)"].max())
pyiast.plot_isotherm(ch3ch3_isotherm)

# ## Perform IAST at same mixture conditions as binary GCMC simulations

# In[9]:

n_mixtures = df_mixture.shape[0]
iast_component_loadings = np.zeros(
    (2, n_mixtures))  # store component loadings here

for i in range(n_mixtures):
    y_ethane = df_mixture['y_ethane'].iloc[i]
    partial_pressures = 65.0 * np.array([y_ethane, 1.0 - y_ethane])
    iast_component_loadings[:, i] = pyiast.iast(
        partial_pressures, [ch3ch3_isotherm, ch4_isotherm], verboseflag=False)

# ## Compare pyIAST predictions to binary GCMC

# In the following plot, the points are the dual component GCMC simulation loadings at the respective ethane mole fraction in the gas phase. The lines are the result of the IAST calculation. The IAST calculations match the binary GCMC simulations very well.

# In[10]:

fig = plt.figure()

plt.scatter(
    df_mixture['y_ethane'],
    df_mixture['EthaneLoading(mmol/g)'],
    color=color_key["ethane"],
    label='Ethane',
    s=50)
plt.scatter(
    df_mixture['y_ethane'],
    df_mixture['MethaneLoading(mmol/g)'],
    color=color_key["methane"],
    label='Methane',
    s=50)

plt.plot(
    df_mixture['y_ethane'],
    iast_component_loadings[0, :],
    color=color_key['ethane'],
    label='Ethane, IAST')
plt.plot(
    df_mixture['y_ethane'],
    iast_component_loadings[1, :],
    color=color_key['methane'],
    label='Methane, IAST')

plt.xlabel('Mole fraction ethane in gas phase, $y_{C_2H_6}$')
plt.ylabel('Gas uptake (mmol/g)')

plt.ylim(ymin=0)
plt.xlim([0, .12])

plt.legend(loc='center left', prop={'size': 12})

plt.tight_layout()
plt.savefig("IAST_validation.pdf", format='pdf')
plt.savefig("IAST_validation.png", format='png', dpi=250)
plt.show()

# ## Another sanity check

# We use IAST for a three-component mixture of 5 bar methane (all the same!). This should yield the loading at 15 bar.

# In[11]:

iast_component_loadings = pyiast.iast(
    [5., 5., 5.], [ch4_isotherm, ch4_isotherm, ch4_isotherm])
print("Sum of loadings: ", np.sum(iast_component_loadings))
print("Loading at 15 bar: ", ch4_isotherm.loading(15.))
np.testing.assert_almost_equal(
    np.sum(iast_component_loadings), ch4_isotherm.loading(15.))

# # Reverse IAST

# In reverse IAST, we specify the mole fraction in the adsorbed phase and calculate the mole fractions in the gas phase that will yield these adsorbed phase mole fractions. We will test this using the binary component GCMC methane/ethane simulations in `df_mixture`. We define the mole fraction of ethane in the adsorbed phase, $x_{ethane}$, below. ($x_{methane} = 1- x_{ethane}$).

# In[12]:

df_mixture['x_ethane'] = df_mixture['EthaneLoading(mmol/g)'] / (
    df_mixture['EthaneLoading(mmol/g)'] + df_mixture['MethaneLoading(mmol/g)'])
df_mixture.head()

# Perform reverse IAST at same adsorbed phase composition as in the binary GCMC simulations.

# In[13]:

n_mixtures = df_mixture.shape[0]
gas_mole_fractions = np.zeros(
    (2, n_mixtures))  # store bulk gas mole fractions here
component_loadings = np.zeros((2, n_mixtures))  # store component loadings here
P_total = 65.0  # bar, same units as in xe_isotherm and kr_isotherm
for i in range(n_mixtures):
    x_ethane = df_mixture['x_ethane'].iloc[i]
    gas_mole_fractions[:, i], component_loadings[:, i] = pyiast.reverse_iast(
        [x_ethane, 1.0 - x_ethane], P_total, [ch3ch3_isotherm, ch4_isotherm])

# #### Compare reverse IAST to binary GCMC simulations.

# Relationship between mole fraction of ethane in the bulk gas phase and the resulting mole fraction of ethane in the adsorbed phase.

# In[14]:

plt.figure()

plt.scatter(
    df_mixture['x_ethane'],
    df_mixture['y_ethane'],
    color='b',
    label='binary GCMC simulations',
    s=50)
plt.plot(
    df_mixture['x_ethane'].values,
    gas_mole_fractions[0, :],
    color='b',
    label='pyIAST')

plt.xlabel('$x_{C_2H_6}$')
plt.ylabel('$y_{C_2H_6}$')

plt.xlim([0, .4])
plt.ylim([0, .4])

plt.legend(loc='upper left')

plt.tight_layout()
plt.savefig("y_vs_x.pdf", format='pdf')
plt.show()

# # Compare fits for ethane

# In[15]:

ch3ch3_interpolated = pyiast.InterpolatorIsotherm(
    df_ch3ch3,
    loading_key="Loading(mmol/g)",
    pressure_key="Pressure(bar)",
    fill_value=df_ch3ch3["Loading(mmol/g)"].max())

ch3ch3_langmuir = pyiast.ModelIsotherm(
    df_ch3ch3,
    loading_key="Loading(mmol/g)",
    pressure_key="Pressure(bar)",
    model="Langmuir")

p_plot = np.logspace(
    np.log(df_ch3ch3['Pressure(bar)'].min()), np.log(80), num=200)

fig = plt.figure()

plt.xlabel("Pressure (bar)")
plt.ylabel("Ethane uptake (mmol/g)")

plt.scatter(
    df_ch3ch3['Pressure(bar)'],
    df_ch3ch3['Loading(mmol/g)'],
    marker='o',
    color='r',
    s=40,
    label='Data',
    clip_on=False)

plt.plot(
    p_plot,
    ch3ch3_interpolated.loading(p_plot),
    color='k',
    linewidth=1,
    label='Linear interpolation',
    linestyle='--')

plt.plot(
    p_plot,
    ch3ch3_langmuir.loading(p_plot),
    color='k',
    linewidth=1,
    label='Langmuir fit')
plt.xscale("log")
plt.xlim(10**-4, 10**2)
plt.ylim(0, 22)

plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig("Ethane_isotherm.pdf", format='pdf')

plt.show()

# In[16]:

pyiast._MODELS
