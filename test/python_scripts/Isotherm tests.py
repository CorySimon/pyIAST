# coding: utf-8

# # Test isotherm fitting

# Our strategy here is to generate data points that follow a given isotherm model, then fit an isotherm model to the data using pyIAST, and check that pyIAST identifies the parameters correctly.

# In[1]:

from __future__ import absolute_import
from __future__ import print_function
import pyiast
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# We test all analytical models implemented in pyIAST.

# In[2]:

models = pyiast._MODELS
models

# This dictionary gives the model parameters for which we generate synthetic data to test pyIAST fitting. Note that, because the DSLF model has so many parameters, it is highly likely that such a model will overfit the data. Thus, we expect pyIAST to reach a local minimum for DSLF yet still obtain a reasonable fit with the default starting guess.

# In[3]:

model_params = {
    "Langmuir": {
        "M": 10.0,
        "K": 10.0
    },
    "Quadratic": {
        "M": 10.0,
        "Ka": 10.0,
        "Kb": 10.0**2 * 3
    },
    "BET": {
        "M": 10.0,
        "Ka": 10.0,
        "Kb": .2
    },
    "DSLangmuir": {
        "M1": 10.0,
        "K1": 1.0,
        "M2": 30.0,
        "K2": 30.0
    },  # warning: 1/2 is arbitrary
    "Henry": {
        "KH": 10.0
    },
    "TemkinApprox": {
        "M": 10.0,
        "K": 10.0,
        "theta": -0.1
    }
}

# The loading function generates synthetic data for a given model. We pass it an array of pressures and it returns loading using the given model. Note that the parameters for each model are taken from the above dictionary.

# In[4]:


def loading(P, model):
    """
    Return loading at pressure P using a given model.
    
    :param P: np.array array of pressures
    :param model: string specify model
    """
    if model not in models:
        raise Exception("This model is not implemented in the test suite.")

    if model == "Langmuir":
        M = model_params[model]["M"]
        K = model_params[model]["K"]
        return M * K * P / (1.0 + K * P)

    if model == "Quadratic":
        M = model_params[model]["M"]
        Ka = model_params[model]["Ka"]
        Kb = model_params[model]["Kb"]
        return M * P * (Ka + 2.0 * Kb * P) / (1.0 + Ka * P + Kb * P**2)

    if model == "BET":
        M = model_params[model]["M"]
        Ka = model_params[model]["Ka"]
        Kb = model_params[model]["Kb"]
        return M * Ka * P / (1.0 - Kb * P) / (1.0 - Kb * P + Ka * P)

    if model == "DSLangmuir":
        M1 = model_params[model]["M1"]
        K1 = model_params[model]["K1"]

        M2 = model_params[model]["M2"]
        K2 = model_params[model]["K2"]

        return M1 * K1 * P / (1.0 + K1 * P) + M2 * K2 * P / (1.0 + K2 * P)

    if model == "TemkinApprox":
        M = model_params[model]["M"]
        K = model_params[model]["K"]
        theta = model_params[model]["theta"]

        fractional_langmuir_loading = K * P / (1.0 + K * P)
        return M * (
            fractional_langmuir_loading + theta * fractional_langmuir_loading**
            2 * fractional_langmuir_loading)

    if model == "Henry":
        return model_params[model]["KH"] * P


# ## Test model fits

# Loop through all models, generate synthetic data using parameters in `model_params` and the `loading` function here, then fit model using pyIAST. Plot data and fits, check that pyIAST identified parameters match the model.

# In[5]:

for model in models:
    print("Testing model:", model)

    # Generate synthetic data
    df = pd.DataFrame()
    df['P'] = np.linspace(0, 1, 40)
    df['L'] = loading(df['P'], model)

    # use pyIAST to fit model to data
    isotherm = pyiast.ModelIsotherm(
        df, pressure_key='P', loading_key='L', model=model)
    isotherm.print_params()

    # plot fit
    P_plot = np.linspace(0, 1, 100)

    fig = plt.figure()
    plt.scatter(df['P'], df['L'], label='Synthetic data', clip_on=False)
    plt.plot(P_plot, isotherm.loading(P_plot), label='pyIAST fit')
    plt.xlim([0, 1])
    plt.ylim(ymin=0)
    plt.xlabel('Pressure')
    plt.ylabel('Uptake')
    plt.title(model)
    plt.legend(loc='lower right')
    plt.show()

    # assert parameters are equal
    for param in isotherm.params.keys():
        np.testing.assert_almost_equal(
            isotherm.params[param], model_params[model][param], decimal=3)

# ### Quick visual test on the Interpolator isotherm

# In[7]:

isotherm = pyiast.InterpolatorIsotherm(df, pressure_key='P', loading_key='L')
pyiast.plot_isotherm(isotherm)
