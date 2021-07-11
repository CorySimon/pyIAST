import pytest
from pyiast import ModelIsotherm
import numpy as np
import pandas as pd


def test_isotherm_constructor():
    # make some mock data for a langmuir isotherm
    p = np.linspace(0, 10)
    true_K = 1e-1
    q = true_K * p / (1 + true_K * p)
    df = pd.DataFrame(
        data=np.stack((p, q), axis=1),
        columns=['p', 'uptake']
    )

    isotherm = ModelIsotherm(df, loading_key='uptake', pressure_key='p', model='Langmuir')
    assert set(isotherm.params.keys()) == {'M', 'K'}
    assert np.allclose(
        np.sort(list(isotherm.params.values())),
        [0.1, 1.0],
        rtol=1e-4
    )

    params = {
        'K': 1e-1,
        'M': 1.0
    }
    isotherm = ModelIsotherm(params=params, model='Langmuir')
    test_loadings = isotherm.loading(p)
    assert np.allclose(q, test_loadings)

    with pytest.raises(Exception):
        ModelIsotherm(df=None, params=None)
