from theseus import load_model, setup_model
from me_scripts.hindsight.hindsight import (minimize_maximize, normalize_by_carbons,
                                            calculate_category_map, me_exchange_for_metabolite_name,
                                            _category_growth_coupled, _category_lethal,
                                            _category_non_unique)
                                            # , _category_rate_yield

import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal

def test_minimize_maximize():
    model = load_model('iJO1366')
    model = setup_model(model, 'EX_glc_e')
    out = minimize_maximize(model, 'EX_ac_e', 10, 'EX_glc_e')
    assert len(out) == 5
    assert type(out[3]) == list
    assert type(out[3][0]) == tuple
    out = minimize_maximize(model, 'EX_ac_e', 10, 'EX_glc_e', solver='glpk')


def test_me_exchange_for_metabolite_name():
    assert me_exchange_for_metabolite_name('EX_glc_e', True) == 'D-Glucose_source'
    assert me_exchange_for_metabolite_name('EX_glc_e', False) == 'D-Glucose_sink'
    assert me_exchange_for_metabolite_name('EX_1poh_e', False) == 'EX_1poh_e'


def test_normalize_by_carbons():
    model = load_model('iJO1366')
    out = normalize_by_carbons([('EX_glc_e', 1), ('EX_ac_e', 1)], model=model)
    assert out == [('EX_glc_e', 0.75), ('EX_ac_e', 0.25)]

def test__category_growth_coupled():
    df = pd.DataFrame({'yield_min': [0.3, 0.1]})
    assert_array_equal(_category_growth_coupled(df)[1].values, [True, False])

def test__category_lethal():
    df = pd.DataFrame({'gr': [np.nan, 1e-10, 0.4]})
    assert_array_equal(_category_lethal(df)[1].values, [True, True, False])

def test__category_non_unique():
    df = pd.DataFrame({'yield_max': [np.nan, 0.4, 0.25],
                       'yield_min': [np.nan, 0.5, 1e-3]})
    assert_array_equal(_category_non_unique(df)[1].values, [False, False, True])

# def test__category_rate_yield():
#     df = pd.DataFrame({'Target': ['Ethanol', 'D-Lactate', 'L-Alanine', 'Ethanol']},
#                       index=pd.MultiIndex(levels=[['p1', 'p2'], ['ME', 'iJO1366'], [1999, 2000]],
#                                           labels=[[0, 0, 1, 1], [0, 1, 0, 1], [0, 0, 1, 1]],
#                                           names=[u'Paper', u'Model', u'year']))
#     assert_array_equal(_category_rate_yield(df)[1].values, [True, True, False, False])

def test_calculate_category_map():
    df = pd.DataFrame({'yield_min': [np.nan, np.nan, 0.5, 0.1],
                       'yield_max': [np.nan, np.nan, 0.5, 0.1],
                       'Target': ['E', 'B', 'C', 'D'],
                       'gr': [np.nan, 1e-10, 0.4, 0.5]},
                      index=pd.MultiIndex(levels=[['p1', 'p2'], ['m1', 'm2'], [1999, 2000]],
                                          labels=[[0, 0, 1, 1], [0, 1, 0, 1], [0, 0, 1, 1]],
                                          names=[u'Paper', u'Model', u'year']))
    assert_array_equal(calculate_category_map(df, models=['m1', 'm2']).categories,
                       np.array([['Growth coupled', 'None'],
                                 ['Lethal', 'Lethal']]))
