# -*- coding: utf-8 -*-

from theseus.bigg import download_model
from theseus import load_model
from me_scripts.hindsight.hindsight import *

import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal
import pytest

ex_environment = Environment(['EX_glc__D_e'], [], True, {})
ex_design = Design(None, [], 'EX_ac_e')

@pytest.fixture(scope='function')
def me_model():
    return load_model('ME')

def test_minimize_maximize():
    model = download_model('iJO1366')
    sim_setup = SimulationSetup(model, ex_environment, ex_design, True)
    out = minimize_maximize(sim_setup)
    assert len(out) == 5
    assert type(out[3]) == list
    assert type(out[3][0]) == tuple

def test_minimize_maximize_me(me_model):
    sim_setup = SimulationSetup(me_model, ex_environment, ex_design, True)
    me_model = apply_environment(me_model, ex_environment)
    me_model = apply_design(me_model, ex_design, True)
    out = minimize_maximize(sim_setup)
    assert out.growth_rate > 0.1
    assert out.minmax[0] > 0.1
    assert out.minmax[1] > 0.1
    assert 'EX_ac_e' in [x[0] for x in out.max_secretion]

def test_get_absolute_max():
    sim_setup = SimulationSetup(me_model, ex_environment, ex_design, True)
    me_model = apply_environment(me_model, ex_environment)
    me_model = apply_design(me_model, ex_design, True)
    abs_max, heterologous = get_absolute_max(sim_setup)
    assert abs_max > 5
    assert np.isnan(heterologous)

def test_check_setup(me_model):
    new_design = Design(None, ['b4154', 'b0903', 'b0001'], 'EX_ac_e')
    sim_setup = SimulationSetup(me_model, ex_environment, new_design, True)
    out = check_setup(sim_setup)
    assert out == ['b0001']

def test_get_reaction_knockouts(me_model):
    new_design = Design(None, ['b1779'], 'EX_ac_e')
    kos, greedy = get_reaction_knockouts(me_model, new_design, True)
    assert 'GAPD_FWD_GAPDH-A-CPLX' in kos
    assert 'E4PD_FWD_GAPDH-A-CPLX' in kos
    assert 'E4PD_FWD_GAPDH-A-CPLX' in greedy
