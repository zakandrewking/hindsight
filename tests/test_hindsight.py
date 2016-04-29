# -*- coding: utf-8 -*-

from theseus.bigg import download_model
from theseus import load_model
from me_scripts.hindsight.hindsight import *

import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal

ex_environment = Environment(['EX_glc__D_e'], [], True, {})
ex_design = Design(None, [], 'EX_ac_e')

def test_minimize_maximize():
    model = download_model('iJO1366')
    sim_setup = SimulationSetup(model, ex_environment, ex_design, True)
    out = minimize_maximize(sim_setup)
    assert len(out) == 5
    assert type(out[3]) == list
    assert type(out[3][0]) == tuple

def test_minimize_maximize_me():
    model = load_model('ME')
    sim_setup = SimulationSetup(model, ex_environment, ex_design, True)
    out = minimize_maximize(sim_setup)
    assert len(out) == 5
    assert type(out[3]) == list
    assert type(out[3][0]) == tuple
