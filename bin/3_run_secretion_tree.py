#!/usr/bin/env python -W ignore -u
# -*- coding: utf-8 -*-

from me_scripts.hindsight.secretion_tree import run_secretions_for_knockouts_dataframe
import pandas as pd

df = pd.read_pickle('sims_table.pickle')
out = run_secretions_for_knockouts_dataframe(df, 'results.json', 10, debug=False)
