#!/usr/bin/env python -W ignore -u
# -*- coding: utf-8 -*-

from hindsight.secretion_tree import run_secretions_for_knockouts_dataframe
import pandas as pd

df = pd.read_pickle('../data/sims_table.pickle')
out = run_secretions_for_knockouts_dataframe(df, '../data/secretion_tree', 10, debug=False)
