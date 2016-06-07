#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight.secretion_tree import run_secretions_for_knockouts_dataframe
import pandas as pd

df = pd.read_pickle('../data/sims_table.pickle')
# run_secretions_for_knockouts_dataframe(df, '../data/secretion_tree', 12, debug=False)
directory = '../data/secretion_tree_w_gene_kos'
run_secretions_for_knockouts_dataframe(df, directory, 12, with_gene_kos=True)
