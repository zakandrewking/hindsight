#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight.secretion_tree import run_secretions_for_knockouts_dataframe

import pandas as pd
idx = pd.IndexSlice
from os import makedirs

df = pd.read_pickle('../data/sims_table.pickle')
tree_dir = '../data/secretion_tree'
tree_w_kos_dir = '../data/secretion_tree_w_gene_kos'
for d in tree_dir, tree_w_kos_dir:
    try:
        makedirs(d)
    except OSError:
        pass
run_secretions_for_knockouts_dataframe(df, tree_dir, 18, with_gene_kos=False)
run_secretions_for_knockouts_dataframe(df, tree_w_kos_dir, 18, with_gene_kos=True)
