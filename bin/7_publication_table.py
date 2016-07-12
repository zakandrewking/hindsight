#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight import models_to_compare
from hindsight.pathways import get_designs

import pandas as pd
idx = pd.IndexSlice

all_sims = pd.io.pickle.read_pickle('../data/sims_table.pickle')
all_sims_no_wildtype = all_sims.loc[all_sims.index.map(lambda x: 'wildtype' not in x[0])]

columns_to_export = [
    'pmid',
    'additions',
    'aerobicity',
    'deletions',
    'evolved',
    'simulation',
    'parent',
    'substrate',
    'target',
    'strategies',
]

(all_sims_no_wildtype
 .xs('iJO1366', level='model')
 .loc[:, columns_to_export]
 .to_csv('../data/publication_table.tsv', sep='\t'))
