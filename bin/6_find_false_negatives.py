#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight import find_summary_lethal_reactions

import pandas as pd

print('loading sims')
all_sims = pd.read_pickle('../data/sims_table.pickle')
all_sims_no_wildtype = all_sims.loc[all_sims.index.map(lambda x: 'wildtype' not in x[0])]

print('finding lethal interactions')
res, ne = find_summary_lethal_reactions(all_sims_no_wildtype)
res.to_pickle('../data/false_negatives.pickle')
print('warnings', ne)
