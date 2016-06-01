#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight.hindsight import load_models_to_compare, run_simulation
from hindsight.variables import min_biomass
from parallel_pandas import apply_p

import pandas as pd
idx = pd.IndexSlice
import numpy as np
import cPickle as pickle
import itertools
import json
from os import listdir
from os.path import join, exists

# load the models
# print('Loading models')
loaded_models = load_models_to_compare(m_only=True)
# loaded_models = {'ME': 'placeholder'}

# print('ATPM: %s' % {n: m.reactions.get_by_id('ATPM').lower_bound
#                     for n, m in loaded_models.iteritems()
#                     if hasattr(m, 'reactions')})

# Reindex a new sims DataFrame to store simulations. Returns a copy
# (http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy)
sims_in = pd.io.pickle.read_pickle('../data/literature_table.pickle')

# Add wildtype
wildtype = pd.DataFrame({'target': ['Ethanol'], 'additions': [''], 'aerobicity': ['anaerobic'],
                         'substrate': ['glucose'], 'deletions_b': [[]], 'deletions': ['']},
                     index=['wildtype_ethanol'])
sims_wt = sims_in.append(wildtype)


# Add models
index = pd.MultiIndex.from_tuples(list((v[0], v[1]) for v in itertools.product(sims_wt.index, loaded_models.keys())),
                                  names=['paper', 'model'])
sims = sims_wt.reindex(index=index, level='paper')

# Test a couple
print('Running test')
test = apply_p(sims.loc[['Lee2005-xo_1']], run_simulation, loaded_models=loaded_models, threads=4)
testp = test.loc[:, ['substrate_exchanges', 'supplement_exchanges', 'aerobic',
                     'other_bounds', 'heterologous_pathway', 'target_exchange', 'gene_knockouts',
                     'use_greedy_knockouts', 'target_exchange', 'genes_not_in_model',
                     'reaction_knockouts', 'greedy_knockouts', 'growth_rate', 'min', 'max',
                     'yield_min', 'yield_max', 'max_secretion', 'max_flux', 'absolute_max_product',
                     'heterologous_pathway_flux', 'error']]

import pdb; pdb.set_trace()

# print('Running me test')
# results = []
# index_names = sims.index.names
out_dir = '/Users/zaking/no-backup-data/hindsight/16..05.01_me_sims_debugging/'
# for i, (ind, row) in enumerate(sims.reset_index().iterrows()):
#     path = join(out_dir, '%s.json' % row['paper'])
#     if exists(path):
#         continue
#     print('Running %d of %d' % ((i + 1), len(sims)))
#     out = run_simulation(row, loaded_models=loaded_models)
#     out.to_json(path)
#     results.append(out)

# import pdb; pdb.set_trace()

# Run all Sims

# print('Running all sims')
# sims_ran_m = apply_p(sims, run_simulation, loaded_models=loaded_models, threads=18)

# sims_ran_m.to_pickle(join(out_dir, 'm_temp.pickle'))

sims_ran_m = pd.read_pickle(join(out_dir, 'm_temp.pickle'))

print('Loading ME sims')
results = []
for filename in listdir(out_dir):
    if not filename.endswith('.json'):
        continue
    path = join(out_dir, filename)
    row = pd.read_json(path, typ='series')
    # None to nan
    row = row.fillna(np.nan)
    results.append(row)
df_me = pd.DataFrame(results)
index_names = sims_ran_m.index.names
sims_ran = pd.concat([sims_ran_m.reset_index(), df_me]).set_index(index_names).sort_index()

import pdb; pdb.set_trace()

with open('sims_table.pickle', 'w') as f:
    pickle.dump(sims_ran, f)

# import pdb; pdb.set_trace()

# bad_addition = sims_ran[(sims_ran.growth_rate > min_biomass) & ((sims_ran.absolute_max_product < 0.01) | sims_ran.absolute_max_product.isnull())]
# HTML(bad_addition.loc[:, ['target', 'additions', 'growth_rate', 'absolute_max_product', 'heterologous_pathway_flux']].to_html())

# def sims_with_additions(sims, additions):
#     return sims[sims.additions == additions]

# no_flux_addition = sims_ran[sims_ran.additions.notnull() &
#                             (sims_ran.heterologous_pathway_flux.apply(lambda x: np.nan if x is np.nan else any([abs(fl) < 0.1 for fl in x.itervalues()])))]
# HTML(no_flux_addition[['target', 'additions', 'absolute_max_product', 'heterologous_pathway_flux']].to_html())

# # get all lethal in iJO cases
# lethal_in_ijo = sims_ran.xs('iJO1366', level='model')
# lethal_in_ijo = lethal_in_ijo.loc[lethal_in_ijo.growth_rate < min_biomass]
# HTML(lethal_in_ijo
#      .loc[:, ['target', 'additions', 'deletions', 'growth_rate']]
#      .to_html())
