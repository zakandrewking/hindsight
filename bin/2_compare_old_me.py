#!/usr/bin/env python -m pdb
# -*- coding: utf-8 -*-

from me.problem_and_solver_classes.me_solution import LP_Solution
from me_scripts.db.queries import m_reactions_for_me_reactions
from me_scripts.hindsight.pathways import get_designs, add_all_heterologous_pathways
from me_scripts.hindsight.hindsight import me_exchange_for_metabolite_name

import pandas as pd
idx = pd.IndexSlice
import numpy as np
import re
import json
import os
import cPickle as pickle
import itertools as it
from theseus import carbons_for_exchange_reaction, load_model

# functions

def get_m_secretion(s):
    def get_m_name(name):
        return 'EX_%s_e' % name.replace('_DASH_', '__') if 'EX_' not in name else name
    s = [(get_m_name(k), v) for k,v in s.iteritems()]
    return sorted(s, key=lambda x: x[1], reverse=True)

def get_batch_target_secretion(ser):
    me_rxn = me_exchange_for_metabolite_name(ser['target_exchange'], source=False)
    if type(ser['io']) is float and np.isnan(ser['io']):
        return np.nan
    try:
        return ser['io'][me_rxn]
    except KeyError:
        pass
    return 0

def get_batch_target_yield(ser, t='min'):
    if type(ser['io']) is float and np.isnan(ser['io']):
        return np.nan

    # get carbon uptake
    substrates = ser['substrate_exchange']
    if type(substrates) is not list:
        substrates = [substrates]
    total_c_uptake = 0
    for m_exchange in substrates:
        me_rxn = me_exchange_for_metabolite_name(m_exchange, source=True)
        me_sink_rxn = me_exchange_for_metabolite_name(m_exchange, source=False)

        carbons = carbons_for_exchange_reaction(hete_model.reactions.get_by_id(m_exchange))
        if me_rxn in ser['io']:
            uptake_flux = ser['io'][me_rxn]
            if me_sink_rxn in ser['io']:
                print '{}, {}'.format(ser['citation_key'], ser[u'io'][me_sink_rxn])
                uptake_flux -= ser['io'][me_sink_rxn]
            total_c_uptake += carbons * uptake_flux
        else:
            print '{}: could not find uptake for {}. gr: {}'.format(ser[u'citation_key'], me_rxn, ser['gr'])
            pass

    # get carbon secretion
    out_carbons = carbons_for_exchange_reaction(hete_model.reactions.get_by_id(ser['target_exchange']))
    out_c_flux = out_carbons * ser[t]
    return out_c_flux / total_c_uptake

def nan_or_zero(x):
    return bool(np.isnan(x) or (x < 1e-4))

# Old ME results

loaded = {}
for directory in ['/Users/zaking/data/hindsight/15..03.18_native_sims/',
                  '/Users/zaking/data/hindsight/15..03.18_plasmid_sims/']:
    for path in os.listdir(directory):
        if path.endswith('.json'): continue
        for path2 in os.listdir(os.path.join(directory, path)):
            loaded[path2] =  LP_Solution.import_from_json(os.path.join(directory, path, path2))

def add_plasmid_flux(m_flux, me_flux, plasmid_reactions, only_exchange=False):
    """Add me_flux reactions to m_flux dictionary if they are in the plasmid_reactions list"""
    for reaction in plasmid_reactions:
        if only_exchange and not any([s in reaction for s in ['EX_', '_sink', '_source']]):
            continue
        try:
            m_flux[reaction] = me_flux[reaction]
        except KeyError:
            pass
    return m_flux

# get all the reactions that are in a plasmid design
all_plasmid_reactions = set(it.chain.from_iterable([d[1].keys() for d in
                                                    get_designs().itervalues()]))

i = []; o = []
for k, v in loaded.iteritems():
    if v.status != 'OPTIMAL':
        print 'Bad status %s in %s' % (v.status, k)
        continue
    citation_key = re.finditer(r'^([a-zA-Z0-9]+(?:_[0-9])?(?:_ethanol)?)', k).next().group(0)
    # fix the Wang2012a issue
    if citation_key == 'Wang2012a_1':
        continue
    elif citation_key == 'Wang2012a_2':
        citation_key = 'Wang2012a'
    case = 'ME_atoB_irrev' if 'atoB_irrev' in k else 'ME'
    sim = re.finditer(r'((?:batch)|(?:fixed_growth_Minimize)|(?:fixed_growth_Maximize))', k).next().group(0)
    i.append((citation_key, case, sim))
    o.append({'objective': v.lp_problem_definition.basic_model_parameters.growth_rate_in_per_hour if sim=='batch' else v.objective_value,
              'me_io': {k: v for k, v in v.primal_dictionary.iteritems() if ('_sink' in k or '_source' in k or 'EX_' in k)},
              'm_flux': add_plasmid_flux(v.get_metabolic_flux(), v.primal_dictionary, all_plasmid_reactions)})
index = pd.MultiIndex.from_tuples(i, names=['paper', 'model', 'sim'])
me_in = pd.DataFrame(o, index=index)

# unstack
me_results = me_in.unstack('sim')
# rearrange cols
new_col = me_results.columns.set_levels([[u'objective', u'me_io', u'm_flux'], [u'batch', u'fixed_growth_Minimize', u'fixed_growth_Maximize']])
me_results = me_results.reindex_axis(new_col, axis=1)
me_results['gr'] = me_results['objective']['batch']
me_results['min'] = me_results['objective']['fixed_growth_Minimize']
me_results['max'] = me_results['objective']['fixed_growth_Maximize']
#me_results['max_flux'] = me_results['flux']['fixed_growth_Maximize']
me_results['flux'] = me_results['m_flux']['batch']
me_results['io'] = me_results['me_io']['batch']

# remove multiindex
me_results_print = me_results[['gr', 'min', 'max', 'flux', 'io']]
me_results_print.columns = pd.Index([x[0] for x in me_results_print.columns])

## Combine with M sims

with open('m_sims_table.pickle', 'r') as f:
    m_sims = pickle.load(f)

# combine m with me
m_sims2 = m_sims.reset_index(level=0, drop=True)
all_sims = pd.concat([m_sims2, me_results_print], axis=0)
all_sims = all_sims.sortlevel()

# combine m with me
m_sims2 = m_sims.reset_index(level=0, drop=True)
all_sims3 = pd.concat([m_sims2, me_results_print], axis=0)
all_sims3 = all_sims3.sortlevel()
# reorder
all_sims3 = all_sims3.reindex(pd.MultiIndex.from_product([all_sims3.index.levels[0],
                                                          ['E. coli core', 'iJR904',
                                                           'iAF1260', 'iAF1260b',
                                                           'iAF1260b_gene_ko', 'iJO1366',
                                                           'iJO1366_gene_ko',
                                                           'iJO1366_atoB_irrev', 'ME',
                                                           'ME_atoB_irrev']],
                                                         names=['Paper', 'Model']))

# propogate info
cols_to_propagate = ['year', 'target_exchange', 'Deletions', 'Aerobicity',
                     'Additions', 'Target', 'c_byproduct_order', 'Evolved', 'citation_key',
                     'substrate_exchange', 'PMID', 'authors', 'title', 'strategies', 'In silico prediction',
                     'Native?', 'Parent strain', 'Substrate']
all_sims3[cols_to_propagate] = all_sims3[cols_to_propagate].groupby(level='Paper').fillna(method='backfill')
all_sims3[cols_to_propagate] = all_sims3[cols_to_propagate].groupby(level='Paper').fillna(method='pad')
# add year to index
all_sims3 = all_sims3.set_index('year', append=True)
all_sims3 = all_sims3.sort_index()

### Temporary solution to Single_Exchange_FVA bug

hete_model = add_all_heterologous_pathways(load_model('iJO1366'))

# year being NaN causes trouble for idx'ing
all_sims = all_sims3.reset_index(level='year')
all_sims.loc[idx[:, 'ME'], 'min'] = all_sims.loc[idx[:, 'ME'], :].apply(get_batch_target_secretion, axis=1)
all_sims.loc[idx[:, 'ME'], 'max'] = all_sims.loc[idx[:, 'ME'], :].apply(get_batch_target_secretion, axis=1)
all_sims.loc[idx[:, 'ME'], 'yield_min'] = all_sims.loc[idx[:, 'ME'], :].apply(get_batch_target_yield, axis=1, t='min')
all_sims.loc[idx[:, 'ME'], 'yield_max'] = all_sims.loc[idx[:, 'ME'], :].apply(get_batch_target_yield, axis=1, t='max')
all_sims = all_sims.set_index('year', append=True)
all_sims = all_sims.sort_index()

## Check for sims where iJO grows and ME dies (and vice versa)

ijo_grows_vs_me = (all_sims
 .sort_index()
 .loc[idx[:, ['iJO1366', 'ME'], :], :]
 .groupby(level='Paper')
 .filter(lambda x: (not nan_or_zero(x.xs('iJO1366', level='Model').gr.values[0])) and
                   (nan_or_zero(x.xs('ME', level='Model').gr.values[0]))))

ijo_grows_vs_me.loc[:, ['gr', 'min', 'max', 'Additions', 'Deletions', 'Deletions_b']]
ijo_grows_vs_me.loc[:, 'reaction_knockouts'].values


### TODO Check these plasmids, and check no-knockout simulations

me_grows_vs_ijo = (all_sims
 .sort_index()
 .loc[idx[:, ['iJO1366', 'ME'], :], :]
 .groupby(level='Paper')
 .filter(lambda x: (nan_or_zero(x.xs('iJO1366', level='Model').gr.values[0])) and
                   (not nan_or_zero(x.xs('ME', level='Model').gr.values[0]))))

me_grows_vs_ijo.loc[:, ['gr', 'min', 'max', 'Deletions', 'reaction_knockouts', 'Additions']]

## Save sims
# all_sims.to_pickle('old_me_m_sims.pickle')
