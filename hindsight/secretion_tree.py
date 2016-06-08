"""Calculate the tree of all possible secretions if fermentation products are
sequentially knocked out.

"""

from hindsight import (load_models_to_compare, get_secretion, setup_for_series,
                       apply_design, apply_environment, me_optimize_growth,
                       Design, Environment, SimulationSetup)
from hindsight.variables import min_biomass, SetUpModelError
from parallel_pandas import apply_p
from theseus import load_model

import pandas as pd
idx = pd.IndexSlice
import numpy as np
from collections import namedtuple
import sys
from os.path import join, exists

Secretions = namedtuple('Secretions', [
    'knockouts', # iterable
    'growth_rate', # float
    'exchange_reactions', # iterable
    'fluxes', # iterable with same order as exchange reactions
])
Secretions.__repr__ = lambda self: '<Secretions {self.knockouts}, {self.growth_rate:.3f}, {self.exchange_reactions}, {self.fluxes} >'.format(self=self)

class FoundReaction(Exception):
    pass

def run_secretions_for_knockouts_dataframe(df, outdir, threads,
                                           with_gene_kos=False, debug=False):
    """Runs secretion tree for the whole DataFrame.

    Arguments
    ---------

    df: A DataFrame generated by 1_run_sims.py.

    outdir: A directory to save results.

    threads: The number of threads.

    with_gene_kos: If True, first apply the gene knockouts from the design.

    debug: If True, then just run a single simulation.

    """
    loaded_models = load_models_to_compare(m_only=True)
    if debug:
        res = []
        df_run = df.loc[idx[: , ['ME']], :]
        for index, row in df_run.reset_index().iterrows():
            res.append(run_secretions_for_knockouts_series(row,
                                                           loaded_models=loaded_models,
                                                           with_gene_kos=with_gene_kos,
                                                           outdir=outdir))
        out_df = pd.DataFrame(res).set_index(df.index.names)
    else:
        out_df = apply_p(df.loc[idx[:, loaded_models.keys()], :],
                         run_secretions_for_knockouts_series,
                         loaded_models=loaded_models, outdir=outdir,
                         with_gene_kos=with_gene_kos,
                         threads=threads)
    out_df.to_pickle(join(outdir, 'secretion_tree_results.pickle'))

def run_secretions_for_knockouts_series(series, loaded_models=None,
                                        with_gene_kos=None, outdir=None):
    outfile = join(outdir, '%s_%s.json' % (series['paper'], series['model']))
    if exists(outfile):
        print('Already exists: %s' % outfile)
        return pd.read_json(outfile, typ='series')
    print('Running: %s' % outfile)

    # set up
    setup = setup_for_series(series, loaded_models, True)

    if not with_gene_kos:
        # do not knock out the genes in the model
        setup = SimulationSetup(
            setup.model,
            setup.environment,
            Design(setup.design.heterologous_pathway, [], setup.design.target_exchange),
            setup.use_greedy_knockouts
        )

    ignore_exchanges = ['EX_h2_e', 'EX_o2_e', 'EX_co2_e']
    try:
        found, data = secretions_for_knockouts(setup,
                                                  return_if_found=setup.design.target_exchange,
                                                  ignore_exchanges=ignore_exchanges)
        out_series = pd.Series({
            'target_exchange': setup.design.target_exchange,
            'can_secrete': found,
            'data': data,
        })
    except SetUpModelError as e:
        return pd.Series({
            'paper': series['paper'],
            'model': series['model'],
            'error': e,
        })

    out_series['paper'] = series['paper']
    out_series['model'] = series['model']

    out_series.to_json(outfile)
    return out_series

def secretions_for_knockouts(setup, knockouts=[], max_depth=10, depth=0,
                             ignore_exchanges=[], return_if_found=None,
                             growth_cutoff=min_biomass, flux_cutoff=0.1):
    """Accepts a SimulationSetup and a set of knockouts.

    Returns a tree of secretions using nested dictionaries.

    Arguments
    ---------

    setup: SimulationSetup.

    knockouts: A list of reaction IDs to knock out.

    max_depth: The maximum depth to search.

    depth: The current depth.

    ignore_exchanges: Exchanges to not knock out.

    return_if_found: A reaction ID that, if found, will raise FoundReaction exception.

    growth_cutoff: Below this growth rate, the simulation is considered lethal.

    flux_cutoff: The minimum flux required to raise for return_if_found.

    """
    # check depth
    if depth > max_depth:
        return False, None

    # always copy the model
    model = setup.model
    if model.id == 'ME':
        model = load_model(model.id)
    else:
        model = model.copy()

    # copy environment for changes, knock out the reactions by adding them to
    # other_bounds
    environment = Environment(setup.environment.substrate_exchanges,
                              setup.environment.supplement_exchanges,
                              setup.environment.aerobic,
                              dict({ko: (0, 0) for ko in knockouts},
                                   **setup.environment.other_bounds))

    # set up model. Have to do this every time because the ME model cannot be
    # copied
    model = apply_environment(model, environment)
    model = apply_design(model, setup.design, setup.use_greedy_knockouts)

    # solve the problem
    sol = me_optimize_growth(model) if model.id == 'ME' else model.optimize()

    if sol.f is None or sol.f <= growth_cutoff:
        return False, None
    else:
        secretion = dict(get_secretion(model, sol.x_dict, sort=False))
        if return_if_found is not None and return_if_found in secretion and secretion[return_if_found] > flux_cutoff:
            can_secrete = True
            children = None
        else:
            children_raw = {new_knockout: secretions_for_knockouts(setup, knockouts + [new_knockout],
                                                                   max_depth, depth + 1,
                                                                   ignore_exchanges, return_if_found,
                                                                   growth_cutoff, flux_cutoff)
                            for new_knockout in secretion.iterkeys()
                            if new_knockout not in ignore_exchanges}
            can_secrete = any(v[0] for v in children_raw.itervalues())
            children = {k: v[1] for k, v in children_raw.iteritems()}
        return can_secrete, {
            'knockouts': knockouts,
            'growth_rate': sol.f,
            'secretion': secretion,
            'children': children,
        }
