# -*- coding: utf-8 -*-

from me_scripts.hindsight.pathways import (add_heterologous_pathway,
                                           get_designs,
                                           exchange_for_metabolite_name)
from me_scripts.hindsight.variables import (min_biomass, default_sur, max_our,
                                            NotFoundError, SetUpModelError)

import numpy as np
import pandas as pd
import itertools as it
from collections import namedtuple
from ipy_progressbar import ProgressBar
import cobra
from cobra.manipulation.delete import find_gene_knockout_reactions
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from minime.solve.algorithms import binary_search, solve_at_growth_rate
from theseus import carbons_for_exchange_reaction, load_model, setup_model
from theseus.bigg import download_model

# --------------------------------------------------
# Models
# --------------------------------------------------

models_to_compare = ['e_coli_core', 'iJR904', 'iAF1260', 'iAF1260b',
                     'iJO1366', 'iML1515', 'ME']
m_models_to_compare = ['e_coli_core', 'iJR904', 'iAF1260', 'iAF1260b',
                       'iJO1366', 'iML1515']
private_models = ['iML1515', 'ME']

def download_or_load_model(name):
    if name in private_models:
        return load_model(name)
    else:
        return download_model(name, host='http://zak.ucsd.edu:8888/api/v2/')

def load_models_to_compare(m_only=False):
    if m_only:
        return {n: download_or_load_model(n) for n in m_models_to_compare}
    else:
        return {n: download_or_load_model(n) for n in models_to_compare}

def find_me_reactions(me_model, m_reaction_id):
    """Returns a list of ME reaction IDs."""
    return [r.id for r in me_model.reactions if r.id.startswith(m_reaction_id + '_')]

def me_optimize_growth(me):
    return binary_search(me, min_mu=0, max_mu=1.1, mu_accuracy=1e-5,
                         compiled_expressions=me.expressions, debug=False)

def me_optimize_target(me, growth_rate):
    return solve_at_growth_rate(me, growth_rate,
                                compiled_expressions=me.expressions,
                                debug=False)

# --------------------------------------------------
# Simulation
# --------------------------------------------------

Environment = namedtuple('Environment', [
    'substrate_exchanges',  # list
    'supplement_exchanges', # list
    'aerobic',              # Boolean
    'other_bounds',         # dict, values are (lb, ub)
])

Design = namedtuple('Design', [
    'heterologous_pathway', # str or None
    'gene_knockouts',       # list, by b number
    'target_exchange',      # str
])

SimulationSetup = namedtuple('SimulationSetup', [
    'model',       # COBRA model, with heterologous pathways
    'environment', # Environment
    'design',      # Design
    'use_greedy_knockouts', # Boolean, for genes
])

ModelComparison = namedtuple('ModelComparison', [
    'models',
    'simulation_setup',
])

MinMaxSolution = namedtuple('MinMaxSolution', [
    'growth_rate',    # float
    'minmax',
    'yield_minmax',
    'max_secretion',
    'max_flux',
])

def apply_environment(model, environment):
    """Modify the model to match the environment."""
    # if there are multiple substrates, divide evenly
    sur = 1000 if model.id == 'ME' else default_sur / len(environment.substrate_exchanges)
    this_max_our = 1000 if model.id == 'ME' else max_our

    for substrate_exchange in environment.substrate_exchanges:
        if not substrate_exchange in model.reactions:
            raise SetUpModelError('Substrate exch. not found: %s' % substrate_exchange)

    model = setup_model(model, environment.substrate_exchanges,
                        aerobic=environment.aerobic, sur=sur,
                        max_our=this_max_our)

    other_bounds = environment.other_bounds
    # deal with me reactions
    if model.id == 'ME':
        new_bounds = {}
        for r_id, bs in other_bounds.items():
            for me_r_id in find_me_reactions(model, r_id):
                if model.reactions.get_by_id(me_r_id).reverse:
                    new_bounds[me_r_id] = (0, -bs[0])
                else:
                    new_bounds[me_r_id] = (0, bs[1])
        other_bounds = new_bounds

    for r_id, bs in environment.other_bounds.items():
        try:
            reaction = model.reactions.get_by_id(r_id)
        except KeyError:
            SetUpModelError('Other bound reaction not found: %s' % r_id)
        reaction.lower_bound, reaction.upper_bound = bs
    return model

def apply_design(model, design, use_greedy_knockouts):
    """Modify the model to match the design."""
    # add non-native pathway
    if design.heterologous_pathway is not None:
        try:
            model = add_heterologous_pathway(model, design.heterologous_pathway)
        except NotFoundError:
            raise SetUpModelError('bad addition: %s' % design.heterologous_pathway)
    # add fhl for hydrogen production
    if design.target_exchange == 'H2':
        model = add_heterologous_pathway(model, 'fhl')
    reaction_kos, _ = get_reaction_knockouts(model, design, use_greedy_knockouts)
    # get me reactions
    if model.id == 'ME':
        reaction_kos = it.chain(*[find_me_reactions(model, r) for r in reaction_kos])
    # perform knockouts
    for r in reaction_kos:
        model.reactions.get_by_id(r).knock_out()
    # check target
    if not design.target_exchange in model.reactions:
        raise SetUpModelError('Target exchange not found: %s' % design.target_exchange)
    return model

def get_absolute_max(sim_setup, copy=False):
    """Find the absolute max production (theoretical yield), and check for usage of
    the nonnative pathways.

    """
    model = sim_setup.model.copy() if copy else sim_setup.model
    biomass_reaction = next(iter(model.objective.keys()))

    # set min biomass
    biomass_reaction.lower_bound = min_biomass

    # objective is target production
    model.change_objective(sim_setup.design.target_exchange)

    sol = me_optimize_target(me, min_biomass) if model.id == 'ME' else model.optimize()

    if sol.f is None:
        return np.nan, np.nan

    # for ME model, sol.f is still growth rate
    absolute_max_product = sol.x_dict[sim_setup.design.target_exchange]

    # get flux through non-native pathways
    designs = get_designs()
    pathway = sim_setup.design.heterologous_pathway
    if pathway in designs:
        reactions = designs[pathway][1]
        flux = model.get_metabolic_flux() if model.id == 'ME' else sol.x_dict
        additions_pathway_flux = {r: flux[r] for r in reactions if not r.startswith('EX_')}
    else:
        additions_pathway_flux = np.nan

    return absolute_max_product, additions_pathway_flux

def minimize_maximize(sim_setup, biomass_fraction=0.9999,
                      calculate_yield=True, copy=False):
    """Run FVA at max growth rate. Uses pFBA for max fluxes.

    """
    model = sim_setup.model.copy() if copy else sim_setup.model
    biomass_reaction = next(iter(model.objective.keys()))

    substrate_exchanges = sim_setup.environment.substrate_exchanges
    target_exchange = sim_setup.design.target_exchange

    sol = me_optimize_growth(model) if model.id == 'ME' else model.optimize()

    growth_rate = 0.0 if sol.f is None else sol.f
    if growth_rate < min_biomass:
        yield_minmax = (np.nan, np.nan)
        minmax = (np.nan, np.nan)
        max_secretion = np.nan
        max_flux = np.nan
        heterologous_pathway_flux = np.nan
    else:
        if model.id == 'ME':
            max_target = sol.x_dict[target_exchange]
            minmax = (max_target, max_target)
            max_secretion = get_secretion(model, max_sol.x_dict)
            max_flux = model.get_metabolic_flux()
        else:
            # set the minimum biomass production
            biomass_reaction.lower_bound = biomass_fraction * growth_rate
            # max and min
            model.change_objective(target_exchange)
            max_sol = optimize_minimal_flux(model)
            min_sol = model.optimize(objective_sense='minimize')
            minmax = (min_sol.f, max_sol.f)
            # get secretions
            max_secretion = get_secretion(model, max_sol.x_dict)
            max_flux = max_sol.x_dict
        if calculate_yield:
            # get target carbons
            carbons = carbons_for_exchange_reaction(model.reactions.get_by_id(target_exchange))
            # loop through targets and add carbons
            def c_yield(sub):
                c = carbons_for_exchange_reaction(model.reactions.get_by_id(sub))
                return c * default_sur / len(substrate_exchanges)
            substrate_carbon_uptake = sum(c_yield(sub) for sub in substrate_exchanges)
            yield_minmax = tuple([x * carbons / substrate_carbon_uptake for x in minmax])
        else:
            yield_minmax = (np.nan, np.nan)

    return MinMaxSolution(growth_rate, minmax, yield_minmax, max_secretion,
                          max_flux)

def check_setup(sim_setup):
    """Return genes that are in the design but not the model."""
    design = sim_setup.design
    model = sim_setup.model
    return [g for g in design.gene_knockouts if g not in model.genes]

def get_reaction_knockouts(model, design, use_greedy_knockouts):
    """Determine the reactions that will be removed by gene knockouts, also
    considering greedy knockouts.

    Returns (reaction knockouts, greedy knockouts)

    """
    reaction_knockouts = set()
    gene_reaction_knockouts = set()
    for g in design.gene_knockouts:
        if g not in model.genes:
            continue
        gene_obj = model.genes.get_by_id(g)
        if use_greedy_knockouts:
            reactions = gene_obj.reactions
            for r in reactions:
                r.knock_out()
            reaction_knockouts = reaction_knockouts.union([x.id for x in reactions])
            rg = [x.id for x in find_gene_knockout_reactions(model, [model.genes.get_by_id(g)])]
            gene_reaction_knockouts = gene_reaction_knockouts.union(rg)
        else:
            r = gene_obj.remove_from_model(model)
            reaction_knockouts = reaction_knockouts.union(r)
    greedy_knockouts = reaction_knockouts.difference(gene_reaction_knockouts)
    return list(reaction_knockouts), list(greedy_knockouts)

def error_series(series, err):
    """Error as series"""
    return series.append(pd.Series({'error': str(err)}))

def run_simulation(series, gene_ko=False, loaded_models=None,
                   use_greedy_knockouts=True):
    """Run a simulation on a DataFrame row.

    Example:

        df.apply(run_simulation, axis=1, loaded_models=loaded_models)

    or:

        from me_scripts.parallel_pandas import apply_p
        apply_p(df, run_simulation, loaded_models=loaded_models)

    The DataFrame should have a multi-index with (paper, model). The function
    uses the keys:

        ['additions', 'substrate', 'target', 'aerobicity', 'deletions_b']

    If any of the following columns are already in the DataFrame, it will throw
    an error:

        [ TODO ]

    """

    # copy the model
    model_id = series['model']
    if model_id == 'ME':
        # this is necessary because I can't copy ME models right now
        model = load_model(model_id)
    else:
        model = loaded_models[model_id].copy()

    # get the substrates and supplements
    substrate_exchanges, supplement_exchanges = exchange_for_metabolite_name(series['substrate'])

    # aerobicity
    aerobic = series['aerobicity'].strip().lower() == 'aerobic'

    # heterologous pathway
    additions = series['additions']
    heterologous_pathway = None if additions.strip() == '' else additions

   # knockouts
    gene_knockouts = series['deletions_b']

    # target_exchange
    try:
        target_exchange = exchange_for_metabolite_name(series['target'])[0][0]
    except NotFoundError:
        return error_series(series, 'Bad target name: %s' % target)

    # other bounds
    other_bounds = {}
    if target_exchange == 'EX_h2_e':
        if 'FHL' in model.reactions:
            other_bounds['FHL'] = (0, 1000)

    # environment is generic for any model
    environment = Environment(substrate_exchanges, supplement_exchanges,
                              aerobic, other_bounds)
    design = Design(heterologous_pathway, gene_knockouts, target_exchange)
    setup = SimulationSetup(model, environment, design, use_greedy_knockouts)

    # check the setup
    genes_not_in_model = check_setup(setup)
    reaction_knockouts, greedy_knockouts = get_reaction_knockouts(model, design,
                                                                  use_greedy_knockouts)

    # update the model
    model = apply_environment(model, environment)
    model = apply_design(model, design, use_greedy_knockouts)

    # run minimize_maximize
    try:
        min_max_solution = minimize_maximize(setup)
    except SetUpModelError as e:
        return error_series(series, e)

    # run theoretical yield
    try:
        absolute_max_product, heterologous_pathway_flux = get_absolute_max(setup)
    except SetUpModelError as e:
        return error_series(series, e)

    new_series = pd.Series({
        'substrate_exchanges': environment.substrate_exchanges,
        'supplement_exchanges': environment.supplement_exchanges,
        'aerobic': environment.aerobic,
        'other_bounds': environment.other_bounds,
        'heterologous_pathway': design.heterologous_pathway,
        'target_exchange': design.target_exchange,
        'gene_knockouts': design.gene_knockouts,
        'use_greedy_knockouts': use_greedy_knockouts,
        'target_exchange': design.target_exchange,
        'genes_not_in_model': genes_not_in_model,
        'reaction_knockouts': reaction_knockouts,
        'greedy_knockouts': greedy_knockouts,
        'growth_rate': min_max_solution.growth_rate,
        'min': min_max_solution.minmax[0],
        'max': min_max_solution.minmax[1],
        'yield_min': min_max_solution.yield_minmax[0],
        'yield_max': min_max_solution.yield_minmax[1],
        'max_secretion': min_max_solution.max_secretion,
        'max_flux': min_max_solution.max_flux,
        'absolute_max_product': absolute_max_product,
        'heterologous_pathway_flux': heterologous_pathway_flux,
        'error': np.nan,
    })

    # convert to columns
    return series.append(new_series)

def calculate_envelopes(model_comparison):
    """Calculate production envelopes for wildtype vs design."""
    target = model_comparison.simulation_setup.design.target_exchange
    envs = []
    for model in model_comparison.models:
        envs.append(calculate_production_envelope(model.copy(), target, min_biomass,
                                                  label=model.id))
    return envs

def quick_production_envelope(model, target, min_biomass=0, axis=None,
                              plot_kwargs={}):
    """Plot an envelope for the model."""
    p = calculate_production_envelope(model, target, min_biomass)
    ax = plot_production_envelope(p, axis=axis, plot_kwargs=plot_kwargs)
    if hasattr(model, 'knockouts'):
        ax.set_title('{}: Envelope for {}'.format(' '.join(model.knockouts), target))
    else:
        ax.set_title('Envelope for {}'.format(target))

# --------------------------------------------------
# Analysis
# --------------------------------------------------

def get_secretion(model, flux_dict, flux_threshold=1e-3, sort=True,
                  non_carbon_secretion=['EX_h2_e']):
    """Get all secretions with fluxes greater than flux_threshold.

    Looks for reactions that secrete carbon or are in non_carbon_secretion.

    """
    secretions = [(k, v) for k, v in flux_dict.iteritems()
                  if 'EX_' in k and v > flux_threshold and
                  (carbons_for_exchange_reaction(model.reactions.get_by_id(k)) > 0 or
                   k in non_carbon_secretion)]
    if sort:
        return sorted(secretions, key=lambda x: x[1], reverse=True)
    else:
        return secretions
