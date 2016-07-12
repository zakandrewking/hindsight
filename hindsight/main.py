# -*- coding: utf-8 -*-

from hindsight.pathways import (add_heterologous_pathway, get_designs,
                                exchange_for_metabolite_name)
from hindsight.variables import (min_biomass, default_sur, max_our,
                                 supplement_uptake, NotFoundError,
                                 SetUpModelError)

import numpy as np
import pandas as pd
idx = pd.IndexSlice
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

models_to_compare = [
    'e_coli_core', 'iJR904', 'iAF1260', 'iAF1260b', 'iJO1366',
    # 'iML1515',
    'ME',
]
m_models_to_compare = [
    'e_coli_core', 'iJR904', 'iAF1260', 'iAF1260b', 'iJO1366',
    # 'iML1515',
]
private_models = ['iML1515', 'ME']

def download_or_load_model_me_placeholder(name):
    if name == 'ME':
        return 'placeholder'
    elif name in private_models:
        return load_model(name)
    else:
        return download_model(name, host='http://bigg.ucsd.edu/api/v2/')

IJO1366 = download_or_load_model_me_placeholder('iJO1366')

def load_models_to_compare(m_only=False):
    if m_only:
        return {n: download_or_load_model_me_placeholder(n) for n in m_models_to_compare}
    else:
        return {n: download_or_load_model_me_placeholder(n) for n in models_to_compare}

def find_me_reactions(me_model, m_reaction_id):
    """Returns a list of ME reaction IDs."""
    return [r.id for r in me_model.reactions if r.id.startswith(m_reaction_id + '_')]

def me_optimize_growth(me):
    return binary_search(me, min_mu=min_biomass, max_mu=1.1, mu_accuracy=1e-5,
                         compiled_expressions=me.expressions)

def me_optimize_target(me, growth_rate):
    return solve_at_growth_rate(me, growth_rate,
                                compiled_expressions=me.expressions)

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

    for sub in environment.substrate_exchanges + environment.supplement_exchanges:
        if not sub in model.reactions:
            raise SetUpModelError('Substrate exch. not found: %s' % sub)

    model = setup_model(model, environment.substrate_exchanges,
                        aerobic=environment.aerobic, sur=sur,
                        max_our=this_max_our)

    # explicit other_bounds override supplement_exchanges
    other_bounds = environment.other_bounds
    for supp in environment.supplement_exchanges:
        if supp not in other_bounds:
            other_bounds[supp] = (-supplement_uptake, 1000)

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

    for r_id, bs in environment.other_bounds.iteritems():
        try:
            reaction = model.reactions.get_by_id(r_id)
        except KeyError:
            print('Other bound reaction not found: %s' % r_id)
        else:
            reaction.lower_bound, reaction.upper_bound = bs
    return model

def apply_design(model, design, use_greedy_knockouts,
                 recompile_expressions=True):
    """Modify the model to match the design.

    model:

    design:

    use_greedy_knockouts:

    recompile_expressions: If True, then recompile_expressions when new ME
    reactions are added.

    """
    # add non-native pathway
    if design.heterologous_pathway is not None:
        try:
            model = add_heterologous_pathway(model, design.heterologous_pathway,
                                             recompile_expressions=recompile_expressions)
        except NotFoundError:
            raise SetUpModelError('bad addition: %s' % design.heterologous_pathway)
    # add fhl for hydrogen production
    if design.target_exchange == 'H2':
        model = add_heterologous_pathway(model, 'fhl',
                                         recompile_expressions=recompile_expressions)
    reaction_kos, _ = get_reaction_knockouts(model, design, use_greedy_knockouts)
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

    sol = me_optimize_target(model, min_biomass) if model.id == 'ME' else model.optimize()

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

def carbon_yield(model, target_exchange, substrate_exchanges,
                    supplement_exchanges, x_dict):
    """Calculates the carbon yield"""
    if not x_dict:
        return None
    # get target carbons
    carbons = carbons_for_exchange_reaction(model.reactions.get_by_id(target_exchange))
    # loop through targets and add carbons
    def c_yield(sub):
        c = carbons_for_exchange_reaction(model.reactions.get_by_id(sub))
        return c * -x_dict.get(sub, 0.0)
    substrate_carbon_uptake = sum(c_yield(sub) for sub in substrate_exchanges + supplement_exchanges)
    return x_dict.get(target_exchange, 0.0) * carbons / substrate_carbon_uptake

def minimize_maximize(sim_setup, biomass_fraction=0.9999,
                      calculate_yield=True, copy=False):
    """Run FVA at max growth rate. Uses pFBA for max fluxes.

    """
    model = sim_setup.model.copy() if copy else sim_setup.model
    biomass_reaction = next(iter(model.objective.keys()))

    substrate_exchanges = sim_setup.environment.substrate_exchanges
    supplement_exchanges = sim_setup.environment.supplement_exchanges
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
            max_secretion = get_secretion(model, sol.x_dict)
            max_flux = model.get_metabolic_flux()
            min_flux = model.get_metabolic_flux()
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
            min_flux = min_sol.x_dict
        if calculate_yield:
            yield_minmax = tuple([carbon_yield(model, target_exchange, substrate_exchanges, supplement_exchanges, x_dict)
                                  for x_dict in (max_flux, min_flux)])
        else:
            yield_minmax = (np.nan, np.nan)

    return MinMaxSolution(growth_rate, minmax, yield_minmax, max_secretion,
                          max_flux)

def check_setup(sim_setup):
    """Return genes that are in the design but not the model."""
    design = sim_setup.design
    model = sim_setup.model
    return [g for g in design.gene_knockouts if g not in model.genes and
            (model.id != 'ME' or 'RNA_%s' % g not in model.metabolites)]

def get_reaction_knockouts(model, design, use_greedy_knockouts):
    """Determine the reactions that will be removed by gene knockouts, also
    considering greedy knockouts.

    Returns (reaction knockouts, greedy knockouts)

    """
    is_me_model = model.id == 'ME'
    if is_me_model:
        # just use iJO genes
        me_model = model
        model = IJO1366
    # copy model because we have to change it to determine the reaction
    # knockouts
    model = model.copy()
    reaction_knockouts = set()
    gene_reaction_knockouts = set()
    for g in design.gene_knockouts:
        if g not in model.genes and (not is_me_model or 'RNA_%s' % g not in
                                     model.metabolites):
            continue
        gene_obj = model.genes.get_by_id(g)
        if use_greedy_knockouts:
            reactions = gene_obj.reactions
            reaction_knockouts = reaction_knockouts.union([x.id for x in reactions])
            r = gene_obj.remove_from_model(model)
            if r:
                gene_reaction_knockouts = gene_reaction_knockouts.union(r)
        else:
            r = gene_obj.remove_from_model(model)
            if r:
                reaction_knockouts = reaction_knockouts.union(r)
    greedy_knockouts = reaction_knockouts.difference(gene_reaction_knockouts)

    if is_me_model:
        reaction_knockouts = it.chain(*[find_me_reactions(me_model, r) for r in reaction_knockouts])
        greedy_knockouts = it.chain(*[find_me_reactions(me_model, r) for r in greedy_knockouts])

    return list(reaction_knockouts), list(greedy_knockouts)

def error_series(series, err):
    """Error as series"""
    return series.append(pd.Series({'error': str(err)}))

def setup_for_series(series, loaded_models, use_greedy_knockouts):
    """Get a SimulationSetup for the series."""
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
    return SimulationSetup(model, environment, design, use_greedy_knockouts)

def run_simulation(series, loaded_models=None, use_greedy_knockouts=True):
    """Run a simulation on a DataFrame row.

    Example:

        df.apply(run_simulation, axis=1, loaded_models=loaded_models)

    or:

        from parallel_pandas import apply_p
        apply_p(df, run_simulation, loaded_models=loaded_models)

    The DataFrame should have a multi-index with (paper, model). The function
    uses the keys:

        ['additions', 'substrate', 'target', 'aerobicity', 'deletions_b']

    If any of the following columns are already in the DataFrame, it will throw
    an error:

        [ TODO ]

    """

    setup = setup_for_series(series, loaded_models, use_greedy_knockouts)

    # check the setup
    genes_not_in_model = check_setup(setup)
    reaction_knockouts, greedy_knockouts \
        = get_reaction_knockouts(setup.model, setup.design, setup.use_greedy_knockouts)

    # update the model
    try:
        # set up
        apply_design(setup.model, setup.design, setup.use_greedy_knockouts)
        apply_environment(setup.model, setup.environment)
        # run minimize_maximize
        min_max_solution = minimize_maximize(setup)
        # run theoretical yield
        absolute_max_product, heterologous_pathway_flux = get_absolute_max(setup)
    except SetUpModelError as e:
        return error_series(series, e)

    new_series = pd.Series({
        'substrate_exchanges': setup.environment.substrate_exchanges,
        'supplement_exchanges': setup.environment.supplement_exchanges,
        'aerobic': setup.environment.aerobic,
        'other_bounds': setup.environment.other_bounds,
        'heterologous_pathway': setup.design.heterologous_pathway,
        'target_exchange': setup.design.target_exchange,
        'gene_knockouts': setup.design.gene_knockouts,
        'use_greedy_knockouts': use_greedy_knockouts,
        'target_exchange': setup.design.target_exchange,
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
                  if k.startswith('EX_') and v > flux_threshold and
                  (carbons_for_exchange_reaction(model.reactions.get_by_id(k)) > 0 or
                   k in non_carbon_secretion)]
    if sort:
        return sorted(secretions, key=lambda x: x[1], reverse=True)
    else:
        return secretions

def normalize_by_carbons(l, model=None):
    """Get a list of secretions by total carbon molecules / time.

    Arguments
    ---------

    l: A list of secretions (e.g. generated by get_secretion).

    model: A model containing all possible exchange reactions.

    """
    if type(l) is float and np.isnan(l):
        return l
    not_in_model = [(a,b) for a, b in l if a not in model.reactions]
    # if len(not_in_model) > 0:
    # print 'Ignoring reaction not in %s: %s' % (model, not_in_model)
    l = [(a, b*carbons_for_exchange_reaction(model.reactions.get_by_id(a))) for a, b in l
            if a in model.reactions]
    m = sum([x[1] for x in l])
    return [(x[0], float(x[1])/m) for x in l]

# ---------------------------------------------
# Lethal genotypes
# ---------------------------------------------

LethalInteraction = namedtuple('LethalInteraction', ['model', # str
                                                     'reactions', # list of str IDs
                                                     'growth_rate',
                                                     'example_paper',
                                                     'example_genotype'])

def find_summary_lethal_reactions(sims, max_combinations=5, debug_limit=None, print_solution=False):
    """Returns a summary dataframe of all the lethal combinations and cases where
    max_combinations was not enough to find the lethal interaction.

    """
    out = [] # list of LethalInteraction's
    not_enough_combinations = []
    debug_count = 0
    loaded_models = load_models_to_compare(m_only=True)
    for _, row in sims.reset_index().iterrows():
        if row['model'] not in loaded_models:
            continue
        res, ne = find_lethal_reactions(row, loaded_models, max_combinations,
                                        print_solution=print_solution)
        if res is None:
            continue
        not_enough_combinations += ne
        out += res
        if debug_limit is not None and debug_count >= debug_limit:
            break
        debug_count += 1

    if len(out) == 0:
        print 'No results'
        return None, None

    # [namedtuple] to DataFrame
    df = pd.DataFrame.from_records(out, columns=out[0]._fields,
                                   index=['model', 'example_paper'])
    return df, not_enough_combinations

def find_lethal_reactions(series, models, max_combinations,
                          print_solution=False):
    def get_gr(model):
        sol = model.optimize()
        return 0.0 if sol.f is None else sol.f

    out = []; could_not_finds = []; not_enough_combinations = []

    try:
        setup = setup_for_series(series, models, True)
        design_no_del = Design(setup.design.heterologous_pathway,
                                [],
                                setup.design.target_exchange)
        # set up
        apply_design(setup.model, design_no_del, setup.use_greedy_knockouts)
        apply_environment(setup.model, setup.environment)
        m = setup.model
    except SetUpModelError as e:
        return None, None

    # knock out reactions associated with genes
    all_reactions = set()
    for g in set(setup.design.gene_knockouts):
        try:
            to_remove = [x.id for x in m.genes.get_by_id(g).reactions]
        except KeyError:
            could_not_finds.append('(Could not find gene %s in model %s)' % \
                                    (g, series['model']))
            continue
        all_reactions = all_reactions.union(to_remove)

    # try reaction combinations
    count = 1; found_lethal = False
    while not found_lethal and count <= max_combinations:
        for rs in it.combinations(all_reactions, count):
            m2 = m.copy()
            for r in rs:
                m2.reactions.get_by_id(r).knock_out()
            gr = get_gr(m2)
            if gr < min_biomass:
                found_lethal = True
                if print_solution:
                    print '%s, %s: %.3f' % (series['model'], '+'.join(rs), gr)
                out.append(LethalInteraction(series['model'], rs, gr,
                                             series['paper'], setup.design.gene_knockouts))
        count += 1

    # not enough combinations
    if not found_lethal:
        if print_solution:
            print '%s: Could not find a lethal combination with combinations of %d' % \
                (series['model'], max_combinations)
        not_enough_combinations.append([series['model'], series['paper'], setup.design.gene_knockouts])

    if print_solution:
        if len(could_not_finds) > 0:
            print
        for f in could_not_finds:
            print f
    return out, not_enough_combinations
