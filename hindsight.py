from me_scripts.hindsight.pathways import (add_heterologous_pathway,
                                           get_substrate_dictionary,
                                           get_product_dictionary, get_designs,
                                           exchange_for_metabolite_name)
from me_scripts.hindsight.variables import (get_min_biomass, get_max_our,
                                            NotFoundError, SetUpModelError)

from theseus import carbons_for_exchange_reaction, load_model, setup_model
from me.get_db_cursor import get_db_cursor

import numpy as np
import pandas as pd
import cPickle as pickle
from collections import namedtuple, defaultdict
import itertools
from itertools import izip, combinations
from copy import copy as copy_fn
import cobra.io
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import json
from os.path import join
from warnings import warn
from scipy.stats import spearmanr
import sys

try:
    from ipy_progressbar import ProgressBar
except ImportError:
    warn('ipy_progressbar not installed')

try:
    from cobra.manipulation import delete_model_genes
    from cobra.manipulation.delete import find_gene_knockout_reactions
except ImportError:
    warn('Using old version of COBRApy')
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from me_scripts.general_plots import (calculate_production_envelope,
                                          plot_production_envelope)
except ImportError:
    warn('Matplotlib not installed')
try:
    import seaborn as sns
except ImportError:
    warn('Seaborn not installed')


models_to_compare = ['E. coli core', 'iJR904', 'iAF1260', 'iAF1260b', 'iJO1366', 'ME']
m_models_to_compare = ['E. coli core', 'iJR904', 'iAF1260', 'iAF1260b', 'iJO1366']


def get_yield_palette(infeasible_color=(0.7, 0.7, 0.7, 1.0)):
    color_palette = sns.cubehelix_palette(9, start=1.1, gamma=1.0, rot=0.7,
                                          hue=1.0, dark=0.2, light=0.9)
    cmap = (mpl
            .colors
            .LinearSegmentedColormap
            .from_list('cubehelix', [list(infeasible_color)] + color_palette))
    bounds = [-1.1, -0.01, 0.01] + [(x + 1) / 5.0 for x in range(5)]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return color_palette, cmap, bounds, norm


def plot_yield_colorbar(ax, cmap, bounds, norm):
    ax.set_axis_bgcolor((0, 1, 0, 0))
    cb2 = mpl.colorbar.ColorbarBase(ax,
                                    cmap=cmap,
                                    norm=norm,
                                    boundaries=bounds,
                                    spacing='uniform',
                                    orientation='vertical')
    cb2.set_ticks([-0.545, 0] + bounds[2:])
    cb2.set_ticklabels(['infeasible', 0] + ['%d %%' % (x*100) for x in bounds][2:])


def print_single_kos(sim_setup):
    print 'Substrate: %s\n' % sim_setup.case.substrate_exchange
    for r in sim_setup.case.reaction_knockouts:
        m = sim_setup.wildtype_model.copy()
        reaction = m.reactions.get_by_id(r)
        reaction.knock_out()
        print r, m.optimize(), reaction.gene_reaction_rule


def get_secretion(model, flux_dict, flux_threshold=1e-3, sort=True):
    """Get all secretions with fluxes greater than flux_threshold."""
    secretions = [(k, v) for k, v in flux_dict.iteritems()
                  if 'EX_' in k and v > flux_threshold and
                  (carbons_for_exchange_reaction(model.reactions.get_by_id(k)) > 0 or k == 'EX_h2_e')]
    if sort:
        return sorted(secretions, key=lambda x: x[1], reverse=True)
    else:
        return secretions


def save_case_for_escher(name, model, reaction_data,
                         directory='/Users/zaking/lab/hindsight/figs/maps/'):
    cobra.io.save_json_model(model, join(directory, '%s_model.json' % name))
    with open(join(directory, '%s_rdata.json' % name), 'w') as f:
        json.dump(reaction_data, f)


def get_absolute_max(model, target_exchange_id, additions,
                     min_biomass=get_min_biomass()):
    """Find the absolute max production, and check for usage of the nonnative pathways."""

    # set min biomass
    biomass_reaction = [x for x in model.reactions
                        if x.objective_coefficient!=0][0]
    biomass_reaction.lower_bound = min_biomass

    model.change_objective(target_exchange_id)
    solution = model.optimize()
    if solution.f is None:
        return np.nan, np.nan

    absolute_max_product = model.solution.f

    # get flux through non-native pathways
    designs = get_designs()
    if additions in designs:
        reactions = designs[additions][1]
        additions_pathway_flux = {r: solution.x_dict[r] for r in reactions
                                  if not r.startswith('EX_')}
    else:
        additions_pathway_flux = np.nan

    return absolute_max_product, additions_pathway_flux


def minimize_maximize(model, target_exchange_id, substrate_reaction_ids,
                      min_biomass=get_min_biomass(), biomass_fraction=0.9999,
                      calculate_yield=True):
    """Perform a single flux variablity analysis."""
    biomass_reaction = [x for x in model.reactions
                        if x.objective_coefficient!=0][0]

    model.optimize()
    if model.solution.f is None:
        gr = 0
        yield_minmax = (np.nan, np.nan)
        minmax = (np.nan, np.nan)
        max_secretion = np.nan
        max_flux = np.nan
    elif model.solution.f < min_biomass:
        gr  = model.solution.f
        yield_minmax = (np.nan, np.nan)
        minmax = (np.nan, np.nan)
        max_secretion = np.nan
        max_flux = np.nan
    else:
        gr = model.solution.f
        # set the minimum biomass production
        biomass_reaction.lower_bound = biomass_fraction * gr
        # max and min
        model.change_objective(target_exchange_id)
        model.optimize(objective_sense='minimize')
        minmax = (model.solution.f,)
        model.optimize(objective_sense='maximize')
        minmax = minmax + (model.solution.f,)
        # get secretions
        max_secretion = get_secretion(model, model.solution.x_dict)
        max_flux = model.solution.x_dict
        if calculate_yield:
            # get sur
            if not isinstance(substrate_reaction_ids, list):
                substrate_reaction_ids = [substrate_reaction_ids]
            sur = abs(model.reactions.get_by_id(substrate_reaction_ids[0]).lower_bound)
            # get target carbons
            carbons = carbons_for_exchange_reaction(model.reactions.get_by_id(target_exchange_id))
            # loop through targets and add carbons
            substrate_carbon_uptake = 0
            for substrate_reaction_id in substrate_reaction_ids:
                c = carbons_for_exchange_reaction(model.reactions.get_by_id(substrate_reaction_id))
                substrate_carbon_uptake = substrate_carbon_uptake + (c * sur)
                yield_minmax = tuple([x*carbons/(substrate_carbon_uptake) for x in minmax])
        else:
            yield_minmax = (np.nan, np.nan)
    return gr, minmax, yield_minmax, max_secretion, max_flux


def get_gene_names(cursor, b_list):
    results = cursor.execute("select id, synonym from genome_region_synonym where id in (%s);" % ",".join(["'"+x+"'" for x in b_list]))
    out = {}
    for x in cursor.fetchall():
        if x[0] not in out: out[x[0]] = []
        out[x[0]].append(x[1])
        return out


def normalize_by_carbons(l, model=None):
    try:
        not_in_model = [(a,b) for a, b in l if a not in model.reactions]
        # if len(not_in_model) > 0:
        # print 'Ignoring reaction not in %s: %s' % (model, not_in_model)
        l = [(a, b*carbons_for_exchange_reaction(model.reactions.get_by_id(a))) for a, b in l
             if a in model.reactions]
        m = sum([x[1] for x in l])
        return [(x[0], float(x[1])/m) for x in l]
    except:
        return l


def get_plasmid_translation_fluxes(model, flux_dictionary):
    """Read the model, containing a plasmid, and return flux values for all
    plasmid translation elongation reactions in the flux_dictionary.

    """
    translation = (model
                   .plasmid
                   .get_heterologous_peptide_translation_reactions(get_db_cursor()))
    elongation = [x for x in translation[0].keys() if 'elongation' in x]
    return {k: (flux_dictionary[k] if k in flux_dictionary else 0)
            for k in elongation}


def _get_design_exchanges():
    return set(itertools.chain.from_iterable(([[r for r in x[1].iterkeys() if 'EX_' in r]
                                               for x in get_designs().itervalues()])))

def me_exchange_for_metabolite_name(name, source=True, cursor=None):
    if cursor is None:
        cursor = get_db_cursor()

    if name in _get_design_exchanges():
        return name

    d = {'EX_glc_e': 'D-Glucose', 'EX_xyl__D_e': 'D-Xylose', 'EX_o2_e': 'O2',
         'EX_glyc_e': 'Glycerol', 'EX_sucr_e': 'Sucrose', 'EX_ac_e': 'Acetate',
         'EX_pyr_e': 'Pyruvate', 'EX_h2_e': 'H2', 'EX_for_e': 'Formate', 'EX_but_e': 'Butyrate (n-C4:0)',
         'EX_succ_e': 'Succinate', 'EX_ac_e': 'Acetate', 'EX_lac__D_e': 'D-Lactate',
         'EX_lac__L_e': 'L-Lactate', 'EX_etoh_e': 'Ethanol', 'EX_ala__L_e': 'L-Alanine',
         'EX_mal__L_e': 'L-Malate'}
    for k, v0 in d.iteritems():
        v = v0 + ('_source' if source else '_sink')
        if k.lower() == name.lower().strip():
            cursor.execute("select reaction_name from reaction where reaction_name = '%s'" % v)
            if cursor.fetchone() is None:
                raise Exception('bad exchange: %s' % name)
            return v
    raise Exception('bad exchange: %s' % name)


def knock_out_reactions(model, reaction_list, copy=True, leaky=False):
    if copy: model = model.copy()
    if not hasattr(model, 'knockouts'):
        model.knockouts = set()
    model.knockouts = model.knockouts.union(reaction_list)
    if leaky:
        cutoff = 0.1
        for r in reaction_list:
            reaction = model.reactions.get_by_id(r)
            reaction.lower_bound = -cutoff
            reaction.upper_bound = cutoff
    else:
        for r in reaction_list:
            model.reactions.get_by_id(r).knock_out()
    return model


def solve_and_print_secretion(model):
    sol = optimize_minimal_flux(model)
    print 'Objective: %.3f' % sol.f
    print 'Secretion: %s' % ', '.join(['%s: %.3f' % x for x in
                                       get_secretion(model, sol.x_dict)])
    return sol


def find_designs_to_couple(sims, paper, target_exchange='EX_succ_e',
                           models=m_models_to_compare):
    out = []
    for model in models:
        found = None
        # set up the sim
        setup = set_up_simulation(sims, paper, model, print_summary=False)
        # see if secretion is possible from the design
        res = find_alternate_optimal_secretions(setup.design_model, only_optimal=False,
                                                unselected_reactions=['EX_co2_e', target_exchange])
        for r in res:
            for rxn, val in r.secretion:
                if rxn == target_exchange and val > 1:
                    found = ('design', paper, model, res[0], r)
                    break
            if found: break
        if found:
            out.append(found)
            continue
        res = find_alternate_optimal_secretions(setup.wildtype_model, only_optimal=False,
                                                unselected_reactions=['EX_co2_e', target_exchange])
        for r in res:
            for rxn, val in r.secretion:
                if rxn == target_exchange and val > 1:
                    found = ('wildtype', paper, model, res[0], r)
                    break
            if found: break
        if found:
            out.append(found)
        else:
            out.append(('none', paper, model))
    return out


# -------------------------------------------------------------------------------
# For production envelopes:
# -------------------------------------------------------------------------------

SimulationSetup = namedtuple('SimulationSetup', ['case', # a row in the simulations DataFrame
                                                 'wildtype_model',
                                                 'design_model'])

def get_case(sims, paper, model='all', columns='all', first=False,
             print_summary=False):
    o = sims.xs(paper, level='Paper')
    if model != 'all':
        o = o.xs(model, level='Model')
    if columns != 'all':
        o = o.loc[:, columns]
    p = o.iloc[0]
    if first:
        o = p
    if print_summary:
        print('Target: %s' % p.Target)
        print('Target exchange: %s' % p.target_exchange)
        print('Substrate: %s' % p.Substrate)
        print('Substrate exchange: %s' % p.substrate_exchange)
        print('Deletions: %s' % p.Deletions)
    return o


def get_target_exchange_id(model, target):
    """Get the target exchange reaction."""
    exchange_id = exchange_for_metabolite_name(target, get_product_dictionary())

    # make sure the exchange is now in the model
    try:
        model.reactions.get_by_id(exchange_id)
    except KeyError:
        raise NotFoundError('target %s not in model %s' % (exchange_id, str(model)))

    return exchange_id


def get_model(model, additions, substrate, target, aerobicity, deletions_b,
              gene_ko=False, copy=False):
    """Make a model based on the inputs from notes."""
    if copy:
        model = model.copy()

    # add non-native pathway
    try:
        model = add_heterologous_pathway(model, additions)
    except NotFoundError:
        raise SetUpModelError('bad addition: %s' % additions)
    # add fhl for hydrogen production
    if target == 'H2':
        model = add_heterologous_pathway(model, 'fhl')

    # get the target
    try:
        target_exchange_id = get_target_exchange_id(model, target)
    except NotFoundError:
        raise SetUpModelError('bad target: %s' % target)

    # get the substrate
    try:
        substrates = exchange_for_metabolite_name(substrate, get_substrate_dictionary())
    except NotFoundError:
        raise SetUpModelError('bad substrate: %s' % substrate)

    # get the supplementations
    if isinstance(substrates, dict):
        substrate_exchange_ids = substrates['substrates']
        supplementation_exchange_ids = substrates['supplementations']
    else:
        substrate_exchange_ids = substrates
        supplementation_exchange_ids = []

    # check substrate(s)
    if model.id != 'Ecoli_core_model':
        # don't check supplementations in core model
        try:
            [model.reactions.get_by_id(x) for x in supplementation_exchange_ids]
        except KeyError:
            raise SetUpModelError('supplementations %s not in model %s' % (supplementation_exchange_ids, str(model)))
    try:
        if isinstance(substrate_exchange_ids, list):
            [model.reactions.get_by_id(x) for x in substrate_exchange_ids]
        else:
            model.reactions.get_by_id(substrate_exchange_ids)
    except KeyError:
        raise SetUpModelError('substrates %s not in model %s' % (substrate_exchange_ids, str(model)))

    # assign aerobicity and sur
    if aerobicity.strip().lower()=='aerobic':
        aerobic = True
        max_our = get_max_our()
    elif aerobicity.strip().lower()=='anaerobic':
        aerobic = False
        max_our = 0
    elif aerobicity.strip().lower()=='microaerobic':
        aerobic = False
        max_our = 0
    else:
        raise SetUpModelError('bad aerobicity %s' % aerobicity)

    # if there are multiple substrates, divide evenly
    sur = 10
    if isinstance(substrate_exchange_ids, list):
        sur = sur / len(substrate_exchange_ids)

    specific_bounds = {}
    # get the supplementations
    for supp in supplementation_exchange_ids:
        specific_bounds[supp] = (-10, 0)

    # bounds for H2
    if target_exchange_id == 'EX_h2_e':
        if 'FHL' in model.reactions:
            specific_bounds['FHL'] = (0, 1000)

    # set up the model
    model = setup_model(model, substrate_exchange_ids, aerobic=aerobic, sur=sur,
                        max_our=max_our)
    for k, v in specific_bounds.iteritems():
        try:
            model.reactions.get_by_id(k).lower_bound = v[0]
            model.reactions.get_by_id(k).upper_bound = v[1]
        except KeyError:
            SetUpModelError('reaction %s not found in model %s' % (k, str(model)))

    # knockouts
    genes = deletions_b
    reaction_knockouts = set()
    not_in_model = set()
    gene_reaction_knockouts = set()
    for g in genes:
        try:
            gene_obj = model.genes.get_by_id(g)
        except KeyError:
            not_in_model.add(str(g))
            continue
        if gene_ko:
            r = gene_obj.remove_from_model(model)
            reaction_knockouts = reaction_knockouts.union(r)
        else:
            reactions = gene_obj.reactions
            for r in reactions:
                r.lower_bound = 0; r.upper_bound = 0
            reaction_knockouts = reaction_knockouts.union([x.id for x in reactions])
            rg = [x.id for x in find_gene_knockout_reactions(model, [model.genes.get_by_id(g)])]
            gene_reaction_knockouts = gene_reaction_knockouts.union(rg)
    greedy_knockouts = reaction_knockouts.difference(gene_reaction_knockouts)

    return (model, substrate_exchange_ids, specific_bounds, target_exchange_id,
            reaction_knockouts, greedy_knockouts, not_in_model)

def set_up_simulation(sims, paper, model_name, gene_ko=False,
                      wildtype_add_pathway=True, print_summary=True):
    case = get_case(sims, paper, model_name, first=True,
                    print_summary=print_summary)
    model = load_model(model_name)

    wiltype_add = case['Additions'] if wildtype_add_pathway else 'none'
    out = get_model(model, wiltype_add, case['Substrate'], case['Target'],
                    case['Aerobicity'], [], gene_ko=gene_ko, copy=True)
    wildtype_model = out[0]

    out = get_model(model, case['Additions'], case['Substrate'], case['Target'],
                    case['Aerobicity'], case['Deletions_b'], gene_ko=gene_ko)
    design_model = out[0]
    return SimulationSetup(case, wildtype_model, design_model)

def calculate_envelopes(simulation_setup):
    target = simulation_setup.case.target_exchange
    envs = []
    for model_type in ['wildtype_model', 'design_model']:
        this_model = getattr(simulation_setup, model_type)
        # calculate the envelope
        env_data = calculate_production_envelope(this_model.copy(), target, get_min_biomass(),
                                                 label=model_type.split('_')[0])
        envs.append(env_data)
    return envs


def quick_production_envelope(model, target, min_biomass=0, axis=None,
                              plot_kwargs={}):
    p = calculate_production_envelope(model, target, min_biomass)
    ax = plot_production_envelope(p, axis=axis, plot_kwargs=plot_kwargs)
    if hasattr(model, 'knockouts'):
        ax.set_title('{}: Envelope for {}'.format(' '.join(model.knockouts), target))
    else:
        ax.set_title('Envelope for {}'.format(target))


# -------------------------------------------------------------------------------
# Min max vs. oxygen
# -------------------------------------------------------------------------------

MinMaxData = namedtuple('MinMaxData', ['min_values',
                                       'max_values',
                                       'our_range',
                                       'label'])

def calculate_min_max_vs_o2(simulation_setup, our_range=list(np.arange(0, 28, 2))):

    def calc(model, our_range, target_exchange):
        mint = []; maxt = []
        for our in our_range:
            m = model.copy()
            m.reactions.get_by_id('EX_o2_e').lower_bound = -our
            _, minmax, _, _, _ = minimize_maximize(m, target_exchange,
                                                   None, None,
                                                   calculate_yield=False)
            mint.append(minmax[0])
            maxt.append(minmax[1])
        return mint, maxt

    d = []
    for model_type in ['wildtype_model', 'design_model']:
        the_model = getattr(simulation_setup, model_type).copy()
        min_max = calc(the_model, our_range,
                       simulation_setup.case.target_exchange)
        d.append(MinMaxData(*min_max, our_range=our_range,
                            label=model_type.split('_')[0]))
    return d

def plot_min_max_vs_o2(axis, min_maxes):
    for min_max, color in izip(min_maxes, sns.color_palette("muted", 4)):
        axis.plot(min_max.our_range, min_max.min_values, color=color, linestyle='--')
        axis.plot(min_max.our_range, min_max.max_values, color=color, label=min_max.label)
    axis.set_xlabel("OUR")
    axis.set_ylabel("Production rate")
    axis.set_xlim(left=0, right=1.0*max(np.concatenate([min_max.our_range for min_max in min_maxes])))
    max_val = max(np.concatenate([min_max.max_values for min_max in min_maxes]))
    axis.set_ylim(-0.04 * max_val, 1.2 * max_val)
    axis.legend(loc='best')
    axis.set_title('Variability as a function of OUR')

# -------------------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------------------

ModelComparison = namedtuple('ModelComparison', ['models',
                                                 'paper',
                                                 'min_secretions', # np.array
                                                 'max_secretions', # np.array
                                                 'growth_rates', # np.array
                                                 'target_exchange',
                                                 'is_yield',
                                                 'cases',
                                                 'label'])
ModelComparison.__repr__ = lambda self: '<ModelComparison %s, %s>' % (self.paper, str(self.models))

def get_model_comparison(sims, paper, models=models_to_compare, plot_yield=False):
    cases = get_case(sims, paper, model='all')
    cases = cases.reset_index(level='year').loc[models, :]
    if plot_yield:
        min_secretions = cases['yield_min'].values
        max_secretions = cases['yield_max'].values
    else:
        min_secretions = cases['min'].values
        max_secretions = cases['max'].values
    growth_rates = cases['gr'].values
    target_exchange = cases.iloc[0].target_exchange
    return ModelComparison(models, paper, min_secretions, max_secretions,
                           growth_rates, target_exchange, plot_yield, cases,
                           paper)

try:
   default_colors = sns.color_palette("muted", 2)
except NameError:
    default_colors = ['r', 'g']
def plot_model_comparison(axis, model_comparison, colors=default_colors):
    """Plot production variability and growth rates for each model.

    """
    xs = np.array(range(len(model_comparison.models)))

    # plot growth_rates
    ax2 = axis.twinx()
    gr_mask = np.isfinite(model_comparison.growth_rates)
    ax2.plot(xs[gr_mask], model_comparison.growth_rates[gr_mask],
             c=colors[1], label='Growth rates')
    ax2.set_ylabel('Growth rate', rotation=-90, labelpad=20)
    max_val = max(model_comparison.growth_rates[gr_mask])
    ax2.set_ylim(-0.04 * max_val, 1.2 * max_val)
    ax2.grid(False)
    ax2.get_yaxis().set_tick_params(which='both', direction='out', length=5,
                                    colors=colors[1])
    for tl in ax2.get_yticklabels():
        tl.set_color(colors[1])

    # plot production min and max
    y1_mask = np.isfinite(model_comparison.max_secretions)
    y2_mask = np.isfinite(model_comparison.min_secretions)
    axis.plot(xs[y1_mask], model_comparison.max_secretions[y1_mask],
              c=colors[0], label='Maximum')
    axis.plot(xs[y2_mask], model_comparison.min_secretions[y2_mask],
              linestyle='--', c=colors[0], linewidth=6, label='Minimum')
    axis.set_xticklabels([''] + model_comparison.models)
    axis.set_xlabel('Model')
    axis.set_ylabel('Product yield' if model_comparison.is_yield else
                    'Production rate')
    axis.set_xlim(-0.2, xs[-1] + 0.2)
    max_val = max(model_comparison.max_secretions[y1_mask])
    axis.set_ylim(-0.04 * max_val, 1.2 * max_val)

    axis.legend(loc='best')

# -------------------------------------------------------------------------------
# Alternative solutions
# -------------------------------------------------------------------------------

class AlternativeOptimalSolution(namedtuple('AlternativeOptimalSolution', ['f',
                                                                           'knockouts',
                                                                           'secretion',
                                                                           'flux_threshold'])):
    """A data structure for alternative optimal secretion solutions.

    """
    def __repr__(self):
        return """kos = %s  (f=%.3f)

    -->    %s

""" % (', '.join(self.knockouts), self.f, ', '.join(['%s: %.2f' % x for x in self.secretion]))


def find_alternate_optimal_secretions(model, only_optimal=True,
                                      selected_reactions=[],
                                      unselected_reactions=[], flux_threshold=2,
                                      copy=True, loop_limit=100):
    """Return all alternate optimal secretions by interatively knocking out other
    options.

    """

    if copy:
        model = model.copy()

    out = []
    all_kos = []
    last_f = None
    l = 0
    while True:
        # KO and optimize
        if len(all_kos) > 0:
            model.reactions.get_by_id(all_kos[-1]).knock_out()
        soln = model.optimize()

        # check for break conditions
        if (soln.status == 'infeasible' or
            (only_optimal and last_f and abs(last_f - soln.f) > 0.01)):
            break

        # remember biomass objective
        last_f = soln.f

        # filter fluxes
        io = get_secretion(model, soln.x_dict, flux_threshold)
        out.append(AlternativeOptimalSolution(soln.f, copy_fn(all_kos), io, flux_threshold))

        # get the max secretion
        io_rxns = [x[0] for x in io]; found = False
        for i, selected_reaction in enumerate(selected_reactions):
            if selected_reaction in io_rxns and selected_reaction not in unselected_reactions:
                next_ko = selected_reaction
                selected_reactions.pop(i)
                found = True
        if not found:
            best = None
            for k, v in io:
                if (k not in unselected_reactions) and (best is None or v < best[1]):
                    best = (k, v)
            if best:
                next_ko = best[0]
            else:
                break
        all_kos.append(next_ko)

        # loop_limit
        if l > loop_limit:
            break
        l += 1
    return out

# -------------------------------------------------------------------------------
# Heat maps
# -------------------------------------------------------------------------------

HeatMapData = namedtuple('HeatMapData', ['x_labels',
                                         'y_labels',
                                         'values',
                                         'title'])


def calculate_heat_map(df, column='yield_min', sort='year',
                       title='Target secretion', models=models_to_compare):
    df = (df
          .reset_index()
          .set_index(['Model', 'Paper', 'year'])
          .loc[:, 'yield_min']
          .unstack('Model')
          .loc[:, models]
          .sortlevel('year', ascending=False))
    y_labels = ["%s" % x[0] for x in list(df.index)]
    return HeatMapData(models, y_labels, df.values, title)


def plot_heat_map(heat_map_data, figsize=[10,15]):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot2grid((1,8), (0, 0), colspan=7)
    ax2 = plt.subplot2grid((1,8), (0, 7), colspan=1)

    _, cmap, bounds, norm = get_yield_palette()

    plot_yield_colorbar(ax2, cmap, bounds, norm)

    # plot
    ax.pcolormesh(heat_map_data.values, cmap=cmap, norm=norm)

    # y axis
    ax.set_ylim(0, len(heat_map_data.y_labels))
    ax.set_yticks([x + 0.5 for x in range(len(heat_map_data.y_labels))])
    ax.set_yticklabels(heat_map_data.y_labels)

    # x axis
    ax.set_xlim(0, len(heat_map_data.x_labels))
    ax.set_xticks([x + 0.5 for x in range(len(heat_map_data.x_labels))])
    ax.set_xticklabels(heat_map_data.x_labels)

    ax.set_title(heat_map_data.title)


# -------------------------------------------------------------------------------
# Category maps
# -------------------------------------------------------------------------------

CategoryMapData = namedtuple('CategoryMapData', ['x_labels',
                                                 'y_labels',
                                                 'categories',
                                                 'title'])


def _category_model_limitation(sims):
    # TODO open data file 'cannot_be_coupled.json'
    model_limitation = (
        ((sims.index.get_level_values('Model') == 'E. coli core') & (sims.loc[:, 'Target'] == 'L-Lactate')) |
        (sims.index.get_level_values('Paper') == 'Atsumi2008_3') |
        (sims.index.get_level_values('Paper') == 'Dellomonaco2011') |
        (sims.index.get_level_values('Paper') == 'Atsumi2008_2') |
        (sims.index.get_level_values('Paper') == 'Jian2010_1') |
        (sims.index.get_level_values('Paper') == 'Jian2010_2') |
        (sims.index.get_level_values('Paper') == 'Jian2010_3') |
        (sims.index.get_level_values('Paper') == 'Jung2010_2') |
        (sims.index.get_level_values('Paper') == 'Jung2011_2') |
        (sims.index.get_level_values('Paper') == 'Zhang2011') |
        (sims.index.get_level_values('Paper') == 'Park2012e') |
        ((sims.index.get_level_values('Paper') == 'Wang2012a') & (sims.index.get_level_values('Model') == 'E. coli core')) |
        ((sims.index.get_level_values('Paper') == 'Wang2012a') & (sims.index.get_level_values('Model') == 'iAF1260b')) |
        ((sims.index.get_level_values('Model') == 'E. coli core') & (sims.index.get_level_values('Paper') == 'Ma2013')) |
        ((sims.index.get_level_values('Model') == 'E. coli core') & (sims.index.get_level_values('Paper') == 'Jian2010_1')) |
        ((sims.index.get_level_values('Paper') == 'Saini2014') & (sims.index.get_level_values('Model') == 'iAF1260b'))
        # no ME here
    ) & (sims.index.get_level_values('Model') != 'ME')
    return 'Cannot be growth coupled', model_limitation

# def _category_rate_yield(sims):
#     rate_yield = ((sims.index.get_level_values('Model') == 'ME') & (sims.loc[:, 'Target'] == 'Ethanol'))
#     return 'Rate-yield tradeoff', rate_yield

def _category_parameterization(sims):
    parameterization = (
        (sims.index.get_level_values('Paper') == 'Zhou2008') | # TODO check
        (sims.index.get_level_values('Paper') == 'Atsumi2008_3') |
        (sims.index.get_level_values('Paper') == 'Dellomonaco2011') |
        (sims.index.get_level_values('Paper') == 'Atsumi2008_2') |
        (sims.index.get_level_values('Paper') == 'Jian2010_1') |
        (sims.index.get_level_values('Paper') == 'Jian2010_2') |
        (sims.index.get_level_values('Paper') == 'Jian2010_3') |
        (sims.index.get_level_values('Paper') == 'Jung2010_2') |
        (sims.index.get_level_values('Paper') == 'Jung2011_2') |
        (sims.index.get_level_values('Paper') == 'Zhang2011') |
        (sims.index.get_level_values('Paper') == 'Park2012e') |
        (sims.index.get_level_values('Paper') == 'Stols1997a') |
        (sims.index.get_level_values('Paper') == 'Stols1997b') |
        (sims.index.get_level_values('Paper') == 'Donnelly1998') |
        (sims.index.get_level_values('Paper') == 'Vemuri2002') |
        (sims.index.get_level_values('Paper') == 'Lee2005a_1') |
        (sims.index.get_level_values('Paper') == 'Sanchez2005') |
        (sims.index.get_level_values('Paper') == 'Sanchez2005a') |
        (sims.index.get_level_values('Paper') == 'Jantama2008') |
        (sims.index.get_level_values('Paper') == 'Blankschien2010') |
        (sims.index.get_level_values('Paper') == 'Singh2011') |
        (sims.index.get_level_values('Paper') == 'Ma2013') |
        (sims.index.get_level_values('Paper') == 'Dellomonaco2011') |
        (sims.index.get_level_values('Paper') == 'Wang2012a') |
        (sims.index.get_level_values('Paper') == 'Lim2013a') |
        (sims.index.get_level_values('Paper') == 'Saini2014') |
        (sims.index.get_level_values('Paper') == 'Zelic2003') |
        (sims.index.get_level_values('Paper') == 'Trinh2011') |
        (sims.loc[:, 'Deletions'].apply(lambda s: 'dld' in [x.strip() for x in s.split(', ')])) |
        (sims.loc[:, 'Target'] == 'Ethanol')
    ) & (sims.index.get_level_values('Model') == 'ME')
    return 'ME-model parameterization', parameterization

def _category_insufficient_knockouts(sims):
    insufficient_knockouts = (
        (sims.index.get_level_values('Paper') == 'Chang1999_2') |
        (sims.index.get_level_values('Paper') == 'Dien2001_1') |
        (sims.index.get_level_values('Paper') == 'Dien2001_2') |
        (sims.index.get_level_values('Paper') == 'Mazumdar2014') |
        (sims.index.get_level_values('Paper') == 'Atsumi2008_1') |
        (sims.index.get_level_values('Paper') == 'Shen2013') |
        (sims.index.get_level_values('Paper') == 'Yan2009_1') |
        ((sims.index.get_level_values('Paper') == 'Stols1997a') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Stols1997b') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Donnelly1998') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Vemuri2002') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Lee2005a_1') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Sanchez2005') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Sanchez2005a') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Paper') == 'Ma2013') & (sims.index.get_level_values('Model') != 'ME')) |
        ((sims.index.get_level_values('Model') != 'ME') & (sims.loc[:, 'Target'] == 'D-Lactate')) |
        ((sims.index.get_level_values('Paper') == 'Wang2012a') & (sims.index.get_level_values('Model') == 'iJO1366')) |
        ((sims.index.get_level_values('Paper') == 'Lim2013a') & (sims.index.get_level_values('Model') == 'iJO1366')) |
        ((sims.index.get_level_values('Paper') == 'Saini2014') & (sims.index.get_level_values('Model') == 'iJO1366')) |
        ((sims.index.get_level_values('Paper') == 'Trinh2011') & (sims.index.get_level_values('Model') != 'ME')) |
        (sims.loc[:, 'Deletions'].apply(lambda s: 'dld' in [x.strip() for x in s.split(', ')]) & (sims.index.get_level_values('Model') != 'ME'))
    )
    return 'Insufficient or detrimental knockouts', insufficient_knockouts

def _category_non_unique(sims, cutoff=0.15):
    non_unique = ((sims.loc[:, 'yield_min'] < 0.8 * cutoff) &
                  (sims.loc[:, 'yield_max'] >= cutoff))
    return 'Alternative optima', non_unique

def _category_lethal(sims, cutoff=0.05):
    # deal iwth nan's
    gr = (sims.loc[:, 'gr'].fillna(0) < cutoff)
    return 'Lethal genotypes', gr

def _category_growth_coupled(sims, cutoff=0.15, h2_cutoff=2):
    growth_coupled = ((sims.loc[:, 'yield_min'] >= cutoff) |
                      (sims.loc[:, 'Target'] == 'H2') & (sims.loc[:, 'min'] >= h2_cutoff))
    return 'Growth coupled (>{:.0f}% C yield)*'.format(cutoff*100), growth_coupled

_category_fns = [
    _category_model_limitation,
    # _category_rate_yield,
    _category_parameterization,
    _category_insufficient_knockouts, # could also be rate-yield tradeoff. Can also mean too many KOs.
    _category_non_unique, # could also be insufficient knockouts
    _category_lethal, # could also be isozymes
    _category_growth_coupled, # highest priority
]
_none_category_name = 'Uncategorized'

def calculate_category_map(df, sort='year', title='Failure model categories',
                           models=models_to_compare):
    # fill in an empty row
    categories = pd.DataFrame({'category': _none_category_name}, index=df.index)
    # apply the functions
    for category_fn in _category_fns:
        name, result = category_fn(df)
        categories[result] = name
    categories = (categories
                  .reset_index()
                  .set_index(['Model', 'Paper', 'year'])
                  .loc[:, 'category']
                  .unstack('Model')
                  .loc[:, models]
                  .sortlevel('year', ascending=False))
    y_labels = ["%s" % x[0] for x in list(categories.index)]
    return CategoryMapData(models, y_labels, categories.values, title)


def _category_colors(unique_categories):
    palette = sns.color_palette('Paired', len(unique_categories) + 1, desat=0.9)
    palette = palette[3:4] + palette[:2] + palette[4:-1] + palette[-1:]+ palette[2:3]

    cmap = (mpl
            .colors
            .LinearSegmentedColormap
            .from_list('categories', palette))
    bounds = [x + 0.5 for x in range(len(unique_categories) + 1)]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return palette, cmap, bounds, norm


def plot_category_map(category_map_data, figsize=[10,15]):
    # fig
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot2grid((1,8), (0, 0), colspan=7)
    ax2 = plt.subplot2grid((1,8), (0, 7), colspan=1)

    # colors
    unique_categories = np.unique(category_map_data.categories)
    palette, cmap, bounds, norm = _category_colors(unique_categories)

    # nums
    category_lookup = {cat: i + 0.5 for i, cat in enumerate(unique_categories)}
    get_category_num = np.vectorize(lambda x: category_lookup[x])
    category_nums = get_category_num(category_map_data.categories)

    # plot
    ax.pcolormesh(category_nums, cmap=cmap, norm=norm)

    # y axis
    ax.set_ylim(0, len(category_map_data.y_labels))
    ax.set_yticks([x + 0.5 for x in range(len(category_map_data.y_labels))])
    ax.set_yticklabels(category_map_data.y_labels)

    # x axis
    ax.set_xlim(0, len(category_map_data.x_labels))
    ax.set_xticks([x + 0.5 for x in range(len(category_map_data.x_labels))])
    ax.set_xticklabels(category_map_data.x_labels)

    ax.set_title(category_map_data.title)

    # legend
    cb2 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
                                    norm=norm,
                                    boundaries=bounds,
                                    spacing='uniform',
                                    orientation='vertical')
    cb2.set_ticks([x + 0.5 for x in bounds])
    cb2.set_ticklabels(unique_categories)


ModelGrowthCoupledCategories = namedtuple('ModelGrowthCoupledCategories', ['models',
                                                                           'count_arrays', # len(count_arrays) == # categories
                                                                           'target', # can be 'all'
                                                                           'categories',
                                                                           'total']) # boolean
ModelGrowthCoupledCategories.__repr__ = lambda self: ('<ModelGrowthCoupledCategories('
                                                      'target={self.target}, models={self.models}, '
                                                      'total={self.total}, categories={self.categories}, '
                                                      'count_arrays={self.count_arrays})>'
                                                      .format(self=self))


def calculate_model_growth_coupled_categories(sims, target='all',
                                              models=models_to_compare,
                                              category_list=None,
                                              first_category=None,
                                              last_category=None):
    if target != 'all':
        sims = sims[sims.loc[:, 'Target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupledCategories(models, [[]]*len(models), target, [], 0)
    design_total = len(sims.xs(models[0], level='Model'))
    # fill in an empty row
    categories_df = pd.DataFrame({'category': _none_category_name}, index=sims.index)
    # apply the functions
    cats = [] # keep track of category order
    for category_fn in _category_fns:
        name, result = category_fn(sims)
        categories_df[result] = name
        cats.append(name)
    cats.reverse()
    if _none_category_name in list(categories_df.loc[:, 'category']):
        cats.append(_none_category_name)
    # list of ordered categories
    if not category_list:
        category_list = cats
        # category_list = list(categories_df.loc[:, 'category'].unique())
        if first_category:
            category_list.insert(0, category_list.pop(category_list.index(first_category)))
        if last_category:
            category_list.append(category_list.pop(category_list.index(last_category)))
    count_df = (categories_df
                .loc[:, 'category']
                .groupby(level='Model')
                .value_counts()
                .unstack(level=1)
                .T
                .loc[category_list, models]
                .fillna(0))
    count_arrays = count_df.values
    categories = list(count_df.index)

    return ModelGrowthCoupledCategories(models, count_arrays, target, categories,
                                        design_total)


def plot_model_growth_coupled_categories(model_growth_coupled_categories,
                                         axis=None, show_legend=True,
                                         gc_line=False, y_lim_max_label=False):
    if axis is None:
        _, axis = plt.subplots()

    palette, _, _, _ = _category_colors(model_growth_coupled_categories.categories)

    model_range = range(len(model_growth_coupled_categories.models))
    bottom = np.zeros(len(model_growth_coupled_categories.models))
    all_handles = []; all_labels = []
    for i, (category, level, color) in enumerate(zip(model_growth_coupled_categories.categories,
                                                     model_growth_coupled_categories.count_arrays,
                                                     palette)):
        bars = axis.bar([x + 0.1 for x in model_range], level,
                        bottom=bottom, color=color, label=category, linewidth=1)
        for bar in bars:
            # if i == 0:
                # bar.set_hatch('///')
            bar.set_edgecolor((0.3,0.3,0.3,1.0)) #[(x*1.3 if x*1.3 < 1 else 1) for x in color[:3]] + [1.0])
        if i == 0 and gc_line:
            axis.plot([x + 0.5 for x in model_range], level, label='Model accuracy',
                      linewidth=3, linestyle='-', marker='s', color=(0.1, 0.2, 0.1, 1),
                      markersize=10, fillstyle='full')

        bottom = bottom + level

    total = model_growth_coupled_categories.total

    if show_legend:
        # reverse the legend order
        handles, labels = axis.get_legend_handles_labels()
        handles.reverse(); labels.reverse()
        axis.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

    axis.set_xlim(0, model_range[-1] + 1)
    axis.set_xticks([x + 0.5 for x in model_range])
    axis.set_xticklabels(model_growth_coupled_categories.models)
    axis.set_ylim(0, total)
    if y_lim_max_label:
        yticks = [0, total]
    else:
        yticks = np.arange(0, total + 1, 10)
    axis.set_yticks(yticks)
    axis.set_yticks(np.arange(0, total + 1, 1), minor=True)
    if model_growth_coupled_categories.target == 'all':
        axis.set_title('Designs categories (n={})'.format(total), y=1.08)
    else:
        axis.set_title(('{} (n = {})' .format(model_growth_coupled_categories.target, total)),
                       y=1.08)
    axis.set_ylabel('Number of designs')
    axis.set_xlabel('Model')

    return axis


# ---------------------------------------------
# lethal genotypes
# ---------------------------------------------

LethalInteraction = namedtuple('LethalInteraction', ['model', # str
                                                     'reactions', # list of str IDs
                                                     'growth_rate',
                                                     'example_paper',
                                                     'example_genotype'])

def find_summary_lethal_reactions(sims, max_combinations=5,
                                  models=m_models_to_compare, debug_limit=None):
    """Returns a summary dataframe of all the lethal combinations and cases where
    max_combinations was not enough to find the lethal interaction.

    """
    to_string = lambda x: ','.join(sorted(x))
    out = [] # list of LethalInteraction's
    not_enough_combinations = []
    dc = 0
    for paper in ProgressBar(list(sims.reset_index()['Paper'].unique())):
        res, ne = find_lethal_reactions(sims, paper, max_combinations, models, False)
        not_enough_combinations += ne
        out += res
        if debug_limit is not None and dc >= debug_limit:
            break
        dc += 1
    df = pd.DataFrame.from_records(out, columns=out[0]._fields, index=['model', 'example_paper'])
    return df, not_enough_combinations

def find_lethal_reactions(sims, paper, max_combinations=3,
                          models=m_models_to_compare, print_solution=False,
                          gr_lim=1e-3):
    def get_gr(model):
        sol = model.optimize()
        return 0.0 if sol.f is None else sol.f

    out = []; could_not_finds = []; not_enough_combinations = []
    for model_name in models:
        case = get_case(sims, paper, model=model_name, first=True)
        model = load_model(model_name)
        try:
            m, _, _, _, _, _, _ = get_model(model, case.Additions, case.Substrate,
                                            case.Target, case.Aerobicity,
                                            case.Deletions_b, copy=True)
        except SetUpModelError as e:
            continue

        if get_gr(m) > gr_lim:
            continue

        # get the model with no deletions
        m, _, _, _, _, _, _ = get_model(model, case.Additions, case.Substrate,
                                        case.Target, case.Aerobicity, 'none')

        # knock out reactions associated with genes
        all_reactions = set()
        for g in set(case.Deletions_b):
            try:
                to_remove = [x.id for x in m.genes.get_by_id(g).reactions]
            except KeyError:
                could_not_finds.append('(Could not find gene %s in model %s)' % \
                                       (g, model_name))
                continue
            all_reactions = all_reactions.union(to_remove)

        # try reaction combinations
        count = 1; found_lethal = False
        while not found_lethal and count <= max_combinations:
            for rs in combinations(all_reactions, count):
                m2 = model.copy()
                for r in rs:
                    m2.reactions.get_by_id(r).knock_out()
                gr = get_gr(m2)
                if gr < gr_lim:
                    found_lethal = True
                    if print_solution:
                        print '%s, %s: %.3f' % (model_name, '+'.join(rs), gr)
                    else:
                        out.append(LethalInteraction(model_name, rs, gr, paper,
                                                     case.Deletions))
            count += 1

        # not enough combinations
        if not found_lethal:
            if print_solution:
                print '%s: Could not find a lethal combination with combinations of %d' % \
                    (model_name, max_combinations)
            else:
                not_enough_combinations.append([model_name, paper, case.Deletions])

    if print_solution:
        if len(could_not_finds) > 0:
            print
        for f in could_not_finds:
            print f
    else:
        return out, not_enough_combinations


# -------------------------------------------------------------------------------
# Model GrowthCoupled
# -------------------------------------------------------------------------------

ModelGrowthCoupled = namedtuple('ModelGrowthCoupled', ['models',
                                                       'counts',
                                                       'target', # can be 'all'
                                                       'total'])
ModelGrowthCoupled.__repr__ = lambda self: ('<ModelGrowthCoupled(target="%s", models=%s, total=%s, counts=%s)>' %
                                            (self.target, self.models, self.total, self.counts))


def calculate_model_growth_coupled(sims, target='all', threshold=1, use_yield=False,
                                   models=models_to_compare):
    if target != 'all':
        sims = sims[sims.loc[:, 'Target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupled(models, [0]*len(models), target, 0)
    design_total = len(sims.xs(models[0], level='Model'))
    col = ('yield_min' if use_yield else 'min')
    greater_than = (sims.loc[:, col] > threshold).reset_index()
    greater_than = greater_than[greater_than.Model.isin(models)]
    is_growth_coupled = greater_than.groupby('Model').apply(lambda x: len(x[x[col] == True]))
    is_growth_coupled = is_growth_coupled.loc[models]
    models, counts = zip(*is_growth_coupled.iteritems())
    if len(models) != len(counts):
        raise Exception('Counts not available for all models')
    return ModelGrowthCoupled(models, counts, target, design_total)


try:
   default_colors = sns.color_palette("muted", 2)
except NameError:
    default_colors = ['r', 'g']
def plot_model_growth_coupled(model_growth_coupled, mgc2=None, axis=None,
                              y_tick_step=5, color=default_colors, labels=None):
    if axis is None:
        _, axis = plt.subplots()

    model_range = range(len(model_growth_coupled.models))
    try:
        l1 =  labels[0]
    except TypeError, KeyError:
        l1 = None
    axis.bar([x + 0.1 for x in model_range], model_growth_coupled.counts,
             color=default_colors[0], label=l1)
    total = model_growth_coupled.total
    if mgc2 is not None:
        try:
            l2 =  labels[1]
        except TypeError, KeyError:
            l2 = None
        axis.bar([x + 0.1 for x in model_range], mgc2.counts,
                 bottom=model_growth_coupled.counts, color=default_colors[1],
                 label=l2)
        total += mgc2.total

    axis.set_xlim(0, model_range[-1] + 1)
    axis.set_xticks([x + 0.5 for x in model_range])
    axis.set_xticklabels(model_growth_coupled.models)
    axis.set_ylim(0, total)
    axis.set_yticks(np.arange(0, total + 1, y_tick_step))
    if model_growth_coupled.target == 'all':
        axis.set_title('Designs where target molecule is growth-coupled in silico (n = {})'.format(total),
                       y=1.08)
    else:
        axis.set_title('{} (n = {})'.format(model_growth_coupled.target, total),
                       y=1.08)
    axis.set_ylabel('Number of designs')
    axis.set_xlabel('Model')

    return axis


ModelGrowthCoupledGradient = namedtuple('ModelGrowthCoupled', ['models',
                                                               'count_arrays', # len(count_arrays) == total
                                                               'target', # can be 'all'
                                                               'total',
                                                               'max_yield']) # boolean
ModelGrowthCoupledGradient.__repr__ = lambda self: ('<ModelGrowthCoupledGradient(target="%s", models=%s, total=%s, count_arrays=%s)>' %
                                                    (self.target, self.models, self.total, self.count_arrays))


def calculate_model_growth_coupled_gradient(sims, target='all',
                                            models=models_to_compare,
                                            use_max_yield=False):
    if target != 'all':
        sims = sims[sims.loc[:, 'Target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupledGradient(models, [[]]*len(models), target, 0)
    design_total = len(sims.xs(models[0], level='Model'))
    # get lists for the yield values
    yield_lists = (sims
                   .loc[:, ('yield_max' if use_max_yield else 'yield_min')]
                   .groupby(level='Model')
                   .agg(lambda x: list(x.sort(inplace=False, ascending=False))))
    yield_lists = yield_lists.loc[models]
    models, counts = zip(*yield_lists.iteritems())
    if len(models) != len(counts):
        raise Exception('Counts not available for all models')
    return ModelGrowthCoupledGradient(models, counts, target, design_total,
                                      use_max_yield)


def plot_model_growth_coupled_gradient_colorbar(ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(3, 4))
    _, cmap, bounds, norm = get_yield_palette((0, 0, 0, 0.0))
    plot_yield_colorbar(ax, cmap, bounds, norm)


def plot_model_growth_coupled_gradient(model_growth_coupled_gradient, axis=None,
                                       y_tick_step=5):
    if axis is None:
        _, axis = plt.subplots()

    palette, cmap, bounds, norm = get_yield_palette((0, 0, 0, 0.0))
    sm_fn = mpl.cm.ScalarMappable(norm, cmap).to_rgba
    model_range = range(len(model_growth_coupled_gradient.models))
    for i, level in enumerate(zip(*model_growth_coupled_gradient.count_arrays)):
        bars = axis.bar([x + 0.1 for x in model_range], [1]*len(model_range),
                        bottom=[i]*len(model_range))
        for bar, val in zip(bars, level):
            if np.isnan(val):
                bar.set_color((0, 0, 0, 0))
            else:
                bar.set_color(sm_fn([val])[0])
                # bar.set_edgecolor((0.3, 0.3, 0.3))

    total = model_growth_coupled_gradient.total

    axis.set_xlim(0, model_range[-1] + 1)
    axis.set_xticks([x + 0.5 for x in model_range])
    axis.set_xticklabels(model_growth_coupled_gradient.models)
    axis.set_ylim(0, total)
    axis.set_yticks(np.arange(0, total + 1, y_tick_step))
    if model_growth_coupled_gradient.target == 'all':
        axis.set_title('Designs where target molecule is growth-coupled in silico (n = {})'.format(total),
                       y=1.08)
    else:
        axis.set_title(('{} {} (n = {})'
                        .format(model_growth_coupled_gradient.target,
                                ('max' if model_growth_coupled_gradient.max_yield else 'min'),
                                total)),
                       y=1.08)
    axis.set_ylabel('Number of designs')
    axis.set_xlabel('Model')

    return axis


# -------------------------------------------------------------------------------
# OptKnock design comparisons
# -------------------------------------------------------------------------------

def to_dataframe(ser, column_name, index_name, copy=False):
    """Make the series a dataframe with the given names for index and single column."""
    if copy:
        ser = ser.copy()
    ser.index.name = index_name
    return pd.DataFrame(ser, columns=[column_name]).reset_index()

def get_gpr_df(gene_counts_df, model):
    # make a gene-reaction dataframe
    gpr = pd.DataFrame(gene_counts_df.gene)
    # get reactions as lists
    gpr['reactions'] = gpr.gene.apply(lambda g: get_reactions_for_gene(g, model))
    # stack the genes
    gpr = gpr.set_index('gene')
    gpr = gpr['reactions'].apply(lambda x: pd.Series(x)).stack()
    gpr.index = gpr.index.droplevel(-1)
    # make it a dataframe
    gpr = pd.DataFrame(gpr, columns=['reaction']).reset_index()
    return gpr

def drop_na_reactions_or_single_row(df):
    """Drop nan reactions or reactions with 0 counts in gpr rows where a
    gene already has a matching reaction with counts.

    """
    new_df = df.dropna(subset=['reaction'])
    new_df = new_df[new_df.reaction_count_equiv > 0]
    if len(new_df) == 0:
        new_df = df.head(1)
        new_df['reaction'] = np.nan
        new_df['label'] = new_df['gene']
    return new_df

def get_matching_equiv(reaction, equivs):
    try:
        return (rs for rs in equivs if reaction in rs).next()
    except StopIteration:
        return np.nan

def get_reactions_for_gene(gene_id, a_model):
    try:
        gene = a_model.genes.get_by_id(gene_id)
    except KeyError:
        return []
    return [x.id for x in gene.reactions]

def calculate_reaction_gene_comparison(reaction_counts, gene_counts, model,
                                       equivalent_reactions=[('ACKr', 'PTAr')]):
    # (1) deal with equivalent reactions
    # make it a DataFrame
    reaction_counts = to_dataframe(reaction_counts, 'reaction_count', 'reaction')
    # add new rows for missing reactions
    all_reactions = (x for sublist in equivalent_reactions for x in sublist)
    missing_reactions = (x for x in all_reactions
                         if not reaction_counts['reaction'].apply(lambda y: x == y).any())
    reaction_counts = pd.concat([reaction_counts,
                                 pd.DataFrame({'reaction': list(missing_reactions),
                                               'reaction_count': 0})],
                                ignore_index=True)
    # find the equivalent reactions
    reaction_counts['equivalent_reactions'] = (reaction_counts['reaction']
                                               .apply(lambda r: get_matching_equiv(r, equivalent_reactions)))
    # add a new column for the summed counts
    groups = reaction_counts.groupby('equivalent_reactions').sum()
    reaction_counts = (reaction_counts
                       .merge(groups, how='outer', left_on='equivalent_reactions',
                              right_index=True, suffixes=['', '_equiv']))
    reaction_counts['reaction_count_equiv'].loc[reaction_counts['reaction_count_equiv'].isnull()] = reaction_counts['reaction_count']

    # (2) make the genes dataframe
    gene_counts.index.name = 'gene'
    gene_counts = pd.DataFrame(gene_counts, columns=['gene_count']).reset_index()

    # (3) get the gpr dataframe
    gpr = get_gpr_df(gene_counts, model)

    # (4) merge the dataframes
    result = (gpr.merge(gene_counts, how='outer', on='gene')
             .merge(reaction_counts[['reaction', 'reaction_count_equiv', 'equivalent_reactions']], how='outer', on='reaction'))
    result['label'] = result['equivalent_reactions'].fillna(result['reaction']).fillna('').apply(str) + ', ' + result['gene'].fillna('')

    # (5) drop nan reactions in gpr rows where a gene already has a matching reaction with counts
    result = result.groupby('gene').apply(lambda x: drop_na_reactions_or_single_row(x))

    return result

def plot_reaction_gene_comparison(comp, label_cutoffs=(25, 100), disp_fn=None):
    x, y, labels = izip(*comp[['gene_count', 'reaction_count_equiv', 'label']].fillna(0).itertuples(index=False))
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    for xa, ya, la in zip(x, y, labels):
        if xa > label_cutoffs[0] or ya > label_cutoffs[1]:
            xa += 0.2
            if disp_fn is not None:
                xa, ya, la = disp_fn(xa, ya, la)
            ax.text(xa, ya, la, fontsize=12)
    print 'r={} p={}'.format(*spearmanr(x, y))
    return ax


# -------------------------------------------------------------------------------
# Secretion trees
# -------------------------------------------------------------------------------

Secretions = namedtuple('Secretions', ['knockouts', # iterable
                                       'growth_rate', # float
                                       'exchange_reactions', # iterable
                                       'fluxes', # iterable with same order as exchange reactions
                                       ])
Secretions.__repr__ = lambda self: '<Secretions {self.knockouts}, {self.growth_rate:.3f}, {self.exchange_reactions}, {self.fluxes} >'.format(self=self)


class FoundReaction(Exception):
    pass


def secretions_for_knockouts(model, knockouts=[], max_depth=10, depth=0,
                             ignore_exchanges=[], raise_if_found=None,
                             growth_cutoff=get_min_biomass(), flux_cutoff=0.1,
                             debug=False):
    """Accepts a model and a set of knockouts.

    Returns a tree of secretions using nested dictionaries.

    Does not keep data if the raise_if_found option is True.

    Arguments
    ---------

    model: The cobra model.

    knockouts: A list of reaction IDs to knock out.

    max_depth: The maximum depth to search.

    depth:

    ignore_exchanges:

    raise_if_found:

    growth_cutoff:

    flux_cutoff:

    debug:

    """
    # check depth
    if depth > max_depth:
        return None

    # always copy the model
    model = model.copy()

    # knock out the reactions
    for ko in knockouts:
        model.reactions.get_by_id(ko).knock_out()

    # solve the problem
    solution = model.optimize()
    if solution.f is None or solution.f <= growth_cutoff:
        if debug:
            print knockouts, solution.f
            sys.stdout.flush()
        return None
    else:
        secretion = Secretions(knockouts, solution.f, *zip(*get_secretion(model, solution.x_dict, sort=False)))
        if debug:
            print secretion
            sys.stdout.flush()
        if raise_if_found:
            try:
                index = secretion.exchange_reactions.index(raise_if_found)
            except ValueError:
                pass
            else:
                if secretion.fluxes[index] > flux_cutoff:
                    raise FoundReaction(str(secretion))
        return {'data': (None if raise_if_found else secretion),
                'children': {new_knockout: secretions_for_knockouts(model, knockouts + [new_knockout],
                                                                    max_depth, depth + 1,
                                                                    ignore_exchanges, raise_if_found,
                                                                    growth_cutoff, flux_cutoff,
                                                                    debug)
                             for new_knockout in secretion.exchange_reactions
                             if new_knockout not in ignore_exchanges}}


def print_secretion(sec):
    return '[f: %.2f | %s]' % (sec.growth_rate,
                               ', '.join('%s: %.2f' % x for x in zip(sec.exchange_reactions, sec.fluxes)))


def print_tree(tree, depth=0, format_data=lambda x: str(x)):
    disp = ('   ' * depth)
    out = '--> ' + format_data(tree['data']) + '\n'
    for k, v in tree['children'].iteritems():
        if type(v) is dict:
            out += (disp + '|\n' + disp + ('\--%s%s' % (k, print_tree(v, depth + 1, format_data))))
    #out += '\n'
    return out
