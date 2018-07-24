# -*- coding: utf-8 -*-

from hindsight import models_to_compare
from hindsight.variables import (min_biomass, growth_coupled_cutoff,
                                 growth_coupled_h2_cutoff)

import pandas as pd
idx = pd.IndexSlice
import numpy as np
from collections import namedtuple
import cPickle as pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


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
        sims = sims[sims.loc[:, 'target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupled(models, [0]*len(models), target, 0)
    design_total = len(sims.xs(models[0], level='model'))
    col = ('yield_min' if use_yield else 'min')
    greater_than = (sims.loc[:, col] > threshold).reset_index()
    greater_than = greater_than[greater_than.model.isin(models)]
    is_growth_coupled = greater_than.groupby('model').apply(lambda x: len(x[x[col] == True]))
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

# -------------------------------------------------------------------------------
# Model GrowthCoupled Gradient
# -------------------------------------------------------------------------------

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
        sims = sims[sims.loc[:, 'target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupledGradient(models, [[]]*len(models), target, 0)
    design_total = len(sims.xs(models[0], level='model'))
    # get lists for the yield values
    yield_lists = (sims
                   .loc[:, ('yield_max' if use_max_yield else 'yield_min')]
                   .groupby(level='model')
                   .agg(lambda x: list(x.sort_values(inplace=False, ascending=False))))
    yield_lists = yield_lists.loc[models]
    models, counts = zip(*yield_lists.iteritems())
    if len(models) != len(counts):
        raise Exception('Counts not available for all models')
    return ModelGrowthCoupledGradient(models, counts, target, design_total,
                                      use_max_yield)

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
# Category maps
# -------------------------------------------------------------------------------

CategoryMapData = namedtuple('CategoryMapData', ['x_labels',
                                                 'y_labels',
                                                 'categories',
                                                 'title'])

def get_tree_table(sims):
    def _load_pickle(path):
        with open(path, 'r') as f:
            return pickle.load(f)
    secretion_tree_results = _load_pickle('../data/secretion_tree/secretion_tree_results.pickle')
    secretion_tree_w_gene_kos_results = _load_pickle('../data/secretion_tree_w_gene_kos/secretion_tree_results.pickle')
    tree_table = (sims
                  .merge(secretion_tree_results, left_index=True, right_index=True, how='left')
                  .merge(secretion_tree_w_gene_kos_results, left_index=True, right_index=True, how='left', suffixes=['', '_w_gene_kos'])
                  .loc[:, ['can_secrete', 'can_secrete_w_gene_kos']])
    return tree_table

def get_sampling_table(sims):
    """Get the metabolic flux of the target exchange."""
    samp = pd.read_pickle('../data/sampling_table.pickle')
    return (sims
            .merge(samp, how='left', left_index=True, right_index=True)
            .loc[:, ['target_exchange_fluxes', 'target_exchange']])

# Secretion Tree categories:
#
#                secretion tree:       y                    n             nan   (all not growth coupled and not non-unique)
# with kos: y                     insufficient            ERROR          ERROR
#           n                      detrimental       model limitation    ERROR
#           nan                       ERROR               ERROR           none

def _ijo_to_me(ser):
    ser.loc[idx[:, 'ME']] = ser.loc[idx[:, 'iJO1366']].values
    return ser

def _category_model_limitation(sims, tree_table=None, **kwargs):
    """This category means that the target cannot be growth-coupled with any set of
    iterative fermentation KOs."""
    model_limitation = ((tree_table['can_secrete'] == False) &
                        (tree_table['can_secrete_w_gene_kos'] == False))
    return 'Target byproduct cannot be growth coupled in the model (after exhaustive search)', _ijo_to_me(model_limitation)

def _category_insufficient_knockouts(sims, tree_table=None, **kwargs):
    """This category means that the target can be growth-coupled with some set of
    iterative fermentation KOs starting from the designed strain KOs."""
    insufficient = ((tree_table['can_secrete'] == True) &
                    (tree_table['can_secrete_w_gene_kos'] == True))
    return 'Experimental KO(s) insufficient for in silico growth coupling', _ijo_to_me(insufficient)

def _category_detrimental_knockouts(sims, tree_table=None, **kwargs):
    """This category means that the target cannot be growth-coupled with any set of
    iterative fermentation KOs starting from the designed strain KOs."""
    detrimental = ((tree_table['can_secrete'] == True) &
                   (tree_table['can_secrete_w_gene_kos'] == False))
    return 'Experimental KO(s) prohibit in silico growth coupling', _ijo_to_me(detrimental)

def _category_error(sims, tree_table=None, **kwargs):
    """An error occurred somewhere along the way."""
    error = (
        ((tree_table['can_secrete'] == False) & (tree_table['can_secrete_w_gene_kos'] == True)) |
        ( tree_table['can_secrete'].isnull() & ~tree_table['can_secrete_w_gene_kos'].isnull()) |
        (~tree_table['can_secrete'].isnull() &  tree_table['can_secrete_w_gene_kos'].isnull())
    )
    return 'ERROR', error

def _category_me_uncategorized(sims, **kwargs):
    return 'Uncategorized ME simulations', sims.index.map(lambda x: x[1] == 'ME')

def _category_uncategorized(sims, **kwargs):
    return 'Target byproduct is not growth coupled (not lethal in silico)', ~sims.target_exchange.isnull()

def _category_parameterization(sims, sampling_table=None, **kwargs):
    def any_g_coupled(ser):
        l = ser['target_exchange_fluxes']
        if type(l) is float and np.isnan(l):
            return False
        else:
            return any((d['growth_rate'] >= min_biomass and
                        (d['exchange_yield'] > growth_coupled_cutoff or
                         (ser['target_exchange'] == 'EX_h2_e' and d['exchange_flux'] > growth_coupled_h2_cutoff)))
                       for d in l)
    return (
        'With some parameter sets, target byproduct is growth coupled (>{:.0f}% C yield)*'.format(growth_coupled_cutoff * 100),
        sampling_table.apply(any_g_coupled, axis=1)
    )

def _category_non_unique(sims, cutoff=0.15, **kwargs):
    non_unique = ((sims.loc[:, 'yield_min'] < 0.8 * cutoff) &
                  (sims.loc[:, 'yield_max'] >= cutoff))
    return 'Alternative optimal growth coupled solutions', non_unique

def _category_lethal(sims, cutoff=min_biomass, **kwargs):
    # deal iwth nan's
    gr = (sims.loc[:, 'growth_rate'].fillna(0) < cutoff)
    return 'Experimental KO(s) are lethal in silico ', gr

def _category_growth_coupled(sims, **kwargs):
    growth_coupled = (
        (sims.loc[:, 'yield_min'] >= growth_coupled_cutoff) |
        ((sims.loc[:, 'target_exchange'] == 'EX_h2_e') & (sims.loc[:, 'min'] >= growth_coupled_h2_cutoff))
    )
    return 'Target byproduct is growth coupled (>{:.0f}% C yield)*'.format(growth_coupled_cutoff * 100), growth_coupled

# low to high priority
_category_fns = [
    # _category_error,
    # _category_insufficient_knockouts,
    # _category_detrimental_knockouts,
    # _category_model_limitation, # TODO make sure none of these overlap with
    #                             # insufficient_knockouts
    # _category_me_uncategorized,
    _category_uncategorized,
    _category_lethal, # could also be isozymes
    _category_non_unique, # could also be insufficient knockouts
    _category_parameterization,
    _category_growth_coupled, # highest priority
]
_none_category_name = 'Uncategorized'

def calculate_category_map(df, sort='year', title='Failure model categories',
                           models=models_to_compare):

    # fill in an empty row
    categories = pd.DataFrame({'category': _none_category_name}, index=df.index)
    # apply the functions
    tree_table = get_tree_table(df)
    sampling_table = get_sampling_table(df)
    for category_fn in _category_fns:
        name, result = category_fn(df, tree_table=tree_table,
                                   sampling_table=sampling_table)
        categories[result] = name
    categories = (categories
                  .reset_index()
                  .set_index(['model', 'paper', 'year'])
                  .loc[:, 'category']
                  .unstack('model')
                  .loc[:, models]
                  .sortlevel('year', ascending=False))
    y_labels = ["%s" % x[0] for x in list(categories.index)]
    return CategoryMapData(models, y_labels, categories.values, title)

def _category_colors(unique_categories):
    palette = sns.color_palette('Paired', len(unique_categories) + 1, desat=0.9)
    # #         greens             blues         dark orange      purple          red           pink
    # palette = palette[3:1:-1] + palette[:2] + palette[-2:-1] + palette[-1:] + palette[5:6] + palette[4:5]
    palette = palette[3:1:-1] + palette[1::-1] + palette[4:5]

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
        sims = sims[sims.loc[:, 'target'] == target]
    if len(sims) == 0:
        return ModelGrowthCoupledCategories(models, [[]]*len(models), target, [], 0)
    design_total = len(sims.xs(models[0], level='model'))
    # fill in an empty row
    categories_df = pd.DataFrame({'category': _none_category_name}, index=sims.index)
    # apply the functions
    tree_table = get_tree_table(sims)
    sampling_table = get_sampling_table(sims)
    cats = [] # keep track of category order
    for category_fn in _category_fns:
        name, result = category_fn(sims, tree_table=tree_table,
                                   sampling_table=sampling_table)
        categories_df[result] = name
        cats.append(name)
    cats.reverse()

    # Comment out to ignore Uncategorized:
    show_none_category = False
    if show_none_category:
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
                .groupby(level='model')
                .value_counts()
                .unstack(level=1)
                .T
                .loc[category_list, models]
                .fillna(0))
    count_arrays = count_df.values
    categories = list(count_df.index)

    return ModelGrowthCoupledCategories(models, count_arrays, target, categories,
                                        design_total)

def rename_models(model_list, year=True):
    repl = {
        'e_coli_core': 'Core model\n(2006)' if year else 'Core model',
        'iJR904': 'iJR904\n(2003)' if year else 'iJR904',
        'iAF1260': 'iAF1260\n(2007)' if year else 'iAF1260',
        'iAF1260b': 'iAF1260b\n(2010)' if year else 'iAF1260b',
        'iJO1366': 'iJO1366\n(2011)' if year else 'iJO1366',
        'ME': 'iOL1650-ME\n(2013)' if year else 'iOL1650-ME',
    }
    return [repl[x] if x in repl else x for x in model_list]

def plot_model_growth_coupled_categories(model_growth_coupled_categories,
                                         axis=None, show_legend=True,
                                         gc_line=False, y_lim_max_label=False):
    if axis is None:
        _, axis = plt.subplots(figsize=(9, 5))

    palette, _, _, _ = _category_colors(model_growth_coupled_categories.categories)

    model_range = range(len(model_growth_coupled_categories.models))
    bottom = np.zeros(len(model_growth_coupled_categories.models))
    all_handles = []; all_labels = []
    level_0 = None
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
            level_0 = level
        if i == 1 and gc_line:
            axis.plot([x + 0.5 for x in model_range], [sum(x) for x in zip(level, level_0)], label=None,
                      linewidth=3, linestyle='dotted', marker='s', color=(0.1, 0.2, 0.1, 1),
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
    axis.set_xticklabels(rename_models(model_growth_coupled_categories.models))
    axis.set_ylim(0, total)
    if y_lim_max_label:
        yticks = [0, total]
    else:
        yticks = np.arange(0, total + 1, 10)
    axis.set_yticks(yticks)
    axis.set_yticks(np.arange(0, total + 1, 1), minor=True)
    if model_growth_coupled_categories.target == 'all':
        axis.set_title('Simulations of all literature strains (n={})'.format(total), y=1.08)
    else:
        axis.set_title(('{} (n = {})' .format(model_growth_coupled_categories.target, total)),
                       y=1.08)
    axis.set_ylabel('Number of designs')
    axis.set_xlabel('Model')

    return axis

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
                                                 'label'])
ModelComparison.__repr__ = lambda self: '<ModelComparison %s, %s>' % (self.paper, str(self.models))

def get_model_comparison(sims, paper, models=models_to_compare, plot_yield=False):
    paper_sims = sims.xs(paper, level=0).loc[models]
    if plot_yield:
        min_secretions = paper_sims['yield_min'].fillna(0).values
        max_secretions = paper_sims['yield_max'].fillna(0).values
    else:
        min_secretions = paper_sims['min'].fillna(0).values
        max_secretions = paper_sims['max'].fillna(0).values
    growth_rates = paper_sims['growth_rate'].fillna(0).values
    target_exchange = paper_sims.iloc[0].target_exchange
    return ModelComparison(models, paper, min_secretions, max_secretions,
                           growth_rates, target_exchange, plot_yield,
                           paper)

def plot_model_comparison(axis, model_comparison,
                          colors=sns.color_palette('muted', 2)):
    """Plot production variability and growth rates for each model.

    """
    xs = np.array(range(len(model_comparison.models)))

    # # plot growth_rates
    # ax2 = axis.twinx()
    # gr_mask = np.isfinite(model_comparison.growth_rates)
    # ax2.plot(xs[gr_mask], model_comparison.growth_rates[gr_mask],
    #          c=colors[1], label='Growth rates')
    # ax2.set_ylabel('Growth rate', rotation=-90, labelpad=20)
    # max_val = max(model_comparison.growth_rates[gr_mask])
    # ax2.set_ylim(-0.04 * max_val, 1.2 * max_val)
    # ax2.grid(False)
    # ax2.get_yaxis().set_tick_params(which='both', direction='out', length=5,
    #                                 colors=colors[1])
    # for tl in ax2.get_yticklabels():
    #     tl.set_color(colors[1])

    # plot production min and max
    y1_mask = np.isfinite(model_comparison.max_secretions)
    y2_mask = np.isfinite(model_comparison.min_secretions)
    axis.plot(xs[y1_mask], model_comparison.max_secretions[y1_mask],
              c=colors[0], label='Max optimal production')
    axis.plot(xs[y2_mask], model_comparison.min_secretions[y2_mask],
              linestyle='--', c=colors[0], linewidth=6, label='Min optimal production')
    axis.set_xticklabels([''] + rename_models(model_comparison.models, year=False))
    axis.set_xlabel('Model')
    axis.set_ylabel('Product yield' if model_comparison.is_yield else
                    'Production rate')
    axis.set_xlim(-0.2, xs[-1] + 0.2)
    max_val = max(model_comparison.max_secretions[y1_mask])
    axis.set_ylim(-0.04 * max_val, 1.2 * max_val)

    axis.legend(loc='upper left')
