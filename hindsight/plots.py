# -*- coding: utf-8 -*-

from me_scripts.hindsight.hindsight import models_to_compare

import numpy as np
from collections import namedtuple

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
