# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import json
from os.path import join
import random
import logging

def sample_series(series, data_directory):
    outfile = join(data_directory, '%s_%s_%d_output.json' % (series['paper'],
                                                             series['model'],
                                                             series['index']))
    results = pd.Series({
        'paper': series['paper'],
        'model': series['model'],
        'growth_rate': growth_rate,
        'flux': flux,
        'keffs': keffs,
    })
    results.to_json(outfile)

def generate_sampling_problems(df, number_samples, random_seed=None,
                               distribution='lognormal', center_for_mu=True,
                               mu=2.48, sigma=3.29, multiplier=5,
                               normalize_by_median=False, start_with=0,
                               nondefault_fit=None):
    """Beginning with a DataFrame of simulation descriptions, generate sampling
    problems for each and return a new DataFrame with a MultiIndex that has a
    new level 'sample'.

    Arguments
    ---------

    df: The DataFrame.

    number_samples: Number of sampling problems for each input.

    random_seed: set to a known value to ensure the same keffs between runs.

    distribution: 'lognormal' or 'uniform' distribution.

    center_for_mu: use log(current value) to center distribution

    mu: mu for lognormal distribution.

    sigma: sigma for lognormal distribution.

    multiplier: multiplier for uniform distribution

    normalize_by_median: for uniform distribution

    start_with:

    nondefault_fit: A dictionary with keys corresponding to keffs and values
    that are tuples defining the distribution shape (loc, shape, ...). Ignored
    if distribution is 'uniform'.

    Warnings
    --------

    TODO: test that random_seed makes for a deterministic sim.
    TODO: account for enzyme size in lognorm sampling

    """

    # set seed
    if random_seed:
        random.seed(random_seed)

    # get the starting keffs
    keffs_start = get_keffs_start()

    # just sample all of them
    keffs_params = {k: (v, True) for k, v in keffs_start.iteritems()}

    # generate the samples
    if distribution == 'lognormal':
        samples = [get_lognormal_keff(keffs_params, mu, sigma,
                                      center_for_mu=center_for_mu,
                                      nondefault_fit=nondefault_fit)
                   for _ in range(number_samples)]
    elif distribution == 'uniform':
        samples = [get_random_keff(keffs_params, multiplier,
                                   normalize_by_median=normalize_by_median)
                   for _ in range(number_samples)]
    else:
        raise Exception('Bad distribution %s. Must be lognormal or uniform' % distribution)

    # make a new DataFrame with the combined MultiIndex
    samples_df = pd.DataFrame([samples for _ in df.iterrows()], index=df.index)
    samples_df.columns.name = 'sample'
    stack = samples_df.stack()
    stack.name = 'keffs'
    # merge them and set the index
    index_names = df.index.names
    merged = df.reset_index().merge(stack.reset_index()).set_index(index_names + ['sample'])

    return merged

def get_lognormal_keff(keffs_params, mu, sigma, center_for_mu=False,
                       nondefault_fit=None):
    """Will return randomly sampled keffs within input dictionary bounds, using a
    normal distribution.

    Arguments
    ---------

    keffs_params: dictionary of tuples (keff, should_sample).

    mu: mu

    sigma: sigma

    center_for_mu: use log(current value) to center distribution

    nondefault_fit: A dictionary with keys corresponding to keffs and values
    that are tuples defining the distribution shape (loc, shape).

    """
    def lognorm_sample(name, value):
        if nondefault_fit is not None and name in nondefault_fit:
            loc, shape = nondefault_fit[name]
        elif value == 0:
            logging.debug('keff == zero: %s' % name)
            return 0
        elif center_for_mu:
            loc, shape = np.log(value), sigma
        else:
            loc, shape = mu, sigma
        return random.lognormvariate(loc, shape)

    # sample values where should_sample is true
    return {k: (lognorm_sample(k, v[0]) if v[1] else v[0])
            for k, v in keffs_params.iteritems()}

def get_random_keff(keffs_params, multiplier, normalize_by_median=False):
    """Will return randomly sampled keffs with uniform distribution.

    keffs_params: dictionary of tuples (keff, should_sample).

    multiplier: Min and max determined my dividing and multiplying by 5.

    normalize_by_median: If true, then the keffs will be normalized by the
                         median of the centers of the bound ranges.

    """
    constant_enzyme_eff = {}
    sampled_enzyme_eff = {}
    sampled_centers = []
    for k, v in keffs_params.iteritems():
        if (v[1]):
            sampled_enzyme_eff[k] = random.uniform(v[0]/multiplier, v[0]*multiplier)
            sampled_centers.append((v[0]/multiplier + v[0]*multiplier) / 2.0)
        else:
            constant_enzyme_eff[k] = v[0]

    if normalize_by_median:
        keff_per_weight = (np.median(sampled_centers) *
                           len(sampled_enzyme_eff) /
                           sum(sampled_enzyme_eff.values()))
        sampled_enzyme_eff = dict([(k, v * keff_per_weight)
                                   for k, v in sampled_enzyme_eff.iteritems()])

    constant_enzyme_eff.update(sampled_enzyme_eff)
    return constant_enzyme_eff

def get_keffs_start():
    with open('../data/keffs_default.json', 'r') as f:
        return json.load(f)
