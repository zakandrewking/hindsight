#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight import carbon_yield
from hindsight.pathways import add_all_heterologous_pathways
from theseus.bigg import download_model

import pandas as pd
idx = pd.IndexSlice
import numpy as np
from sys import argv, exit
from os import listdir
from os.path import join

iJO1366_heterologous = add_all_heterologous_pathways(download_model('iJO1366'))

def get_fl(d, k):
    return d.get(k, 0.0) if d else None

def get_fluxes(ser, df):
    """Get the target exchange fluxes for a series."""
    if ser['model'] != 'ME':
        return np.nan
    if ser['paper'] not in df.index.levels[0]:
        return np.nan
    fluxes = [{'growth_rate': row['growth_rate'],
               'exchange_flux': get_fl(row['metabolic_flux'], ser['target_exchange']),
               'exchange_yield': carbon_yield(iJO1366_heterologous, ser['target_exchange'], ser['substrate_exchanges'], ser['supplement_exchanges'], row['metabolic_flux'])}
              for _, row in df.xs(ser['paper']).iterrows()]
    return fluxes

def main(directory):
    # load the sampling results
    gen = (pd.read_json(join(directory, filename), typ='series')
           for filename in listdir(directory)
           if filename.endswith('.json'))
    df = pd.DataFrame(gen).set_index(['paper', 'sample'])

    # load the sims table to get target exchanges
    samp = pd.read_pickle('../data/sims_table.pickle')

    # get the exchange fluxes
    samp['target_exchange_fluxes'] = (samp
                                      .reset_index()
                                      .apply(get_fluxes, args=(df,), axis=1)
                                      .values)

    # save
    samp.loc[:, ['target_exchange_fluxes']].to_pickle('../data/sampling_table.pickle')

if __name__ == '__main__':
    if len(argv) <= 1:
        print 'Usage: ./5_collect_me_sampling_results.py results_dir'
        exit()
    main(argv[1])
