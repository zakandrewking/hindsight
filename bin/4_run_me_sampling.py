#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Sample kinetic parameters for the ME model."""

from hindsight import generate_sampling_problems, sample_series
from parallel_pandas import apply_slurm

import pandas as pd
idx = pd.IndexSlice
from os import getenv, mkdir
from os.path import join

df = pd.read_pickle('../data/sims_table.pickle')
directory = join(getenv('SCRATCH'), 'hindsight_sampling')
try:
    mkdir(directory)
except OSError:
    pass

# make 500 sims for each ME simulation
df_me = df.loc[idx[:, 'ME'], :].reset_index(level='model')
df_samples = generate_sampling_problems(df_me, 500)

# run parallel
apply_slurm(df_samples, 'hindsight', 'sample_series', directory)
