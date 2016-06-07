#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Sample kinetic parameters for the ME model."""

from hindsight import generate_sampling_problems
from parallel_pandas import apply_slurm

import pandas as pd
from os import getenv, mkdir
from os.path import join

df = pd.read_pickle('../data/sims_table.pickle')
directory = join(getenv('SCRATCH'), 'hindsight_sampling')
try:
    mkdir(directory)
except OSError:
    pass

# make 500 sims
df = generate_sampling_problems(df, 500)

# run parallel
# apply_slurm(df, 'hindsight', 'sample_series', directory)
