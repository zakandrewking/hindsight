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

# uncategorized in ME
uncategorized_papers = [
    'Atsumi2008-zd_2', 'Atsumi2008-zd_3', 'Blankschien2010-hm',
    'Chang1999-qi_1', 'Chang1999-qi_2', 'Dien2001-gq_1',
    'Dien2001-gq_2', 'Donnelly1998-sb', 'Fong2005-ty_2', 'Hu2010-if',
    'Iverson2013-ng', 'Iverson2016-nn', 'Jian2010-ex_1',
    'Jian2010-ex_2', 'Jian2010-ex_3', 'Jung2010-aw_2', 'Jung2011-dx_2',
    'Kim2016-dy', 'Lee2005-xo_2', 'Lim2013-sy', 'Ma2013-nl',
    'Maeda2007-oh', 'Mazumdar2010-pt', 'Park2012-rs', 'Portnoy2008-ce',
    'Saini2014-br', 'Sanchez2005-cz', 'Sanchez2005-mu', 'Shen2013-uo',
    'Singh2011-am', 'Stols1997-kf', 'Stols1997-yg', 'Tran2014-mt',
    'Trinh2011-rx', 'Vemuri2002-xk', 'Volker2014-pd', 'Wang2012-ib',
    'Yan2009-ij_1', 'Zhang2010-ud', 'Zhang2011-zm', 'Zhou2011-nt',
]

# make samples for each ME simulation
df_me = df.loc[idx[uncategorized_papers, 'ME'], :].reset_index(level='model')
df_samples = generate_sampling_problems(df_me, 200)

# run parallel
apply_slurm(df_samples, 'hindsight', 'sample_series', directory)
