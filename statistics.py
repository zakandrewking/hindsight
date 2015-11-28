
import pandas as pd
import numpy as np
import itertools
import math
import re
import brewer2mpl
from os import listdir
from os.path import join
from warnings import warn

import matplotlib.pylab as plt
from scipy.stats import f_oneway, norm
from scipy.linalg import svd

try:
    from me.problem_and_solver_classes.me_solution import LP_Solution
except ImportError:
    warn('Could not import LP_Solution')

DEFAULT_COLORS = brewer2mpl.get_map('Dark2', 'Qualitative', 8).hex_colors
CARBON_CONTENT = {'Succinate': 4, 'Hexanoate (n-C6:0)': 6, 'Pyruvate': 3, 'CO2': 1,
                  'Succinate': 3, 'D-Glucose': 6, 'Dihydroxyacetone': 3}
    
def get_id(k):
    matches = re.findall(r'([0-9]+)\.json', k)
    return int(matches[0]) if len(matches) > 0 else None
    
def by_carbon(s, sub='D-Glucose', carbons=CARBON_CONTENT):
    """Normalize all entries in series s by the substrate uptake s[sub], and by the
    number of carbons, defined in dictionary `carbons`.

    """
    return pd.Series({k: (float(v) / abs(s[sub]) * carbons[k] / carbons[sub])
                      for k, v in s.iteritems()})

def run_and_plot_svd(a_df, n=4):
    U, s, Vh = svd(a_df.values)
    projected = np.dot(a_df.values, Vh)
    comb = list(itertools.combinations(range(n), 2))
    nrows = int(math.ceil(float(len(comb))/3))
    fix, axes = plt.subplots(ncols=3, nrows=nrows, figsize=(15,4*nrows))
    for ax, (i, j) in zip(axes.flatten(), comb):
        ax.scatter(projected[:,i], projected[:,j], c='#eeaaaa', alpha=0.5, s=0.2*projected.shape[0])
        ax.set_xlabel('PC%d' % i)
        ax.set_ylabel('PC%d' % j)
    return projected, U, s, Vh

def print_top_pc(s, a_Vh, a_df, n=3):
    """Identify key components of the secretion profile.

    """
    for i in range(3):
        print ", ".join(": ".join(['%s (%.2g)' % (y, x)]) 
                        for x, y in sorted(zip(a_Vh[i], list(a_df.columns)), 
                                           key=lambda x: abs(x[0]), reverse=True) if abs(x) > 0.05)
        print
    fig, ax = plt.subplots()
    ax.bar(range(s.size), s, log=True)
    _=ax.set_title('Singular values of secretion mix')

def plot_principle_components(proj, labels, n=None, colors=DEFAULT_COLORS):
    fix, axes = plt.subplots(ncols=3, nrows=1, figsize=(10,2))
    for ax, (i, j) in zip(axes.flatten(), [(0,1), (0,2), (1,2)]):
        ax.scatter(proj[:, i], proj[:, j],
                   c=[colors[x] if x < 8 else 'w' for x in labels], alpha=0.2, s=0.2*proj.shape[0])
        ax.set_xlabel('PC%d' % i)
        ax.set_ylabel('PC%d' % j)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if n is not None and (i, j) == (0, 2):
            ax.set_title('n = %d' % n)

def generate_angles(df):
    div = lambda d, k1=None, k2=None: (d[k1] / d[k2] if d[k2] != 0 else 0)
    angle = lambda s: math.atan(s) # abs(math.atan(s) - math.pi/4)
    ss = {'_d_'.join([k1, k2]): df.apply(div, k1=k1, k2=k2, axis=1).apply(angle)
          for k1, k2 
          in itertools.combinations(df.columns, 2)}
    ndf = pd.DataFrame(data=ss)
    return ndf.fillna(ndf.max().max())

def filter_out_unsampled(df, threshold=0.01):
    new_df = df[df.apply(lambda vals: abs((max(vals) - min(vals)) / max(vals)) >= threshold, axis=1)]
    print '%d of %d were sampled by > %d%%' % (len(new_df), len(df), threshold*100)
    return new_df

def run_anova(df, group_inds, p_thresholds=[1e-2, 1e-3, 1e-4]):
    sig_keffs = {}
    s = df.shape[1]
    for n in range(50, s, 50) + [s]:
        # get the relevant indices
        this_df = df.iloc[:,:n]
        # get the groups
        def anova_on_df(ser):
            return f_oneway(*[ser[g].dropna()
                              for g in group_inds.itervalues()])[1]
        p = this_df.apply(anova_on_df, axis=1)
        # sort
        p = p.sort(ascending=True, inplace=False)
        sig_keffs.update({(thresh, n): p[p < thresh] for thresh in p_thresholds})
        # print {k: len(v) for k, v in sig_keffs.iteritems()}
    return pd.DataFrame(sig_keffs)

def plot_significant_keffs(df):
    fix, ax = plt.subplots()
    # for each threshold
    for p in df.T.index.levels[0]:
        this_df = df.T.xs(p, level=0)
        ax.plot(this_df.index, this_df.count(axis=1))
        ax.set_xlabel('number of samples')
        ax.set_ylabel('number of significant keffs')

def plot_significant_shared(df, group_size=3):
    fix, ax = plt.subplots()
    # for each threshold
    out = {}
    for p in df.T.index.levels[0]:
        xs = []; ys = []
        for inds in [[i+x for x in range(group_size)] for i in range(len(df.columns.levels[1]) - group_size + 1)]:
            find_shared = lambda ser: len(set(ser.values))
            xs.append(df.columns.levels[1][max(inds)])
            ys.append(df[[p]].iloc[:,inds].dropna(how='any', axis=0).shape[0])
            out[(p, xs[-1])] = list(df[[p]].iloc[:,inds].dropna(how='any', axis=0).index)
        ax.plot(xs, ys)
        ax.set_xlabel('number of samples')
        ax.set_ylabel('number of significant')
        ax.set_title('Significant keffs shared with last %d groups' % group_size)
    return out

def plot_norm(data):
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(data, bins=40, histtype='bar', normed=True,
                               color=DEFAULT_COLORS[1], alpha=0.4) # '#DD8888'
    print 'log median: %s' % np.median(data)
    loc, scale = norm.fit(data)
    print 'loc: %.2f, scale: %.2f' % (loc, scale)
    dist = norm(loc=loc, scale=scale)
    pdf = dist.pdf(bins)
    ax.plot(bins, pdf, color=DEFAULT_COLORS[2], linewidth=3) # '#333333'
    
def kcats_for_reaction_list(mee, reactions, ec_numbers):
    """Get kcats for the list of reactions.

    Arguments
    ---------

    mee:

    reactions:

    ec_numbers:

    
    Returns
    -------

    (kcats, fit)

    kcats: DataFrame of kcat values with the given levels.

    fit: loc and shape for the best normal fit of the log(data), and the unlogged mean.

    """
    out = {}
    for reaction in reactions:
        try:            
            ec = map(int, ec_numbers[reaction].split('.'))
        except KeyError:
            print 'No ec number for %s' % reaction
            continue
        level = len(ec)
        while level >= 0:
            this_ec = mee.copy()
            for i in range(level):
                this_ec = this_ec[this_ec['EC%d'%(i+1)] == ec[i]]
            this_ec = this_ec.dropna(axis=0, how='any', subset=['kcat']).kcat
            if len(this_ec) >= 3:
                out[(reaction, level)] = pd.Series(list(this_ec))
                break
            level -= 1
    kcats = pd.DataFrame(out)
    for l, name in zip(kcats.columns.levels, ['reaction', 'level']):
        l.name = name
    fit = pd.DataFrame(columns = kcats.columns)
    fit.loc['loc'] = kcats.apply(lambda x: norm.fit(x.dropna().map(np.log))[0], axis=0)
    fit.loc['shape'] = kcats.apply(lambda x: norm.fit(x.dropna().map(np.log))[1], axis=0)
    fit.loc['mean'] = fit.loc['loc'].map(np.exp)
    return kcats, fit

def io_dataframe_for_dir(d, print_shape=False):
    """Return a dataframe of all io for JSON solution files in the given directory.
    
    """
    loaded = {}
    for directory in [d]:
        for path in listdir(directory):
            if not path.endswith('.json'): continue
            loaded[path] =  LP_Solution.import_from_json(join(directory, path))
    io = {get_id(k): s.get_io() for k, s in loaded.iteritems() if s.primal_dictionary is not None}
    df = pd.DataFrame(io).T
    df = df.fillna(0)
    df = df.drop(None, axis=1)
    if print_shape:
        print '%d samples and %s exchanges' % df.shape

        # pull out all the keffs by row
    keffs = {get_id(k): s.lp_problem_definition.basic_model_parameters.nondefault_enzyme_efficiency_in_per_second
             for k, s in loaded.iteritems() if s.lp_problem_definition.basic_model_parameters is not None}
    keffs_df = pd.DataFrame(keffs)
    keffs_log_df = keffs_df.applymap(np.log)
    return df, keffs_df, keffs_log_df
