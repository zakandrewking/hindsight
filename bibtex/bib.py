# -*- coding: utf-8 -*-

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import author, convert_to_unicode
import pandas as pd
import numpy as np
import re

from me_scripts.hindsight import (exchange_for_metabolite_name,
                                  add_all_heterologous_pathways,
                                  no_route_exchanges)

from theseus.bigg import download_model

# ----------------------------------
# BibTeX
# ----------------------------------

def strip_chars(s):
    return s.replace('{_}', '_').replace('{&}', '&')

def customizations(record):
    """Use some functions delivered by the library."""
    # record = type(record)
    record = author(record)
    # record = editor(record)
    # record = journal(record)
    # record = keyword(record)
    # record = link(record)
    # record = page_double_hyphen(record)
    # record = doi(record)
    record = convert_to_unicode(record)
    record['annote'] = strip_chars(record['annote'])
    return record

def load_bib(filepath):
    parser = BibTexParser()
    parser.customization = customizations
    with open(filepath, 'rb') as f:
        return bibtexparser.load(f, parser=parser)

# ----------------------------------
# DataFrame
# ----------------------------------

split_text = '---xxx---'
translate_keys = {
    'ID':      'design',
    'annote':  'notes',
    'pmid':    'pmid',
    # 'doi':     'doi', # parser does not yet recognize doi in BibTeX file
    'year':    'year',
    'author':  'authors',
}
note_keys = ['target', 'additions', 'deletions', 'strain', 'parent', 'aerobicity',
             'simulation', 'evolved', 'see_study', 'substrate',
             'yield', 'g_byproduct_order', 'molar_byproduct_order',
             'byproduct_order', 'strategies', 'gene_details']
IJO1366 = download_model('iJO1366')
IJO1366_HETEROLOGOUS = add_all_heterologous_pathways(IJO1366.copy())
NAME_TO_BNUM = {g.name: g.id for g in IJO1366.genes}
def get_b_genes(gene_names):
    return [NAME_TO_BNUM[name] if name in NAME_TO_BNUM else name
            for name in gene_names]

def duplicate_row(df, selector, num):
    row = df.loc[selector]
    new_index = '%s_%d' % (selector, num)
    return pd.concat([df, pd.DataFrame([row], index=[new_index])], axis=0), new_index

def extract_key(key, s, default=np.nan):
    l = re.findall(r'(?:^|\n)%s:[ \t]*(.*)[ \t]*(?:$|\n)' % key, s)
    if len(l) == 0:
        return default
    elif len(l) == 1:
        return l[0].strip()
    else:
        print('WARNING: Multiple matches for %s in %s' % (key, s))
        return l

def parse_order_string(s):
    out = {}
    for y in [x.strip().rsplit(None, 1) for x in s.split(',')]:
        if len(y) == 2:
            out[y[0].strip()] = float(y[1].strip())
        else:
            raise Exception('Bad order string: %s' % s)
    return out

def get_byproduct_names(series):
    if series.g_byproduct_order not in [np.nan, 'nd']:
        return parse_order_string(series.g_byproduct_order).keys()
    if series.molar_byproduct_order not in [np.nan, 'nd']:
        return parse_order_string(series.molar_byproduct_order).keys()
    return []

def no_val(v):
    return (type(v) is float and np.isnan(v)) or v == 'nd'

def get_metabolite_for_exchange(reaction_id):
    return IJO1366_HETEROLOGOUS.reactions.get_by_id(reaction_id).metabolites.keys()[0]

def parse_byproduct_order(series):
    # functions
    g_to_mol = lambda g, MW: g / MW
    mol_to_cmol = lambda mol, cnum: mol * cnum
    def normalize(l):
        """Normalize by percentage or total."""
        m = sum([x[1] for x in l])
        return [(x[0], x[1]/m) for x in l]
    def sort(l):
        """Sort by el[1]"""
        return sorted(l, key=lambda x: x[1], reverse=True)

    # find in the model
    def get_mw_c_num(met_name):
        """Get reaction ID, molecular weight, and number of carbons"""
        rxn_id, _ = exchange_for_metabolite_name(met_name)
        if len(rxn_id) > 1:
            raise Exception('Not supporting multiple products')
        met = get_metabolite_for_exchange(rxn_id[0])
        MW = met.formula_weight
        cnum = met.elements['C']
        return rxn_id, MW, cnum

    # if no data, return that
    if all([no_val(series[col]) for col in ['g_byproduct_order', 'molar_byproduct_order', 'byproduct_order']]):
        return np.nan

    if not no_val(series.g_byproduct_order):
        out = []
        try:
            parsed = parse_order_string(series.g_byproduct_order)
        except Exception as e:
            print '\nERROR: %s (%s)' % (e, series.name)
            return 'ERROR'
        for k, v in parsed.iteritems():
            try:
                rxn_id, MW, cnum = get_mw_c_num(k)
            except Exception as e:
                print '\nERROR: %s (%s)' % (e, series.name)
                return 'ERROR'
            out.append((rxn_id, mol_to_cmol(g_to_mol(v, MW), cnum)))
        # normalize
        return sort(normalize(out))

    if not no_val(series.molar_byproduct_order):
        out = []
        try:
            parsed = parse_order_string(series.molar_byproduct_order)
        except Exception as e:
            print '\nERROR: %s (%s)' % (e, series.name)
            return 'ERROR'
        for k, v in parsed.iteritems():
            try:
                rxn_id, _, cnum = get_mw_c_num(k)
            except Exception as e:
                print '\nERROR: %s (%s)' % (e, series.name)
                return 'error'
            out.append((rxn_id, mol_to_cmol(v, cnum)))
        # normalize
        return sort(normalize(out))

    print '\nERROR: No byproduct order (%s)' % series.name
    return 'ERROR'

def native_non_native(t, model_native=None, model_non_native=None):
    try:
        target, _ = exchange_for_metabolite_name(t.lower())
    except KeyError:
        print 'Could not find %s in mets dictionary' % t
        return np.nan
    if len(target) > 1:
        raise Exception('Not supporting multiple targets')
    target = target[0]
    if target in model_native.reactions:
        return True
    elif target in model_non_native.reactions:
        return False
    else:
        print 'bad target: %s, %s' % (t, target)
        return np.nan

def db_to_df(bibtexparser_db):
    return (
        pd.DataFrame(bibtexparser_db.entries)
        .rename(columns=translate_keys)
        .loc[:, translate_keys.values()]
    )

def process_df(df, note_as_html=False):
    # add columns
    df['citation_key'] = df['design']
    df['first_author'] = df['authors'].map(lambda x: x[ 0])
    df['last_author']  = df['authors'].map(lambda x: x[-1])

    # check for citation keys
    assert(len(df[df.citation_key.isnull()]) == 0)

    # set index
    df = df.set_index('design')

    # take out html
    def replace_br(t):
        nn = t.replace('\n', ' ').replace('</div>', '').replace('<div>', '\n')
        br = re.sub(r'<br\s*/?>', '\n', nn)
        return re.sub('<[^>]+>', '', br)
    if note_as_html:
        df.loc[:, 'notes'] = df.loc[:, 'notes'].apply(replace_br)

    # split with notes
    for i, row in list(df.iterrows()):
        if split_text in row['notes']:
            for num, note in enumerate(row['notes'].split(split_text)):
                df, n = duplicate_row(df, i, num + 1)
                df.loc[n, 'notes'] = note
            df = df.drop(i)

    # get note columns
    for key in note_keys:
        df[key] = df['notes'].map(lambda x: extract_key(key, x))
        # check for missing notes
        col = df[key]
        try:
            empty_str = (col == '-') | (col.str.strip() == '') | (col.str.contains('ERROR'))
        except TypeError:
            empty_str = [False] * len(col)
        empty = col[col.isnull() | empty_str]
        if len(empty) > 0:
            print '\n%d rows with error or missing %s: %s' % (len(empty), key, ', '.join(map(str, empty.index)))

    # sort by year
    df = df.sort_values('year', ascending=True)

    # get byproduct order
    df['c_byproduct_order'] = df.apply(parse_byproduct_order, axis=1)

    # get b genes
    def split_on_comma_or_space(x):
        if x == 'none':
            return []
        return [y.strip() for y in re.split(r'\s*,\s*|\s+', x)]
    df['deletions_b'] = df['deletions'].map(split_on_comma_or_space).map(get_b_genes)

    # native versus nonnative
    df['native'] = df['target'].apply(native_non_native,
                                      model_native=IJO1366,
                                      model_non_native=IJO1366_HETEROLOGOUS)
    return df

def get_category(df1, cat):
    """Count the various strategies by row in the DataFrame."""
    # (1) cell to columns
    return (pd.DataFrame(df1.loc[:, cat].str.split(',').apply(pd.Series)).stack()
           # (2) strip
          .str.strip()
          # (3) count
          .value_counts())

def get_max_category(df2, cat):
    """Group the rows by paper (citation_key) and find all strategies for each one."""
    return pd.DataFrame(df2
            .groupby('citation_key')[cat]
            .agg(lambda x: reduce(lambda c, n: c.union([x.strip() for x in n.split(',')]), x, set()))
            .apply(lambda x: ', '.join(x)))

def get_strategies(df1):
    """Count the various strategies by row in the DataFrame."""
    # (1) cell to columns
    return (pd.DataFrame(df1.simulation.str.split(',').apply(pd.Series)).stack()
           # (2) strip
          .str.strip()
          # (3) count
          .value_counts())

def get_max_strategies(df2):
    """Group the rows by paper (citation_key) and find all strategies for each one."""
    return pd.DataFrame(df2
            .groupby('citation_key')
            .simulation
            .agg(lambda x: reduce(lambda c, n: c.union([x.strip() for x in n.split(',')]), x, set()))
            .apply(lambda x: ', '.join(x)))
