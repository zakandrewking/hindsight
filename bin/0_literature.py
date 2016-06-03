#!/usr/bin/env python
# -*- coding: utf-8 -*-

from hindsight.bibtex import (load_bib, db_to_df, process_df, get_category,
                              get_strategies, get_max_category,
                              get_max_strategies)

bib = load_bib('../data/hindsight.bib')
df_raw = db_to_df(bib)
df = process_df(df_raw, True)

print get_category(get_max_category(df, 'strategies'), 'strategies')
print get_category(get_max_category(df, 'simulation'), 'simulation')
print get_strategies(get_max_strategies(df))

dups = (df[df.duplicated(['citation_key', 'target', 'deletions', 'substrate', 'additions'])]
        .loc(axis=1)['target', 'citation_key', 'additions', 'deletions'])
assert len(df[df.citation_key.isin(dups.citation_key)].loc(axis=1)['target', 'citation_key', 'additions', 'deletions', 'substrate']) == 0

print 'Number of papers: %d' % len(df.citation_key.unique())
print 'Number of designs: %d' % len(df.drop_duplicates(['citation_key', 'deletions', 'substrate', 'additions']))
print 'Number of targets: %d' % len(df.target.unique())

df.to_pickle('../data/literature_table.pickle')

key_cols = ['citation_key', 'year', 'target', 'simulation', 'strain', 'parent',
            'additions', 'deletions', 'aerobicity', 'substrate', 'evolved',
            'byproduct_order', 'strategies', 'see_study', 'gene_details']
df.loc[:, key_cols].to_csv('../data/literature_table_summary.tsv', sep='\t', encoding='utf8', index=False)
