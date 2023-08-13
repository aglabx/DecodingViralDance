import numpy as np
import pandas as pd
import functools as ft
import re
"""
This script will parse SNP count file, bin, and optionally map to genes


210_G_T 241_C_T         1524    1735    3209236
210_G_T 3037_C_T        1524    15166   3202595
210_G_T 3644_G_T        1524    18007   10908
210_G_T 4181_G_T        1524    20486   2854260
210_G_T 6402_C_T        1524    31188   2856741
210_G_T 7124_C_T        1524    34929   2862987
210_G_T 8986_C_T        1524    42856   2882122
"""

GENE_RANGES = [
    (-1, 265, ['gene=5UTR'], 0),
    (265, 805, ['gene=nsp1'], 0),
    (805, 2719, ['gene=nsp2'], 0),
    (2719, 8554, ['gene=nsp3'], 0),
    (8554, 10054, ['gene=nsp4'], 0),
    (10054, 10972, ['gene=nsp5'], 0),
    (10972, 11842, ['gene=nsp6'], 0),
    (11842, 12091, ['gene=nsp7'], 0),
    (12091, 12685, ['gene=nsp8'], 0),
    (12685, 13024, ['gene=nsp9'], 0),
    (13024, 13441, ['gene=nsp10'], 0),
    (13441, 13468, ['gene=nsp12A'], 0),
    (13467, 16236, ['gene=nsp12B'], 0),
    (13475, 13503, ['gene=orf1ab_SL1'], 0),
    (13487, 13542, ['gene=orf1ab_SL2'], 0),
    (16236, 18039, ['gene=nsp13'], 0),
    (18039, 19620, ['gene=nsp14'], 0),
    (19620, 20658, ['gene=nsp15'], 0),
    (20658, 21552, ['gene=nsp16'], 0),
    (21552, 21562, ['gene=spacer1'], 0),
    (21562, 25384, ['gene=S'], 0),
    (25384, 25392, ['gene=spacer2'], 0),

    (25392, 26220, ['gene=ORF3a'], 0),
    (26220, 26244, ['gene=spacer3'], 0),
    (26244, 26472, ['gene=E'], 0),
    (26472, 26522, ['gene=spacer4'], 0),

    (26522, 27191, ['gene=M'], 0),
    (27191, 27201, ['gene=spacer5'], 0),
    (27201, 27387, ['gene=ORF6'], 0),
    (27387, 27393, ['gene=spacer6'], 0),
    (27393, 27759, ['gene=ORF7a'], 0),
    (27755, 27887, ['gene=ORF7b'], 0),
    (27887, 27893, ['gene=spacer7'], 0),
    (27893, 28259, ['gene=ORF8'], 0),
    (28259, 28273, ['gene=spacer7'], 0),
    (28273, 29533, ['gene=Nucleoprotein'], 0),
    (29533, 29557, ['gene=spacer8'], 0),
    (29557, 29674, ['gene=ORF10'], 0),
    (29608, 29644, ['gene=ORF10_SL1'], 0),
    (29628, 29657, ['gene=ORF10_SL2'], 0),
    (29674, 29903, ['gene=3UTR'], 0),
    (29727, 29768, ['gene=STEMLOOP'], 0),
    ]

def read_and_bin(file):
    df = pd.read_table(file, header=None, nrows=10000)
    df['x'] = df[1].apply(lambda x: int(x.split('_')[0]))
    df['y'] = df[2].apply(lambda x: int(x.split('_')[0]))
    df = df.loc[df['x'] < 30000]
    df = df.loc[df['y'] < 30000]
    rows_0 = [True if re.search(r"_[A-Z]{2,}_", x) else False for x in df[0].values.tolist()]
    rows_1 = [True if re.search(r"_[A-Z]{2,}_", x) else False for x in df[1].values.tolist()]
    rows = rows_0 or rows_1
    indices = [i for i, x in enumerate(rows) if x]
    df = df.drop(index=indices)
    df.drop([0, 1, 2], axis=1, inplace=True)
    df.columns = ['z', 'x','y']
    # bin data, step = 5 nt
    values = set(df['x'].values.tolist())
    df.sort_values(by=['x', 'y'], inplace=True)
    index_list = pd.IntervalIndex([
        pd.Interval(round(i.left), round(i.right), i.closed)
        for i in pd.interval_range(0, 30000, 15000, closed='left')
    ])
    for interval in index_list:
        values = interval
        df['x'] = df['x'].apply(lambda x: values.left if values.left < x <= values.right else x)
        df['y'] = df['y'].apply(lambda x: values.left if values.left < x <= values.right else x)

    #todo - normalize data
    df[f"{file}"] = (df['z'] - min(df['z'])) / (max(df['z']) - min(df['z'])) + 1
    df.to_csv(f'/home/daria/data/processed/{file}.csv')
    return df.groupby(['x', 'y']).agg({f"{file}": 'sum'}).reset_index()

def read_and_bin_by_gene(file, indels = False, split_into_bins = False, groupby=False):
    df = pd.read_table(file, header=None)
    gene_ranges = GENE_RANGES
   # gene_ranges[0] = gene_ranges.apply(lambda x: re.search(r"product\=(.*);",x[9]).group(1) if re.search(r"product\=(.*);",x[9]) else x[7], axis = 1)
    # NC_045512.2   	0	265

    if indels == True:
        #29725_ATTT_A
        rows_0 = [True if re.search(r"_[A-Z]{2,}_", x) else False for x in df[0].values.tolist()]
        rows_1 = [True if re.search(r"_[A-Z]{2,}_", x) else False for x in df[1].values.tolist()]
        rows = rows_0 or rows_1
        df = df.loc[rows]

    df.drop([0, 1], axis=1, inplace=True)
    df.columns = ['x','y','z']
    df = df.loc[df['x'] < 30000]
    df = df.loc[df['y'] < 30000]
    if split_into_bins == True:
        index_list = pd.IntervalIndex([
            pd.Interval(round(i.left), round(i.right), i.closed)
            for i in pd.interval_range(0, 30000, 10, closed='left')
        ])

        for interval in index_list:
            values = interval
            df['x'] = df['x'].apply(lambda x: values.left if values.left < x <= values.right else x)
            df['y'] = df['y'].apply(lambda x: values.left if values.left < x <= values.right else x)

    df['c1'] = df.apply(lambda x: [], axis=1)
    df['c2'] = df.apply(lambda x: [], axis=1)
    for gene_row in gene_ranges:
        #(25392, 26220, ['gene=ORF3a'], 0)
        df['c1'] = df.apply(lambda x: np.append(x['c1'], gene_row[2]) if gene_row[0] < x['x'] <= gene_row[1] else x['c1'], axis=1)
        df['c2'] = df.apply(lambda x: np.append(x['c2'], gene_row[2]) if gene_row[0] < x['y'] <= gene_row[1] else x['c2'], axis=1)

    #df['z'] = np.log2((df['z'] - min(df['z'])) / (max(df['z']) - min(df['z'])) + 1)
    df = df.explode(['c1'])
    df = df.explode(['c2'])
    if groupby == True:
        genewize = df.groupby(['c1', 'c2']).agg({'z': 'sum'}).reset_index()
        df = df.merge(genewize, on=['c1', 'c2'], how='outer').drop(['z_x'], axis = 1)
        df.columns = ['x', 'y', 'c1', 'c2', 'z']
    return df
