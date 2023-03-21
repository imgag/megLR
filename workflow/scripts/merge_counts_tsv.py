#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
from os import path
import pandas as pd
import numpy as np
from functools import reduce

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Merge tab separated files on a given field using pandas.""")
parser.add_argument(
    '-j', metavar='join', type=str, help="Join type (outer).", default="outer")
parser.add_argument(
    '-f', metavar='field', type=str, help="Join on this field (Reference).", default="Reference")
parser.add_argument(
    '-o', metavar='out_tsv', type=str, help="Output tsv (merge_tsvs.tsv).", default="merge_tsvs.tsv")
parser.add_argument(
    '-z', action="store_true", help="Fill NA values with zero.", default=False)
parser.add_argument(
    'tsvs', metavar='input_tsvs', nargs='*', type=str, help="Input tab separated files.")


if __name__ == '__main__':
    args = parser.parse_args()

    dfs = {x: pd.read_csv(x, sep="\t") for x in args.tsvs}

    ndfs = []
    for x, df in dfs.items():
        # Transform counts to integers:
        df = df.rename(columns={'NumReads': 'Count', 'Name': 'Reference'})
        df.Count = np.array(df.Count, dtype=int)
        # Take only non-zero counts:
        df = df[df.Count > 0]
        df = df[["Reference", "Count"]]
        df = df.sort_values(by=["Count"], ascending=False)
        name_search = re.search('(?<=isoform_counts\/)(.*)(?=_counts)', path.dirname(x))       
        if name_search:
            name =  name_search.group(1)
        df = df.rename(columns={'Count': name})
        ndfs.append(df)
    dfs = ndfs
    print(dfs)

    df_merged = reduce(lambda left, right: pd.merge(left, right, on=args.f, how=args.j), dfs)
    if args.z:
        df_merged = df_merged.fillna(0)

    df_merged.to_csv(args.o, sep="\t", index=False)
