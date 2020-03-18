#! /usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import os

def drawHeatmap(input, output, vmax, sortid, notsort):
    nsample = len(input)
    figwidth = 2 + 2*nsample
    fig = plt.figure(figsize=(figwidth, 16))

    if sortid >= nsample:
        print ("Error: --sortid " + str(sortid) + " is larger than the sample number.")
        sys.exit()

    df = pd.read_csv(input[sortid], sep="\t", index_col=0)
    df0 = df.loc[:,'0'].values

    for i, sample in enumerate(input):
        df = pd.read_csv(sample, sep="\t", index_col=0)
        if notsort:
            df_sort = df
        else:
            df["ref"] = df0
#            print(df)
            df_sort = df.sort_values('ref', ascending=False)
            del df_sort['ref']
#            print(df_sort)

        ax = plt.subplot(1, nsample, i+1)
        sns.heatmap(df_sort,
                    vmin=0, vmax=vmax,
                    cmap="Reds",
                    yticklabels=False)

    plt.savefig(output + '.png')

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input matrix (tsv file)", type=str, nargs='+')
    parser.add_argument("-o","--output", help="Output prefix", type=str)
    parser.add_argument("--vmax", help="Max value of heatmap", type=str, default=100)
    parser.add_argument("--sortid", help="Sample number for sorting", type=str, default=0)
    parser.add_argument("--notsort", help="do not sort the rows", action='store_true')

    args = parser.parse_args()
#    print(args)

    print ("Input files:")
    for sample in args.input:
        print ("\t" + sample)
    print ("Output file:")
    print ("\t" + args.output + ".png")

    for sample in args.input:
        if not os.path.isfile(sample):
            print ("Error: " + sample + " is not found.")
            sys.exit()

    drawHeatmap(args.input, args.output, args.vmax, int(args.sortid), args.notsort)
