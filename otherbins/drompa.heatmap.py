#! /usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('profile-K562-aroundTSS.ChIPread.tsv', sep="\t", index_col=0)
df_s = df.sort_values('0', ascending=False)

plt.figure(figsize=(6, 16))
sns.heatmap(df_s,
            vmin=0, vmax=100,
            cmap="Reds",
            yticklabels=False)

plt.savefig('heatmap.png')
