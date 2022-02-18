from genericpath import isdir, isfile
import os
import pandas as pd

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'


df_all = []
for folder in os.listdir(workdir):
    if not os.path.isdir(folder):
        continue
    for file in  os.listdir(workdir + folder):
        if '_summary.tsv' not in file:
            continue
        if not os.path.isfile(workdir + folder + '/' + file):
            continue

        df = pd.read_csv(workdir + folder + '/' + file, sep = '\t')
        df_all.append(df)

result = pd.concat(df_all)

df_group = result.groupby(['file', 'lig'])

df_score = df_group[['score']].sum()

df_score.to_csv(workdir + '_sum_score.tsv', sep = '\t')

scores = []
for g_name, g in df_group:
    df_gg_group = g.groupby(['cg','aa','pos'])
    score = 0
    for gg_name, gg in df_gg_group:
        score += gg[['score']].max()
    scores.append((g_name, score.sum()))

with open(workdir + '_sum_score_rmdu.tsv', 'w') as f:
    f.write('lig\tscore\n')
    for s in scores:
        f.write(str(s[0]) + '\t' +  str(s[1]) + '\n')