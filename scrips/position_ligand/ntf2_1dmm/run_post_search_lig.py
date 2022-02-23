from genericpath import isdir, isfile
import os
import pandas as pd

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_post_search_lig.py
'''

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/wynton_20220218_3/'


df_all = []
for folder in os.listdir(workdir):
    if not os.path.isdir(workdir + folder):
        continue
    print(folder)
    for file in  os.listdir(workdir + folder):
        if '_summary.tsv' not in file:
            continue
        if not os.path.isfile(workdir + folder + '/' + file):
            continue
        print(file)
        df = pd.read_csv(workdir + folder + '/' + file, sep = '\t')
        df_all.append(df)

result = pd.concat(df_all)

result.to_csv(workdir + '_summary_vdms.tsv', sep = '\t')

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
        f.write('\t'.join([str(x) for x in s[0]]) + '\t' +  str(s[1]) + '\n')