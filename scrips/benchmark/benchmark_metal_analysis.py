'''
The function here is for analyzing the effect of metal database.
By searching Fe 3 contact core with Zn database. We first want to know how much could we get. 

'''

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/benchmark/benchmark_metal_analysis.py
'''

import os
import sys
import re
import pandas as pd
import numpy as np

workdir = '/mnt/e/DesignData/ligands/FE_rcsb/20211206/_Seq_core_3contact/out_eval_220111/'

outfile = '/mnt/e/DesignData/ligands/FE_rcsb/20211206/_Seq_core_3contact/out_eval_220111/_all_log_info.tsv'

all_files = []
Not_found = []

df_all = pd.DataFrame()


for dir in os.listdir(workdir):
    
    if not os.path.isdir(workdir + dir):

        continue
    for v in os.listdir(workdir + dir):
        #file_path = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/output_eval_selfcenter__20211118-223136/_best_summary_45_5od1_zn_20211118-223136.tsv'
        #v = '_best_summary_45_5od1_zn_20211118-223136.tsv'

        file_path = workdir + dir + '/' + v
        
        if ('_best') not in v or ('.tsv' not in v):
            continue
        title = '_'.join(v.split('_')[3:6])
        all_files.append(title)
        
        df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)

        if df.shape[0] < 1:
            Not_found.append(title)
            continue

        print(df.shape)

        nature_ind = 0

        df_ext = {}
        df_ext["title"] = title
        df_ext["AAs"] = df["AAs"].iloc[nature_ind]
        df_ext["CluScore"] = df["CluScore"].iloc[nature_ind]
        df_ext["OverlapScoreLn"] = df["OverlapScoreLn"].iloc[nature_ind]
        df_ext["OverlapScore"] = df["OverlapScore"].iloc[nature_ind]
        df_ext["Nature_GeoRmsd"] = df["GeoRmsd"].iloc[nature_ind]
        df_ext["eval_min_rmsd"] = df["eval_min_rmsd"].iloc[nature_ind]

        df_all = df_all.append(df_ext, ignore_index=True)


df_all.to_csv(outfile, sep='\t')


'''
#Filter out the ones are not finished for all the metal_metal_dists.

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason_CYS/mm_dist/'
file_path = workdir + '_all_log_info.tsv'

df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)



for t in np.unique(df['title']):
    if sum(df['title'] == t) < 8:
        df = df[df['title'] != t]

outfile = workdir + '_all_log_info_filtered8.tsv'

df.to_csv(outfile, sep='\t')

'''