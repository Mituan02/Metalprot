'''
The function here is for analyzing the effect of metal-metal-distance on scoring and protein backbone sensitivity. 

After search the ~500 benchmark protein with metal-metal-distance [0.15, 0.25...1.05]
Check the score changes.

'''

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/benchmark/benchmark_mmdist_analysis.py
'''

import os
import sys
import re
import pandas as pd
import numpy as np

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_3contact_CYS/20211129_eval_mmdist/'

outfile = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_3contact_CYS/20211129_eval_mmdist/_all_log_info.tsv'


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
        title = '_'.join(v.split('_')[4:8])
        dist = int(v.split('_')[3])
        
        df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)

        if df.shape[0] < 1:
            continue

        print(df.shape)
        aa_types = [df["ClusterIDs"][0].split('_')[i][0] for i in range(3)]

        nature_ind = df.index[df['eval_is_origin']][0]
        max_clu_idx = df["CluScore"].idxmax()
        max_over_idx = df["OverlapScoreLn"].idxmax()

        df_ext = {}
        df_ext["title"] = title
        df_ext["mmdist"] = dist
        df_ext["aa_types"] = '-'.join(aa_types)
        df_ext["Nature_CluScore"] = df["CluScore"].iloc[nature_ind]
        df_ext["Nature_OverlapScoreLn"] = df["OverlapScoreLn"].iloc[nature_ind]
        df_ext["Nature_OverlapScore"] = df["OverlapScore"].iloc[nature_ind]
        df_ext["Nature_GeoRmsd"] = df["GeoRmsd"].iloc[nature_ind]
        df_ext["nature_eval_contain_origin"] = df["eval_contain_origin"].iloc[nature_ind]

        df_ext["MaxClu_CluScore"] = df["CluScore"].iloc[max_clu_idx]
        df_ext["MaxClu_OverlapScoreLn"] = df["OverlapScoreLn"].iloc[max_clu_idx]
        df_ext["MaxClu_OverlapScore"] = df["OverlapScore"].iloc[max_clu_idx]
        df_ext["MaxClu_GeoRmsd"] = df["GeoRmsd"].iloc[max_clu_idx]
        df_ext["MaxClu_eval_contain_origin"] = df["eval_contain_origin"].iloc[max_clu_idx]

        df_ext["MaxOver_CluScore"] = df["CluScore"].iloc[max_over_idx]
        df_ext["MaxOver_OverlapScoreLn"] = df["OverlapScoreLn"].iloc[max_over_idx]
        df_ext["MaxOver_OverlapScore"] = df["OverlapScore"].iloc[max_over_idx]
        df_ext["MaxOver_GeoRmsd"] = df["GeoRmsd"].iloc[max_over_idx]
        df_ext["MaxOver_eval_contain_origin"] = df["eval_contain_origin"].iloc[max_over_idx]

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