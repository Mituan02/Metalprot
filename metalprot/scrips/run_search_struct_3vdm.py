import os
import sys
import prody as pr

#You can either add the python package path.
sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct, extract_vdm

# Generate queryss
queryss = []


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

print(query_dir)

#Get query pdbs 

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1-1_AAMetalSc_HIS_cluster02'], score_cut = 1, clu_num_cut = 100)

print(len(querys))

queryss.append(querys)

#Get query_2nd pdbs 

query_2nds = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1-1_AAMetalSc_HIS_cluster02'], score_cut = 1, clu_num_cut = 100)

queryss.append(query_2nds)

#Get query_3rd pdbs 

query_3rds = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1-1_AAMetalSc_HIS_cluster02'], score_cut = 1, clu_num_cut = 100)

queryss.append(query_3rds)

# run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_test_3vdm2/'

target_path = workdir + '5od1_zn.pdb'

rmsd_cuts = [0.5, 0.5, 0.5]

dist_cuts = [1, 1, 1]

num_iter = 3

clash_query_query = 2.3

clash_query_target = 2.3

use_sep_aas = [False, False, False]

tolerance = 0.5

ss = search_struct.Search_struct(target_path, outdir, queryss, rmsd_cuts, dist_cuts, num_iter, clash_query_query, clash_query_target, use_sep_aas, tolerance)

#ss.run_search_struct()
ss.run_round_search_structure()