import os
import sys
import prody as pr

#You can either add the python package path.
sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct

# Generate queryss
queryss = []


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

print(query_dir)

#Get query pdbs 

subfolders_with_paths = [f.path for f in os.scandir(query_dir) if f.is_dir() if f.name[0] == '7' and '_' in f.name and 'cluster' in f.name]

querys = []

for subfolder in subfolders_with_paths:
    qs = search_struct.extract_query(subfolder + '/', file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10)
    querys.extend(qs)

queryss.append(querys)

#Get query_2nd pdbs 
# TO DO: currently, cannot superimpose to cluster with phipsi angle.
subfolders_with_paths2 = [f.path for f in os.scandir(query_dir) if f.is_dir() if f.name[0] == '4' and '_' in f.name and 'cluster' in f.name]

query_2nds = []

for subfolder in subfolders_with_paths2:
    qs = search_struct.extract_query(subfolder + '/', file_path = '/_summary.txt', score_cut = 1, clu_num_cut = 50)
    query_2nds.extend(qs)
queryss.append(query_2nds)

#Get query_3rd pdbs 
# TO DO: currently, cannot superimpose to cluster with phipsi angle.
subfolders_with_paths3 = [f.path for f in os.scandir(query_dir) if f.is_dir() if f.name[0] == '4' and '_' in f.name and 'cluster' in f.name]

query_3rds = []

for subfolder in subfolders_with_paths3:
    qs = search_struct.extract_query(subfolder + '/', file_path = '/_summary.txt', score_cut = 1, clu_num_cut = 50)
    query_3rds.extend(qs)
queryss.append(query_3rds)

# run Search_struct

workdir = '/mnt/e/DesignData/ligands/Design_Sam/'

outdir = workdir + 'output_test_3vdm/'

target_path = workdir + '3into4_helix_assembly_renum.pdb'

ss = search_struct.Search_struct(target_path, outdir, queryss, [0.5, 0.5, 0.5], [1, 1, 1], 3)

ss.run_search_struct()