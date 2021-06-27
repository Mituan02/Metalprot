import os
import sys
import prody as pr

#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct, extract_vdm

'''
python /mnt/e/DesignData/ligands/LigandBB/Design_Sam/run_search_struct.py
'''

# Generate queryss
queryss = []

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

print(query_dir)

#Get query pdbs 
querys = extract_vdm.extract_all_centroid(query_dir + '20210421_Archive', summary_name = '_summary.txt', file_name_includes = ['cluster', '7_'], score_cut = 0, clu_num_cut = 10)

queryss.append(querys)

#Get query_2nd pdbs 
# TO DO: currently, cannot superimpose to cluster with phipsi angle, database starting with '5_'.

query_2nds = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M1-1_AAMetalSc'], score_cut = 1, clu_num_cut = 50)

queryss.append(query_2nds)

print(len(queryss[0]))
print(len(queryss[1]))

contact_querys = None

_2nd_querys = None

# run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/Design_Sam/'

outdir = workdir + 'output_3vdm_20210622/'

target_path = workdir + '3into4_helix_assembly_renum.pdb'

rmsd_cuts = [0.5, 0.5, 0.5]

dist_cuts = [1.5, 1.5, 1.5]

num_iter = 3

clash_query_query = 2.3

clash_query_target = 2.3

use_sep_aas = [False, False, False]

tolerance = 0.5

fine_dist_cut = 0.25

win_filter = None
win_filter = list(range(1,10))
win_filter.extend(list(range(47,67)))
win_filter.extend(list(range(98,107)))
win_filter.extend(list(range(108,117)))
win_filter.extend(list(range(154,169)))

ss = search_struct.Search_struct(target_path, outdir, queryss, rmsd_cuts, dist_cuts, num_iter, clash_query_query, clash_query_target, use_sep_aas, 
    tolerance, fine_dist_cut = fine_dist_cut, win_filter = win_filter, contact_querys = contact_querys, secondshell_querys=_2nd_querys)

#ss.run_search_struct()
ss.run_iter_search_structure()

ss.run_search_structure_member()