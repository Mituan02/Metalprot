import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct, extract_vdm, ligand_database
import multiprocessing as mp

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_search_all_structs_3vdm_post_analysis.py
'''

query_bivalences = extract_vdm.extract_all_centroid('/mnt/e/DesignData/ligands/ZN_rcsb/20210527', summary_name = '_summary.txt', file_name_includes = ['M5_2aa_sep_cores_bb_clusters'], score_cut = 1, clu_num_cut = 10)
print(len(query_bivalences))

###------------------------------------------------------------
# Post analysis.

def run_post_analysis(workdir, target_file):

    print('Working on ' + target_file)
    target_path = workdir + target_file

    #target_path = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_test_full_pdbs/5mr8.pdb'
    target = pr.parsePDB(target_path)

    metal_cores = ligand_database.get_metal_core_seq(target, metal_sel = 'name ZN', extend = 0)
    win_extract = []
    for c in metal_cores:
        win_extract.extend(c[1].select('name CA').getResindices())
    win_extract.sort()


    win_filter = search_struct.extract_all_win_filter_by_bivalence(query_bivalences, target, tolerance=0.75, rmsd_cut=0.75)
    win_filter = [w for w in win_filter]
    win_filter.sort()

    #win_filter = None

    win_search = set()
    out_comb_dir = workdir + 'output_' + target_file.split('.')[0] + '/combs'
    if os.path.exists(out_comb_dir):
        for c in os.listdir(out_comb_dir):
            if not c.endswith('.pdb'):
                continue
            x = int(c[:-4].split('_')[-1])
            win_search.add(x)
    win_search = [w for w in win_search] 
    win_search.sort()

    extract_in_search = [True if e in win_search else False for e in win_extract]

    search_in_extract = [True if e in win_extract else False for e in win_search]
    

    return (target_file, os.path.exists(out_comb_dir), win_filter, win_extract, win_search, len(win_search) > 0, all(extract_in_search), all(search_in_extract))



num_cores = int(mp.cpu_count())
pool = mp.Pool(num_cores)

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_test_full_pdbs_sub2/'

target_files = []
for target_file in os.listdir(workdir):
    if target_file.endswith('.pdb'):
        target_files.append(target_file)

results = [pool.apply_async(run_post_analysis, args=(workdir, target_file)) for target_file in target_files]
results = [p.get() for p in results]

with open(workdir + '_summary_all.txt', 'w') as f:
    f.write('target_file\tRun\twin_filter\twin_extract\twin_search\tFind\tFind_all_extract\tNo_extra_Find\n')
    for r in results:
        try:
            f.write(r[0] + '\t')
            f.write(str(r[1]) + '\t')
            f.write(('' if not r[2] else ';'.join([str(x) for x in r[2]])) + '\t')
            f.write(';'.join([str(x) for x in r[3]]) + '\t')
            f.write(';'.join([str(x) for x in r[4]]) + '\t')
            f.write(str(r[5]) + '\t')
            f.write(str(r[6]) + '\t')
            f.write(str(r[7]) + '\t\n')
        except:
            f.write(r[0] + '\t\n')