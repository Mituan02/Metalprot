import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import search_struct, extract_vdm, ligand_database
import multiprocessing as mp

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_search_all_structs_3vdm.py
'''

# Generate queryss
queryss = []


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/'

print(query_dir)

#Get query pdbs 

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1_AAMetal'], score_cut = 0, clu_num_cut = 50)

print(len(querys))

queryss.append(querys)

#Get query_2nd pdbs 

query_2nds = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1_AAMetal'], score_cut = 0, clu_num_cut = 50)

queryss.append(query_2nds)

#Get query_3rd pdbs 

query_3rds = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M1_AAMetal'], score_cut = 0, clu_num_cut = 50)

queryss.append(query_3rds)

contact_querys = None
#contact_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M8_AtomContact4_clusters'], score_cut = 0, clu_num_cut = 2)

_2nd_querys = None
#_2nd_querys = extract_vdm.extract_all_centroid(query_dir + '20210608/', summary_name = '_summary.txt', file_name_includes = ['M7_AA2sMetal-HIS_clusters'], score_cut = 0, clu_num_cut = 0)

#print(len(_2nd_querys))

query_bivalences = extract_vdm.extract_all_centroid('/mnt/e/DesignData/ligands/ZN_rcsb/20210527', summary_name = '_summary.txt', file_name_includes = ['M5_2aa_sep_cores_bb_clusters'], score_cut = 1, clu_num_cut = 10)
print(len(query_bivalences))

# run Search_struct
def run_search(workdir, target_file, queryss, contact_querys, _2nd_querys):

    rmsd_cuts = [0.5, 0.5, 0.5]

    dist_cuts = [1.5, 1.5, 1.5]

    num_iter = 3

    clash_query_query = 2.3

    clash_query_target = 2.3

    use_sep_aas = [False, False, False]

    tolerance = 0.5

    fine_dist_cut = 0.25

    win_filter = None

    validateOriginStruct = True

    print('Working on ' + target_file)
    outdir = workdir + 'output_' + target_file.split('.')[0] + '/'
    target_path = workdir + target_file

    #target_path = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_test_full_pdbs/5mr8.pdb'
    target = pr.parsePDB(target_path)

    metal_cores = ligand_database.get_metal_core_seq(target, metal_sel = 'name ZN', extend = 0)
    win_extract = []
    for c in metal_cores:
        win_extract.extend(c[1].select('name CA').getResindices())
    win_extract.sort()

    try:
        win_filters = search_struct.extract_all_win_filter_by_bivalence(query_bivalences, target, tolerance=0.75, rmsd_cut=0.75)
        win_filter = [w for w in win_filters]
    except:
        win_filters = None 

    ss = search_struct.Search_struct(target_path, outdir, queryss, rmsd_cuts, dist_cuts, num_iter, clash_query_query, clash_query_target, use_sep_aas, 
        tolerance, fine_dist_cut = fine_dist_cut, win_filter = win_filter, contact_querys = contact_querys, secondshell_querys=_2nd_querys, validateOriginStruct = validateOriginStruct)

    try:
        #ss.run_search_struct()
        ss.run_iter_search_structure()
        ##ss.run_win_based_search()
        #ss.run_search_structure_member()

        #ss.run_search_2nshells(outpath = '/mem_combs/', rmsd=0.5)
        ### If only search 2nshell for a specific comb.
        #ss.search_2ndshell(4)
        #ss.write_2ndshell(4)
    except:
        return (target_file, False, win_filter, win_extract, None, False, False, False)

    win_search = set()
    for c in ss.combs:
        for q in c.querys:
            for w in q.win:
                win_search.add(w)
    win_search = [w for w in win_search] 
    win_search.sort()

    extract_in_search = [True if e in win_search else False for e in win_extract]

    search_in_extract = [True if e in win_extract else False for e in win_search]
    
    return (target_file, True, win_filter, win_extract, win_search, True, all(extract_in_search), all(search_in_extract))


num_cores = int(mp.cpu_count())
pool = mp.Pool(num_cores)

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_test_full_pdbs_sub2/'

target_files = []
for target_file in os.listdir(workdir):
    if target_file.endswith('.pdb'):
        target_files.append(target_file)

results = [pool.apply_async(run_search, args=(workdir, target_file, queryss, None, None)) for target_file in target_files]
results = [p.get() for p in results]

with open(workdir + '_summary.txt', 'w') as f:
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

