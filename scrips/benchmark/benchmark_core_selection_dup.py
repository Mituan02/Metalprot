'''
Select cores which are duplicates of the ones used in vdM library.
'''

import os
import shutil

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/benchmark/benchmark_core_selection_dup.py
'''

root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/'

workdir = root_dir + '_Seq_core_3contact_all/'

outdir = root_dir + '_Seq_core_3contact_all_dup/'

os.makedirs(outdir, exist_ok=True)

dup_file = root_dir + 'dup_summary.txt'

dup_dict = {}

with open(dup_file, 'r') as f:
    for line in f.readlines():
        infos = line.split('\t')
        dup_dict[infos[0]] = [infos[i].split('\n')[0] for i in range(1, len(infos))]

dup_sel_dict = {}

for p in os.listdir(workdir):
    key = p.split('.')[0]
    if key in dup_dict.keys():
        if len(dup_dict[key]) <= 1:
            continue
        sel = dup_dict[key][-1]
        dup_sel_dict[sel] = key
        #pdb_sel = root_dir + '_Seq_core_filter/' + sel + '.pdb'
        #shutil.copyfile(pdb_sel, outdir + p)


with open(outdir + 'dup_sel_dict.tsv', 'w') as f:
    for k in dup_sel_dict.keys():
        f.write(k + '\t' + dup_sel_dict[k] + '\n')


