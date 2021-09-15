'''
Categery by the combinations of different type of vdMs.

Specially, we want to kown how many x-x-CYS or x-CYS-CYS core (x != CYS) are contained in the database.
'''

import os
import sys
import prody as pr
from metalprot.database import database_extract as ldb


workdir = "/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/"

cores = ldb.load_cores(workdir + '_Seq_core_date_reps/')

res_num_dict = dict()

core_res_dict = dict()

for core in cores:
    try:
        resnames = sorted(core.full_pdb.select('name CA and resindex ' + ' '.join([str(r) for r in core.contact_aa_resinds])).getResnames())
        name = '-'.join(resnames)

        core_res_dict[core.full_pdb.getTitle()] = name
        if name in res_num_dict.keys():
            res_num_dict[name] += 1
        else:
            res_num_dict[name] = 1

    except:
        print(core.full_pdb.getTitle())

select_aa = ['HIS', 'CYS', 'GLU', 'ASP']

with open(workdir + 'core_aa_summary.tsv', 'w') as f:
    f.write('key\tlength\tcontain_otheraa\tfrequency\n')
    for key in res_num_dict.keys():
        names = key.split('-')
        contain_otheraa = False
        if any([True for name in names if name not in select_aa]):
            contain_otheraa = True
        f.write(key + '\t' + str(len(names)) + '\t' + str(contain_otheraa) + '\t'+ str(res_dict[key]) + '\n')


with open(workdir + 'core_aa_info.tsv', 'w') as f:
    for key in core_res_dict.keys():
        f.write(key + '\t' + core_res_dict[key] + '\n')

