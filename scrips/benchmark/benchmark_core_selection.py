'''
select 3 contact binding core for benchmark.

To run this, please check ..\core_vdm_analysis/core_aa_category.py 
Based 'E:\DesignData\ligands\ZN_rcsb_datesplit\20210624\reason\core_aa_info.tsv'
In folder 'E:\DesignData\ligands\ZN_rcsb_datesplit\20210624\_Seq_core_date_reps\'
Extract the 3-contact binding (HIS, ASP, GLU) core pdbs from and copy it into a new folder.

'''

import os
import shutil

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/benchmark/benchmark_core_selection.py
'''

root_dir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/'

workdir = root_dir + '_Seq_core_date_reps/'

outdir = root_dir + '_Seq_core_3contact_CYS/'

os.makedirs(outdir, exist_ok=True)

core_info_file = root_dir + 'reason/core_aa_info.tsv'

pdb_dict = {}

with open(core_info_file, 'r') as f:
    for line in f.readlines():
        info = line.split('\t')
        if len(info[1].split('-')) != 3:
            continue
        cs = info[1].split('\n')[0]
        sati = True
        containCYS = False
        for c in cs.split('-'):
            if not c in ['HIS', 'ASP', 'GLU', 'CYS']:
                sati = False
            if c == 'CYS':
                containCYS = True
        if sati and containCYS:      
            pdb_dict[info[0]] = cs

for p in os.listdir(workdir):
    if p.split('.')[0] in pdb_dict.keys():
        shutil.copyfile(workdir + p, outdir + p)





