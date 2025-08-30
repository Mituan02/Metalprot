'''
One assumption of the metal vdm based design is that the vdm database represent possibiities of Zn-aa in nature.
However, this is hard to prove. Here we just check the number of ZN structure generate each year to see if the number is decrease.
'''

import os

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_reps/'

date_count_dict = {}
for file in os.listdir(workdir):
    if '.pdb' not in file:
        continue
    date = file.split('_')[0]
    if date in date_count_dict.keys():
        date_count_dict[date] += 1
    else:
        date_count_dict[date] = 1

outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/reason/'

with open(outdir + 'date_vdm_count.tsv', 'w') as f:
    f.write('date\tcount\n')
    for key in date_count_dict.keys():
        f.write(key + '\t' + str(date_count_dict[key]) + '\n')