'''
Is the 2ndshell close to the 1st shell regarding the sequence position?
The COMBS database remove close aa contact.
Searching the 2ndshell may make more sense to use a new database contain close aa contact.
The difference between the close aa contact specified by metal is unknow.

To answer the question, please check '../Metalprot/scrips/database_generate/database_prep_2ndshell_probe.py'
'''
import os
import pandas as pd

#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_seq_core_2ndshell_date_reps_probe2ndshell/'
workdir = '/mnt/e/DesignData/ligands/all/20220710/_seq_core_2ndshell_reps_probe2ndshell/'

aas = ['ASP', 'GLU', 'HIS', 'CYS']

infos = []
for aa in aas:

    for file in os.listdir(workdir + aa + '/'):
        if not '.pdb' in file:
            continue

        info = []
        xs = file.split('.')[0].split('_')
        #>>> Is the year in the title
        if len(xs) == 16: 
            title = xs[0] + '_' + xs[1] + '_' + xs[2] + '_' + xs[3]
            info.append(title)
            info.extend(xs[4:])
            info.append(xs[6] == xs[11])
            info.append(int(xs[7]) - int(xs[12]))
        else:
            title = xs[0] + '_' + xs[1] + '_' + xs[2]
            info.append(title)
            info.extend(xs[3:])
            info.append(xs[5] == xs[10])
            info.append(int(xs[6]) - int(xs[11]))
        infos.append(info)
        
cols = [ 'title', 'hb_type', 'hb', '1st_chid', '1st_resnum', '1st_aa', '1st_atom', '1st_atomtype',
        '2nd_chid', '2nd_resnum', '2nd_aa', '2nd_atom', '2nd_atomtype', 'SameChid', 'seq_dist']
        
pd_infos = pd.DataFrame(infos, columns = cols )

pd_infos.to_csv(workdir + '_2ndshell_seq_dist_summary.tsv', sep='\t')