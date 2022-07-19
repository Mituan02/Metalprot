'''
Some of the pdbs have no date in their title and if we want to add deposit date, we can use the following script.

'''

import os
import pypdb

workdir = '/mnt/e/DesignData/ligands/all/20220710/_Seq_core_2ndshell_CO/'

for pdb in os.listdir(workdir):
    if not '.pdb' in pdb:
        continue
    id = pdb[0:4]
    year = pypdb.get_info(id)['rcsb_accession_info']['deposit_date'][0:4]
    os.rename(workdir + pdb, workdir + year + '_' + pdb)


#>>> If the ZN pdb already has date.
workdir = '/mnt/e/DesignData/ligands/all/20220710/_Seq_core_2ndshell_reps_probe2ndshell/CYS/'

for pdb in os.listdir(workdir):
    if not '.pdb' in pdb:
        continue
    if '_ZN_' in pdb:
        continue
    id = pdb[0:4]
    year = pypdb.get_info(id)['rcsb_accession_info']['deposit_date'][0:4]
    os.rename(workdir + pdb, workdir + year + '_' + pdb)