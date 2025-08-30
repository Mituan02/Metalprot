'''
This is the command I use for the reduce program, where BUILD adds hydrogens, 
FLIP rotates and flips asn gln his hydrogens, 
(I think having both BUILD and FLIP are redundant but I don't have evidence), 
and Quiet suppresses verbose output (optional):
reduce -BUILD -FLIP -Quiet oldpdb.pdb > newpdb.PDB

I just realized the gpu computer doesn’t have reduce anymore. 
Or if it does, I don't know where the path is. 
However, it’s in the wynton amber installation: /wynton/home/grabe/shared/amber/amber18/bin/reduce
I don’t know if this is necessary, but you might have to source amber first: source /wynton/home/grabe/shared/amber/amber18/amber.sh

To install molprobity, check the website here: https://github.com/rlabduke/MolProbity.
'''

import os 

# Define pre-hydrogen and post-hydrogen database paths
db_dir = '/wynton/home/degradolab/lonelu/DesignData/_Seq_core_2ndshell_date_reps/'
db_H_dir = '/wynton/home/degradolab/lonelu/DesignData/_Seq_core_2ndshell_date_reps_withH/'
os.makedirs(db_H_dir, exist_ok = True)
for pdb in os.listdir(db_dir):
    if not pdb.endswith('.pdb'):
        continue
    os.system('/wynton/home/degradolab/lonelu/software/molprobity/molprobity/bin/linux/reduce -BUILD -FLIP -Quiet {} > {}'.format(
        db_dir + pdb, db_H_dir+pdb))

'''
#>>> Run local

import os
workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/68-77-79_DHH/'
pdb = '1ukr_68-77-79_DHH_lig20.pdb'
os.system('~/molprobity/bin/linux/reduce -BUILD -FLIP -Quiet {} > {}'.format(workdir + pdb, workdir + 'a_'+pdb))
'''
