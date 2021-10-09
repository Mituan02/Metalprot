'''
use master and qbits to evaluate the constructed helixs compared with the std coiledcoils.
'''

from smallprot import query
import subprocess
import os
import shlex
import prody as pr
import numpy as np

### query database.

pdb = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_candidates/target_c4__x1_20_x2_5_t.pdb.pdb'
targetList = '/mnt/e/GitHub_Design/Qbits/database/pds_list_2p5.txt'
rmsdCut = 0.5
master_query_top = None
outfile = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_candidates/all_stdout'

query.master_query(pdb, targetList, rmsdCut, 
    topN=master_query_top, outfile=outfile, clobber=False)


# Generate pds of the target.
fit_full_pdb_path = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_delta1_candidates/fit_full.pdb'
cmd = 'createPDS --type target --pdb {}'.format(fit_full_pdb_path)
subprocess.run(shlex.split(cmd))


### align local.

pdb = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_delta1_candidates/C4_16_3fms_NI_1_HIS_2_x_86_y_30_t.pdb.pdb'
targetList = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_delta1_candidates/pds_list_2p5.txt'
rmsdCut = 0.5
master_query_top = None
outfile = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_candidates/stdout'

query.master_query(pdb, targetList, rmsdCut, 
    topN=master_query_top, outfile=outfile, clobber=False)

### Superpose

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4/output_delta1_candidates/'

target = pr.parsePDB(workdir + 'fit_full.pdb')

mobile = pr.parsePDB(workdir + 'C4_16_3fms_NI_1_HIS_2_x_86_y_30_t.pdb.pdb')

#The numbers here need to be changed based on Master results.
aa_indices = list(range(6, 13)) + list(range(48, 55)) + list(range(90, 97)) + list(range(132, 139)) 

#pr.superpose(mobile.select('name CA'), target.select('name CA and resindex ' + ' '.join([str(a) for a in aa_indices])))

#target_sel = target.select('name CA and resindex ' + ' '.join([str(a) for a in aa_indices]))

mobile_coords = mobile.select('bb').getCoords()

target_coords = np.zeros((len(aa_indices)*4, 3), float)

for ind in range(len(aa_indices)):
    for i in range(4):
        target_coords[ind*4 + i] = target.select('bb and resindex ' + str(aa_indices[ind])).getCoords()[i]


pr.calcTransformation(mobile_coords, target_coords).apply(mobile)

pr.writePDB(workdir + 'fit_target2.pdb', mobile)

### Generate mutated one.

target_mutate = target.copy()

