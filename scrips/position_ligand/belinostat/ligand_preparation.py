'''
Belinostat only exist in one pdb structure 5een.
The ligand structure is not ideal regarding the bind with ZN and the double bond. 
Note that the group 'O=C-N-O' binding ZN should be able to be extracted from library to be learned.
MOE is used to stimulate diff structures of the lig.
Finally we add a metal at the ideal position for the searching purpose. 
'''

import os
import prody as pr
import numpy as np
from metalprot.basic import prody_ext

workdir = '/mnt/e/DesignData/Metalloenzyme/ligs/'

moe_indir = workdir + 'meo_50g_amber14eht_md/'
meo_outdir = workdir + 'meo_50g_amber14eht_md_out/'
os.makedirs(meo_outdir, exist_ok=True)

lig = pr.parsePDB(workdir + '50gZN_fix.pdb')

moe_ligs = [pr.parsePDB(moe_indir + file) for file in os.listdir(moe_indir) if '.pdb' in file]

for i in range(len(moe_ligs)):
    moelig = moe_ligs[i]
    mob_coord = [moelig.select('name ' + x)[0].getCoords() for x in ['O4', 'N2', 'C12', 'O3']]
    target_coord = [lig.select('name ' + x)[0].getCoords() for x in ['O4', 'N2', 'C12', 'O3']]
    pr.calcTransformation(np.array(mob_coord), np.array(target_coord)).apply(moelig)
    ag = prody_ext.combine_ags([moelig, lig.select('name ZN').toAtomGroup()], 'x', ['A', 'A'])
    pr.writePDB(meo_outdir + '50g_md_' + str(i), ag)

