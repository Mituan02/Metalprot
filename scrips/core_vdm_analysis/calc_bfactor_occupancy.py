'''
Analyze the bfactor and occupancy value of vdms.
'''

import os
import pickle
import prody as pr
from metalprot.basic import vdmer

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211008/20211009_selfcenter/pickle_noCYS_alignBB/'


with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_vdms = pickle.load(f)

info = []
for _vdm in all_vdms:
    contact_atom = vdmer.get_contact_atom(_vdm.query)
    info.append((_vdm.query.getTitle(), contact_atom.getBeta(), contact_atom.getOccupancy()))

outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211008/reason/'

with open(outdir + 'bfactor_occupancy.tsv', 'w') as f:
    f.write('title\tbfactor\toccupancy\n')
    for o in info:
        f.write(o[0] + '\t' + str(o[1]) + '\t' + str(o[2]) + '\n')
