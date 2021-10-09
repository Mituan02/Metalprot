import os
import prody as pr
from numpy.core.fromnumeric import argmin
from metalprot.database import database_vdMAtomOrder
metal_sel = 'ion or name NI MN ZN CO CU MG FE' 




'''
Developement test.

#test the first one
pdb = pr.parsePDB(workdir + '1992_6enl_ZN_1_mem0.pdb')

contact_ind, contact_name, contact_resind, contact_resname = get_contact_info(pdb)

pdb.select('index ' + str(contact_ind))[0].getCoords()

asp_glu_oxy_shift(pdb)

pdb.select('index ' + str(contact_ind))[0].getCoords()

'''

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210826/AAMetalPhiPsi_GLU_reps/'

pdbs = []
for file in os.listdir(workdir):
    #print(file)
    if '.pdb' not in file:
        continue
    pdb = pr.parsePDB(workdir + file)
    database_vdMAtomOrder.asp_glu_oxy_shift(pdb)
    pdbs.append(pdb)

outdir = workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210826/AAMetalPhiPsi_GLU_reps_shift/'
os.makedirs(outdir, exist_ok=True)

for pdb in pdbs:
    pr.writePDB(outdir + pdb.getTitle(), pdb)