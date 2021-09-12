import prody as pr
import numpy as np
import os

### Default selection

ligands = ['HEM', 'HNI', 'COH', 'HEB', 'FDE', 'ZNH', 'HEC', 'HEA', 'HAS', 'MNH', 'MNR']
    
ligand_sel = 'resname ' + ' '.join(ligands)

metal_sel = 'name NI MN ZN CO CU MG FE' 

### 

workdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/all_pdbs/'


### test case

pdb = pr.parsePDB(workdir + '4ens.pdb')

_lgd = pdb.select(ligand_sel)

#_lgd = pdb.select('not protein and not water')

#repr(_lgd)
_lgd.getResnames()

_lgd.getNames()

_lgd.getResindices()

for resind in np.unique(_lgd.getResindices()):
    mt = _lgd.select('resindex ' + str(resind) + ' and ' + metal_sel)
    print(mt.getIndices())
    print(mt.getNames())


### test all

pdbs = []

for path in os.listdir(workdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(workdir + path)
    pdbs.append(pdb)
    

def check_core(pdb):
    _lgd = pdb.select(ligand_sel)
    if not _lgd:
        print('Failed no ligand: ' + pdb.getTitle())
        return
    for resind in np.unique(_lgd.getResindices()):
        mt = _lgd.select('resindex ' + str(resind) + ' and ' + metal_sel)
        if not mt:
            print('Failed no metal: ' + pdb.getTitle())
            continue

        if len(mt) >= 2:
            print('Metal nums wrong: ' + pdb.getTitle())
            continue
        
        if mt[0].getName() not in metal_sel or mt[0].getResname() not in ligand_sel:
            print(mt.getIndices())
            print(mt.getNames())
            print('Failed others: ' + pdb.getTitle())
            continue
    return

for pdb in pdbs:
    check_core(pdb)

'''
#Failed:
Failed no metal: 1wej  #The hem contains no metal 
Failed no metal: 3fvb 
Failed no ligand: 4ens  #The name is not 'HEM' but somthing like 'CBBHEM'
Failed no ligand: 4ent
Failed no ligand: 4env
Failed no ligand: 4xmd
Failed no ligand: 7bwh
'''

