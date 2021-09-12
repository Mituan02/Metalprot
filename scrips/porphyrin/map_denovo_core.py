'''
Here I want to check the mapping the core supperimpotion of de novo protein 7jrq by Sam M. 
'''
import os
import sys
import numpy as np
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/Metalprot/scripts_other/porphyrin/')
import porphyrin_library
from metalprot import ligand_database

ligands = ['HEM', 'HNI', 'COH', 'HEB', 'FDE', 'ZNH', 'HEC', 'HEA', 'HAS', 'MNH', 'MNR']
    
ligand_sel = 'resname ' + ' '.join(ligands)

metal_sel = 'name NI MN ZN CO CU MG FE' 

#From porphyrin_library
def extract_core(pdb, extend = 4):
    '''
    Extract the porphyrin ligand first. 
    Then extract the metal of the ligand.
    '''
    metal_cores = []
    _lgd = pdb.select(ligand_sel)
    if not _lgd:
        #print('Failed no ligand: ' + pdb.getTitle())
        return

    count = 0
    for resind in np.unique(_lgd.getResindices()):
        mt = _lgd.select('resindex ' + str(resind) + ' and ' + metal_sel)
        if not mt:
            #print('Failed no metal: ' + pdb.getTitle())
            continue

        if mt[0].getName() not in metal_sel or mt[0].getResname() not in ligand_sel:
            #print(mt.getIndices())
            #print(mt.getNames())
            #print('Failed others: ' + pdb.getTitle())
            continue
        
        ni = mt[0]
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb.select('protein and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur'):
            print('Not find: ' + pdb.getTitle())
            continue 
        
        ind_check = set()
        for a_n in all_near.select('nitrogen or oxygen or sulfur'):
            ind = a_n.getResindex()
            if ind in ind_check:
                continue
            ind_check.add(ind)
            ext_inds = ligand_database.extend_res_indices([ind], pdb, extend)
            count += 1
            sel_pdb_prody = pdb.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
            metal_cores.append((pdb.getTitle() + '_' + ni.getResname() + '_'+ str(count), sel_pdb_prody)) 
    
    return metal_cores


workdir = '/mnt/e/DesignData/ligands/porphyrin/others/'

_7jrq = pr.parsePDB(workdir + '7jrq.pdb') # The 'MN1' is changed to 'MN'
ligands.append('SMU') # special porphyrin ligand in 7jrq.
ligand_sel = 'resname ' + ' '.join(ligands)
cores = extract_core(_7jrq)

ligand_database.writepdb(cores, workdir)

_7jrq_core = pr.parsePDB(workdir + '7jrq_SMU_1.pdb')

lib_dir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/M1-1_AAMetalSc_HIS_cluster05/'

core_pdb_reps = []
for path in os.listdir(lib_dir):
    if '.' in path:
        continue
    #print(lib_dir + path)
    for file in os.listdir(lib_dir + path):
        if '.pdb' not in file:
            continue
        #print(lib_dir + path + '/' + file)
        pdb = pr.parsePDB(lib_dir + path + '/' + file)
        core_pdb_reps.append((pdb, path))

#From porphyrin_library
def porphyrin_superimpose_sel(pdb, metal_sel):
    metal_ind = pdb.select(metal_sel)[0].getIndex()
    contact_atom_resind = pdb.select('protein and within 2.83 of index ' + str(metal_ind))[0].getResindex()
    contact_aa_bb_inds = pdb.select('heavy and resindex ' + str(contact_atom_resind)).getIndices()
    sel = pdb.select('index ' + str(metal_ind) + ' ' + ' '.join([str(x) for x in contact_aa_bb_inds]))
    return sel



_7jrq_sel = porphyrin_superimpose_sel(_7jrq_core, metal_sel)
min_rmsd = 5
rmsds = []
min_core = None
min_7jrq = None
for c_c in core_pdb_reps:
    c= c_c[0]
    try:
        core_sel = porphyrin_superimpose_sel(c, metal_sel)
    except:
        print(c_c)

    if len(_7jrq_sel) != len(core_sel):
        print('Failed superimpose: ' + c.getTitle())
        continue
    tr = pr.calcTransformation(_7jrq_sel, core_sel)
    tr.apply(_7jrq_sel)

    rmsd = pr.calcRMSD(_7jrq_sel, core_sel)
    rmsds.append((rmsd, c.getTitle(), c_c[1]))
    if rmsd < min_rmsd:
        min_rmsd = rmsd
        min_core = c
        tr.apply(_7jrq_core)
        min_7jrq = _7jrq_core.copy()

if min_core:
    pr.writePDB(workdir + 'rmsd_' + str(round(min_rmsd, 3)) + '_' + min_core.getTitle(), min_core)
    pr.writePDB(workdir + 'sp_'+ min_7jrq.getTitle(), min_7jrq)

with open(workdir + 'rmsd.tsv', 'w') as f:
    f.write('rmsd\tname\tcluster\n')
    for r in rmsds:
        f.write(str(r[0]) + '\t' + r[1] + '\t' + r[2] + '\n')
    