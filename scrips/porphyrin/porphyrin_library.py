import prody as pr
import numpy as np
import os
from metalprot import ligand_database

### Default selection

ligands = ['HEM', 'HNI', 'COH', 'HEB', 'FDE', 'ZNH', 'HEC', 'HEA', 'HAS', 'MNH', 'MNR']
    
ligand_sel = 'resname ' + ' '.join(ligands)

metal_sel = 'name NI MN ZN CO CU MG FE' 


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



workdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/all_pdbs/'

pdbs = []

for path in os.listdir(workdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(workdir + path)
    pdbs.append(pdb)



all_cores = []

for pdb in pdbs:
    cores = extract_core(pdb)
    if not cores:
        continue
    all_cores.extend(cores)


outdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/contact_cores/'

ligand_database.superimpose_core_and_writepdb(all_cores, all_cores[0], metal_sel, outdir)

#core_pdbs = [c[1] for c in all_cores]
core_pdbs =[]

for path in os.listdir(outdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(outdir + path)
    core_pdbs.append(pdb)


clusters = ligand_database.reduce_dup(core_pdbs, metal_sel)

outdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/contact_core_reps/'

ligand_database.extract_rep_and_writepdb(core_pdbs, clusters, metal_sel, outdir)

ligand_database.write_dup_summary(outdir, core_pdbs, clusters)


### superimpose on the 'metal->bb(N CA C)'

core_pdb_reps = []
for path in os.listdir(outdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(outdir + path)
    core_pdb_reps.append(pdb)

def porphyrin_superimpose_sel(pdb, metal_sel):
    metal_ind = pdb.select(metal_sel)[0].getIndex()
    contact_atom_resind = pdb.select('protein and within 2.83 of index ' + str(metal_ind))[0].getResindex()
    contact_aa_bb_inds = pdb.select('name N CA C and resindex ' + str(contact_atom_resind)).getIndices()
    sel = pdb.select('index ' + str(metal_ind) + ' ' + ' '.join([str(x) for x in contact_aa_bb_inds]))
    return sel

def superimpose_and_writepdb(pdbs, outdir, metal_sel):
    '''
    Superimpose on the metal->bb(N CA C)
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    first = pdbs[0]
    first_sel = porphyrin_superimpose_sel(first, metal_sel)
    for pdb in pdbs:
        pdb_sel = porphyrin_superimpose_sel(pdb, metal_sel)

        if len(first_sel) != len(pdb_sel):
            print('Failed superimpose: ' + pdb.getTitle())
            continue

        pr.calcTransformation(pdb_sel, first_sel).apply(pdb_sel)

        pr.writePDB(outdir + pdb.getTitle(), pdb)


outdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/contact_core_reps_aligned/'

superimpose_and_writepdb(core_pdb_reps, outdir, metal_sel)


###########################################################################


def extract_water_core(pdb, extend = 4):
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
        all_water = pdb.select('water and within 2.83 of index ' + str(ni_index))
        if not all_water or not all_water.select('oxygen'):
            print('Not find: ' + pdb.getTitle())
            continue 
        
        ind_check = set()
        for a_n in all_water.select('oxygen'):
            water_ind = a_n.getIndex()

            all_near_water = pdb.select('protein and within 3.4 of index ' + str(water_ind))
            if not all_near_water or not all_near_water.select('nitrogen or oxygen or sulfur'):
                continue
            for w_a_n in all_near_water.select('nitrogen or oxygen or sulfur'):
                ind = w_a_n.getResindex()
                if ind in ind_check:
                    continue
                ind_check.add(ind)
                ext_inds = ligand_database.extend_res_indices([ind], pdb, extend)
                count += 1
                sel_pdb_prody = pdb.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()) + ' '+ str(a_n.getResindex()) )
                metal_cores.append((pdb.getTitle() + '_' + ni.getResname() + '_'+ str(count), sel_pdb_prody)) 
    
    return metal_cores




workdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/all_pdbs/'

pdbs = []

for path in os.listdir(workdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(workdir + path)
    pdbs.append(pdb)



all_water_cores = []

for pdb in pdbs:
    cores = extract_water_core(pdb)
    if not cores:
        continue
    all_water_cores.extend(cores)


outdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/contact_wather_cores/'

ligand_database.superimpose_core_and_writepdb(all_water_cores, all_water_cores[0], metal_sel, outdir)

#core_pdbs = [c[1] for c in all_cores]
core_water_pdbs =[]

for path in os.listdir(outdir):
    if '.pdb' not in path:
        continue
    pdb = pr.parsePDB(outdir + path)
    core_water_pdbs.append(pdb)


water_clusters = ligand_database.reduce_dup(core_water_pdbs, metal_sel)

outdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/contact_water_core_reps/'

ligand_database.extract_rep_and_writepdb(core_water_pdbs, water_clusters, metal_sel, outdir)

ligand_database.write_dup_summary(outdir, core_water_pdbs, water_clusters)