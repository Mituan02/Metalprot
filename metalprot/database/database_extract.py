import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
from . import core

# Basic function. 

def get_all_pbd_prody(workdir):
    '''
    find all .pdb file in a folder and load them with prody.
    return a list of pdb_prody.
    '''
    pdbs = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            pdbs.append(pdb_prody)
        except:
            print('not sure')   
    return pdbs

def writepdb(cores, outdir):
    '''
    cores = list of (pdb name, pdb_prody)
    '''

    if not os.path.exists(outdir):
            os.mkdir(outdir)
    for c in cores:
        outfile = c[0].split('.')[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

def load_cores(workdir):
    '''
    load pdb and create a Core class for vdM generating.
    '''
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            _core = core.Core(pdb_prody)
            cores.append(_core)
        except:
            print('not sure')   
    return cores


# Prepare rcsb database. // extract seq within +-3 aa for each contact aa. 



def get_metal_core_seq(pdb_prody, metal_sel, extend = 4):
    metal = metal_sel.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb_prody.select('resname HIS GLU ASP CYS and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur') or len(all_near.select('nitrogen or oxygen or sulfur')) < 3:
            continue          
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        # all_near_res = pdb_prody.select('protein and resindex ' + ' '.join([str(ind) for ind in inds]))
        # if not all_near_res or len(np.unique(all_near_res.getResindices())) < 2:
        #     continue     
        # inds_near_res =  all_near_res.getResindices()
        # ext_inds = extend_res_indices(inds_near_res, pdb_prody, extend)
        ext_inds = core.extend_res_indices(inds, pdb_prody, extend)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + metal + '_'+ str(count), sel_pdb_prody))        
    return metal_cores



def get_metal_core_seq_2ndshell(pdb_prody, metal_sel, extend = 4):
    metal = metal_sel.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb_prody.select('protein and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur') or len(all_near.select('nitrogen or oxygen or sulfur')) < 3:
            continue          
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        ext_inds = core.extend_res_indices(inds, pdb_prody, extend)
        _2ndshell_inds = core.get_2ndshell_indices(inds, pdb_prody, ni_index)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ ' '.join([str(ind) for ind in _2ndshell_inds]) + ' ' + str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + metal + '_'+ str(count), sel_pdb_prody))        
    return metal_cores

def extract_all_core_seq_from_path(workdir, metal_sel, extend = 4, extract_2ndshell = False):
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            if extract_2ndshell:
                core = get_metal_core_seq_2ndshell(pdb_prody, metal_sel, extend)
            else:
                core = get_metal_core_seq(pdb_prody, metal_sel, extend)
            if core:
                cores.extend(core)
            pdb_prody = None
        except:
            print('not sure')   
    return cores

def extract_all_core_seq(pdbs, metal_sel, extend = 4):
    cores = []
    for pdb in pdbs:
        core = get_metal_core_seq(pdb, metal_sel, extend)
        if core:
            cores.extend(core)
    return cores

def superimpose_core_and_writepdb(cores, first, metal_sel, outdir):
    '''
    Just supperimpose on the metal. And write pdb.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for c in cores:
        a_coords = first[1].select(metal_sel)[0].getCoords().reshape(1,3)
        b_coords = c[1].select(metal_sel)[0].getCoords().reshape(1,3)

        pr.calcTransformation(b_coords, a_coords).apply(c[1])

        outfile = c[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

# Prepare rcsb database. // Reduce duplication.

def reduce_dup(pdbs, metal_sel):
    '''
    Cluster pdbs based on the structure similarity.
    '''
    clusters = []
    clustered = set()
    for i in range(len(pdbs)):
        if i in clustered: continue

        cluster = [i]
        clustered.add(i)
        for j in range(i + 1, len(pdbs)):
            if j in clustered: continue
            if len(pdbs[i].select('name N CA C O')) != len(pdbs[j].select('name N CA C O')):
                continue
            #print(str(i) + '---' + str(j))
            sequence_equal = False
            if pdbs[i].select('name CA').getSequence() == pdbs[i].select('name CA').getSequence():
                sequence_equal = True
            #TO DO: The calcTransformation here will change the position of pdb. 
            #This will make the output pdb not align well. Current solved by re align.
            pr.calcTransformation(pdbs[j].select('name N CA C O'), pdbs[i].select('name N CA C O')).apply(pdbs[j])
            rmsd = pr.calcRMSD(pdbs[i].select('name N CA C O'), pdbs[j].select('name N CA C O'))

            if (sequence_equal and rmsd < 0.8) or rmsd < 0.25:
                cluster.append(j)
                clustered.add(j)

        clusters.append(cluster)
    return clusters

def write_dup_summary(workdir, pdbs, clusters):
    with open(workdir + 'dup_summary.txt', 'w') as f:
        for clu in clusters:
            line = [pdbs[ind].getTitle() for ind in clu]
            f.write('\t'.join(line) + '\n')

def extract_rep_and_writepdb(pdbs, clusters, metal_sel, outdir):
    '''
    select represent pdb of each cluster and copy into a subfolder.
    TO DO: It is more efficient to copy the file into the subfolder.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sel_pdbs = []

    for clu in clusters:
        sel_pdbs.append((pdbs[clu[0]].getTitle(), pdbs[clu[0]]))

    superimpose_core_and_writepdb(sel_pdbs, sel_pdbs[0], metal_sel, outdir)


