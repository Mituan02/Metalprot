'''
20240501
Kaipeng's'four helix bundle, KH207 etc.
Find the best position for porphyrin.

python /Users/lonelu/GitHub_Design/Metalprot/scrips/porphyrin/position_porphyrin_KH.py
'''

import time
import prody as pr
import os
import numpy as np
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
from sklearn.neighbors import NearestNeighbors


def getOrderedCoords(pdb, sel_str):
    _x = []
    for _y in sel_str:
        _x.append(pdb.select('name ' + _y).getCoords()[0])
    sel_coords = np.array(_x)
    return sel_coords

def rotate_by_ro(lig, rot, rot_range):
    '''
    Note that the FE must be at (0, 0, 0).
    To test:
    #all_ligs = rotate_by_ro(lig, ro1, ro1_range)
    '''
    # transform the lig to z axis.
    rot_coords = getOrderedCoords(lig, rot)
    _rot_coords = (rot_coords[1] - rot_coords[0])/norm(rot_coords[1] - rot_coords[0])
    all_ligs = []
    for i in rot_range:
        _lig = lig.copy()
        rotation = Rotation.from_rotvec(np.radians(i)*_rot_coords)       
        _coords = rotation.apply(_lig.getCoords())
        _lig.setCoords(_coords)
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)
        all_ligs.append(_lig)
    return all_ligs



def generate_rotated_porphyrins(lig, ro1, ro2, ro1_range, ro2_range):
    all_ligs = []
    rot1_ligs = rotate_by_ro(lig, ro1, ro1_range)
    for lig1 in rot1_ligs:
        _ligs = rotate_by_ro(lig1, ro2, ro2_range)
        all_ligs.extend(_ligs)
    return all_ligs

def porphyrin2motif(lig, ligs, sel_string, motif_coord):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''

    lig_sel = getOrderedCoords(lig, sel_string)

    transformation = pr.calcTransformation(lig_sel, motif_coord)

    for lg in ligs:
        transformation.apply(lg)
    return 

def _ligand_clashing_filteredid(target_coords, lig_coords, lig_len, dist):

    nbrs = NearestNeighbors(radius= dist).fit(target_coords)
    adj_matrix = nbrs.radius_neighbors_graph(lig_coords).astype(bool)

    adj_matrix_reshape = adj_matrix.reshape((-1, adj_matrix.shape[1]*lig_len))
    successed = []
    #>>> TO DO: resahpe the adj_matrix (type csr_matrix) will improve the code.
    for i in range(adj_matrix_reshape.shape[0]):
        if adj_matrix_reshape.getrow(i).toarray().any():
            continue
        successed.append(i)

    return successed

def ligand_clashing_filter(ligs, target, lig_sel = 'heavy', tar_sel = 'name N C CA O CB', dist = 3):
    '''
    The ligand clashing: the ligs cannot have any heavy atom within 3 A of a target bb.
    Nearest neighbor is used to calc the distances. 
    '''
    all_coords = []
    ids = []

    for i in range(len(ligs)):
        coords = ligs[i].select(lig_sel).getCoords()
        all_coords.extend(coords)
        ids.extend([i for j in range(len(coords))])

    lig_len = len(ligs[0].select(lig_sel))
    target_coords = target.select(tar_sel).getCoords()
    
    successed_id = _ligand_clashing_filteredid(target_coords, all_coords, lig_len, dist)
    
    filtered_ligs = []
    for i in successed_id:
        filtered_ligs.append(ligs[i])

    return filtered_ligs

################################################################
#Case specific function

def sample_porphyrin_pos(target, tar_sel, lig, ro1, ro2, ro1_range, ro2_range, dist = 2.7, lig_sel = 'not protein and heavy'):

    motif= target.select(tar_sel)
    motif_coord = getOrderedCoords(motif, ['CG', 'ND1', 'CE1', 'CD2', 'NE2'])   
    
    pr.calcTransformation(lig.select('name FE NI CU ZN').getCoords(), np.zeros((1, 3), dtype=float)).apply(lig)  ##Transform the metal to (0, 0, 0)

    all_ligs = generate_rotated_porphyrins(lig, ro1, ro2, ro1_range, ro2_range)
    print('generate {} lig pos'.format(len(all_ligs)))

    sel_string = ['CG', 'ND1', 'CE1', 'CD2', 'NE2']
    porphyrin2motif(lig, all_ligs, sel_string, motif_coord)

    # os.makedirs(workdir + 'all_ligs/', exist_ok=True)
    # for lg in all_ligs:
    #     pr.writePDB(workdir + 'all_ligs/' +  lg.getTitle(), lg)

    filtered_ligs = ligand_clashing_filter(all_ligs, target, lig_sel, dist = dist)
    print('filter to {} lig pos'.format(len(filtered_ligs)))

    return filtered_ligs



################################################################
#vdM sampling function

import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms as search_cg_vdms

pd.set_option("display.max_columns", None)


vdm_cg_aa_atommap_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O1A', 'CGA', 'O2A'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1', 'CG', 'OD2'],
    },    
        (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['O1A', 'CGA', 'O2A'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD2', 'CG', 'OD1'],
    },    
        (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['O1D', 'CGD', 'O2D'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1', 'CD', 'OE2'],
    },    
        (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O1D', 'CGD', 'O2D'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2', 'CG', 'OE1'],
    },    
}


def run_vdm_sample(target, resnum, aa, ligands, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir):

    labels_cgs = {}
    df_cgs = {}
    dist_ind_cgs = {}

    cg_ids = vdm_cg_aa_atommap_dict.keys()
    for cg_id in cg_ids:
        #Load vdm database.
        #Here we specify to check ASP vdMs. Note this can be programed to be easier.
        df_vdm = gvdm_helper.load_old_vdm(path_to_vdm_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
        
        #Transformation.
        pos = target.select('name N CA C and resnum ' + str(resnum))
        tf = pr.calcTransformation(constant.ideal_ala_coords, pos.getCoords())
        df_vdm.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm[['c_x', 'c_y', 'c_z']].to_numpy())

        #Test if the first vdm's coords are changed.
        #x = labels_cgs[cg_id].iloc[0]
        #v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
        #ag = gvdm_helper.df2ag(v, 'test_position_of_1st_vdm')
        #pr.writePDB(outdir + ag.getTitle(), ag)

        #Search vdms.
        search_cg_vdms.search_vdm(df_vdm, ligands, cg_id, vdm_cg_aa_atommap_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.5)
        

    CgCombInfoDict = search_cg_vdms.construct_vdm_write(outdir, ligands, labels_cgs, vdm_cg_aa_atommap_dict, df_cgs, dist_ind_cgs, clash_radius = 2.5)

    #search_cg_vdms.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')
    return


def sample_vdMs(target, ligs, path_to_vdm_database, outdir):
    resnums = [3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 17, 18, 57, 61, 62, 63, 64, 66, 67, 111]
    aas = ['ALA', 'GLY', 'SER', 'THR', 'TYR', 'ASN', 'GLN','ARG', 'LYS', 'HIS']

    for i in range(len(resnums)):
        resnum = resnums[i]
        for aa in aas:
            _outdir = outdir + 'v_' + str(resnum) + '_' + aa + '_'
            
            run_vdm_sample(target, resnum, aa, ligs, path_to_vdm_database, vdm_cg_aa_atommap_dict, _outdir)
    return

################################################################

import multiprocessing as mp

def xxx(target_path, hem_database, tar_sel, ro1, ro2, ro1_range, ro2_range, dist, outdir, path_to_vdm_database):
    target = pr.parsePDB(target_path)
    print('here')
    all_filtered_ligs = []
    for lig_path in os.listdir(hem_database):
        if '.pdb' not in lig_path:
            continue
        
        lig = pr.parsePDB(hem_database + lig_path)
        #Get filtered ligs.
        try:
            filter_ligs = sample_porphyrin_pos(target, tar_sel, lig, ro1, ro2, ro1_range, ro2_range, dist = dist)
        except:
            print('---------------------------------------------------')
            print(lig.getTitle())
        all_filtered_ligs.extend(filter_ligs)
    #Sample vdMs.
    sample_vdMs(target, all_filtered_ligs, path_to_vdm_database, outdir)
        
    return

def main():

    hem_database = '/Users/lonelu/DesignData/ligands_metal/porphyrin/pdbs/M1_AAMetal_HIS_reps/'

    workdir = '/Users/lonelu/DesignData/labmate_runs/KH_heme_240501/'

    outdir = workdir + 'output_20240502_KH212/'
    os.makedirs(outdir, exist_ok=True)

    ro1 = ['FE', 'NA']
    ro1_range = range(-18, 19, 3) 
    ro2 = ['FE', 'NE2']
    ro2_range = range(0, 360, 5) 

    dist = 2.3
    tar_sel = 'resnum 50 and name CG ND1 CE1 CD2 NE2'

    path_to_vdm_database='/Users/lonelu/DesignData/Combs/vdMs/'


    # for target_path in os.listdir(workdir):
    #     if '.pdb' not in target_path:
    #         continue
    #     target = pr.parsePDB(workdir + target_path)

    target_path = workdir + 'KH212.pdb'
    
    xxx(target_path, hem_database, tar_sel, ro1, ro2, ro1_range, ro2_range, dist, outdir, path_to_vdm_database)
    # num_cores = int(mp.cpu_count()- 2)
    # pool = mp.Pool(num_cores)

    # [pool.apply_async(xxx, args=(target_path, hem_database, tar_sel, ro1, ro2, ro1_range, ro2_range, dist, outdir, path_to_vdm_database)) for target_path in [target_path]]
    # pool.close()
    # pool.join()

    return

if __name__=='__main__':
    main()

###########################################################
'''
#Test case

workdir = '/Users/lonelu/DesignData/APEX/Design_210/'

outdir = workdir + 'output_240823/'

ro1 = ['FE', 'NA']
ro1_range = range(-24, 25, 8) 
ro2 = ['FE', 'NE2']
ro2_range = range(0, 360, 20) 

target_path = workdir + 'KH210.pdb'
target = pr.parsePDB(target_path)
tar_sel = 'resnum 50 and name CG ND1 CE1 CD2 NE2'

lig_path = '/Users/lonelu/DesignData/APEX/Design_210/1a2f_hem_oo.pdb'
lig = pr.parsePDB(lig_path)

# all_ligs = generate_rotated_porphyrins(lig, ro1, ro2, ro1_range, ro2_range)

filter_ligs = sample_porphyrin_pos(target, tar_sel, lig, ro1, ro2, ro1_range, ro2_range, dist =1.8, lig_sel = 'chain B or chain C')

os.makedirs(workdir + 'filtered_ligs/', exist_ok=True)
for lg in filter_ligs:
    pr.writePDB(workdir + 'filtered_ligs/' +  lg.getTitle(), lg)

path_to_vdm_database='/Users/lonelu/DesignData/Combs/vdMs/'
sample_vdMs(filter_ligs, path_to_vdm_database, outdir)



'''
