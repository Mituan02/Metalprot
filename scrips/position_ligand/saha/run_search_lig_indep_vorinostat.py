'''
The script is to design ligand-metal binding enzyme using the vdM database. 
The ligand is superimposed on the metal with all possible conformers. 
And the ligs (different conformers) are then searched against the vdM database.
'''

import os
import sys

import prody as pr
import numpy as np
import pickle
import datetime

from metalprot.basic import transformation
from metalprot.combs import search_lig_indep_wrap, position_ligand
'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/saha/run_search_lig_indep_vorinostat.py local 1

'''

class Para:

    ### TaskType
    task_type = 'search_unknow' #>>> 'search_unknow', 'search_eval', 'search_exists'

    #>>> Helix6a
    predefined_win_filters = [('A',17), ('A',21), ('A',132)]
    #predefined_resnums = [7, 14, 17, 21, 58, 61, 62, 65, 68, 69, 72, 77, 78, 80, 81, 84, 85, 88, 89, 91, 92, 95, 132, 135, 138, 139, 142]
    predefined_resnums = [24, 25, 51, 54, 58, 61, 88, 91, 92, 95, 99, 124, 125, 135]
    predefined_chidres = [('A', x) for x in predefined_resnums]

    #>>> Database para
    use_enriched = True
    use_abple=True
    rmsd = 0.75

    vdm_cg_aa_atommap_dict = {
        ('bb_cco_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C6', 'C8', 'O2'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C6', 'C8', 'O2'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C6', 'C8', 'O2'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O2', 'C8', 'N2'],
            'correspond_resname' : 'ASN',
            'represent_name' : 'CG',
            'correspond_names' : ['OD1', 'CG', 'ND2'],
            'aas' : 'STYWHQNDE',
            #'aas' : 'E',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('conh2_1'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O2', 'C8', 'N2'],
            'correspond_resname' : 'GLN',
            'represent_name' : 'CD',
            'correspond_names' : ['OE1', 'CD', 'NE2'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('bb_cnh_0'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C8', 'N2', 'H14'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('bb_cnh_1'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C8', 'N2', 'H14'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('bb_cnh_2'):{
            'cg' : 'bb_cnh',
            'lgd_sel' : ['C8', 'N2', 'H14'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'N', 'H'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('coh_0'):{
            'cg' : 'coh',
            'lgd_sel' : ['N2', 'O3'],
            'correspond_resname' : 'SER',
            'represent_name' : 'CB',
            'correspond_names' : ['CB', 'OG'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('coh_1'):{
            'cg' : 'coh',
            'lgd_sel' : ['N2', 'O3'],
            'correspond_resname' : 'THR',
            'represent_name' : 'CB',
            'correspond_names' : ['CB', 'OG1'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
    }


#>>> Run Search ligands.
def filter_ligs(outdir, ligs, target):
    lig = ligs[0]
    coords = ligs[1]
    labels = ligs[2]

    R, m_com, t_com = transformation.get_rot_trans(lig.select('name O3 ZN').getCoords(), target.select('name OX ZN').getCoords())
    coords = np.dot(np.array(coords) - m_com, R) + t_com

    target_coords = target.select('name N C CA O CB').getCoords()
    lig_size = len(lig.select('heavy'))
    successed_id = position_ligand._ligand_clashing_filteredid(target_coords, coords, lig_size, dist = 3.0)

    coord_reshpe = coords.reshape((-1, lig_size,3))
    filtered_ligs = []
    for i in successed_id:
        lig_cp = lig.select('heavy').copy()
        lig_cp.setCoords(coord_reshpe[i,:,:])
        lig_cp.setTitle(labels[i])
        filtered_ligs.append(lig_cp)

    position_ligand.write_ligands(outdir, filtered_ligs)
    
    return filtered_ligs

def run_all(file, workdir, path_to_database, lig_path, para):
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 

    target_path = workdir + file
    outdir = workdir + para.task_type + '_'+ file.split('.')[0] + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)

    target = pr.parsePDB(target_path)
    ligs = []
    for _f in os.listdir(lig_path):
        if 'coord.pkl' not in _f:
            continue
        print('Load ligs: ' + _f)
        with open(lig_path + _f, 'rb') as f:
            all_ligs = pickle.load(f)
         #>>> Filter ligs.
        _ligs = filter_ligs(outdir, all_ligs, target)
        ligs.extend(_ligs)

    #target, chidres2ind = search_lig_indep.prepare_rosetta_target(outdir, target_path, para.predefined_win_filters)
    
    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)
    print('number of ligs: {}'.format(len(ligs)))

    lig_vdm_dict = search_lig_indep_wrap.run_search(target, ligs, path_to_database, para, para.predefined_chidres, lig_cg_2ndshell = None)
    
    print('lig_vdm_dict size {}'.format(len(list(lig_vdm_dict.keys()))))
    search_lig_indep_wrap.write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict,  para)

    return


def run_local():
    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'
    lig_path = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/SuperLig/'

    #>>> Workdir
    workdir = '/mnt/e/DesignData/Metalloenzyme/HelixZn_sel/helix1_17-21-132/'
    
    para = Para()

    print('Task: ' + para.task_type)

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[2]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para)
    return


def run_wynton():
    path_to_database='/wynton/home/degradolab/lonelu/DesignData/Database/vdMs/'
    workdir = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/SAHA_Vorinostat/helix1_17-21-132/'
    lig_path = '/wynton/home/degradolab/lonelu/DesignData/Metalloenzyme/SAHA_Vorinostat/SuperLig/'
    
    para = Para()
    print('Task: ' + para.task_type)

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[2]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para)
    return

if __name__=='__main__':
    if sys.argv[1] == 'wynton': 
        run_wynton()
    elif sys.argv[1] == 'local':
        run_local()

