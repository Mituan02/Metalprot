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
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/1ukr/run_search_lig_indep.py local 1

'''

class Para:

    ### TaskType
    task_type = 'search_unknow' #>>> 'search_unknow', 'search_eval', 'search_exists'


    #>>> 1ukr
    predefined_win_filters = [('A',37), ('A',66), ('A',164)]

    predefined_resnums = [4, 6, 8, 9, 24, 34, 35, 37, 62, 64, 66, 68, 70, 77, 79, 81, 87, 109, 111, 113, 116, 119, 122, 125, 127, 129, 161, 162, 166, 168, 170]
    predefined_chidres = [('A', x) for x in predefined_resnums]

    #>>> Database para
    use_enriched = True
    use_abple=True
    rmsd = 0.75

    vdm_cg_aa_atommap_dict = {
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['S', 'O1', 'O2'],
            'correspond_resname' : 'ASP',
            'represent_name' : 'CG',
            'correspond_names' : ['CG', 'OD1', 'OD2'],
            'aas' : 'STYWHNQ', 
            'filter_hb' : True,
            'filter_cc' : False
        },   
        ('coo_1'):{
            'cg' : 'coo',
            'lgd_sel' : ['S', 'O1', 'O2'],
            'correspond_resname' : 'GLU',
            'represent_name' : 'CD',
            'correspond_names' : ['CD', 'OE1', 'OE2'],
            'aas' : 'STYWHNQ', 
            #'aas' : 'Y', 
            'filter_hb' : True, 
            'filter_cc' : False 
        },   
        ('bb_cco_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O1'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O1'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O1'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_3'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O2'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_4'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O2'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_5'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C7', 'S', 'O2'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHNQ',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O3', 'C12', 'N2'],
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
            'lgd_sel' : ['O3', 'C12', 'N2'],
            'correspond_resname' : 'GLN',
            'represent_name' : 'CD',
            'correspond_names' : ['OE1', 'CD', 'NE2'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },  
        ('bb_cco_1_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C11', 'C12', 'O3'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C11', 'C12', 'O3'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C11', 'C12', 'O3'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'STYWHQNDE',
            'filter_hb' : True,
            'filter_cc' : False
        }

    #     ('ph_0'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_1'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_2'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_3'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_4'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_5'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_6'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_7'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_8'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_9'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CZ',
    #         'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_10'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     },
    #     ('ph_11'):{
    #         'cg' : 'ph',
    #         'lgd_sel' : ['C7', 'C8', 'C15', 'C13'],
    #         'correspond_resname' : 'PHE',
    #         'represent_name' : 'CG',
    #         'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1'],
    #         'aas' : 'F',
    #         'filter_hb' : False,
    #         'filter_cc' : True
    #     }
    }


#>>> Run Search ligands.

def _filter_ligs_depre(outdir, ligs, target):
    lig = ligs[0] #>>> TO DO: somehow there is a bug for ligs '50g_md_0_*', '50g_md_6_*'
    tf = pr.calcTransformation(target.select('name OX ZN'), lig.select('name O4 ZN'))
    tf_rv = pr.calcTransformation(lig.select('name O4 ZN'), target.select('name OX ZN'))

    target_cp = target.copy()
    tf.apply(target_cp)
    filtered_ligs = position_ligand.ligand_clashing_filter(ligs, target_cp, dist = 2.7)

    # outdir_ori = outdir + 'origin_ligs/' 
    # os.makedirs(outdir_ori, exist_ok=True)
    # position_ligand.write_ligands(outdir_ori, filtered_ligs)
    # pr.writePDB(outdir + target.getTitle() + '_adj.pdb', target_cp)

    [tf_rv.apply(_lig) for _lig in filtered_ligs]
    position_ligand.write_ligands(outdir, filtered_ligs)
    
    return filtered_ligs


def filter_ligs(outdir, ligs, target):
    lig = ligs[0]
    coords = ligs[1]
    labels = ligs[2]

    R, m_com, t_com = transformation.get_rot_trans(lig.select('name O4 ZN').getCoords(), target.select('name OX ZN').getCoords())
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
    lig_path = '/mnt/e/DesignData/Metalloenzyme/belinostat/ligs/meo_50g_amber14eht_md_out/'

    #>>>
    workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/'
    
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
    path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'
    workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'
    #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'
    lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_rdkit.pdb'
    
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

