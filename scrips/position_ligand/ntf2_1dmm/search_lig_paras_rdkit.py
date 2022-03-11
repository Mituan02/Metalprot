import os
import pickle
import numpy as np

def get_file_path(on_wynton):
    if on_wynton:
        with open('/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
            ideal_alanine_bb_only = pickle.load(f)
        ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])

        path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'

        workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'
        #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'

        lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_rdkit.pdb'

    else:
        with open('/mnt/e/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
            ideal_alanine_bb_only = pickle.load(f)
        ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])

        path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'
        #workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'
        #workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_2ndRound/output_F55D_sel/'
        workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r1/output_sel/'
        lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_rdkit.pdb'

    return workdir, path_to_database, ideal_ala_coords, lig_path

### TaskType
task_type = 'search_unknow'
### Target strcture 
#predefined_win_filters = [15, 19, 27]
predefined_win_filters = [85, 87, 100]
lig_cgs = [None, None, None]

### Database para
use_enriched = True
use_abple=True

### Ligand 
ro1 = ['C8', 'C7']
rest1 = ['H7', 'C9', 'O1', 'O3', 'O2', 'FE1']
ro2 = ['C7', 'C6']
rest2 = ['H5', 'H6', 'H7', 'C8', 'C9', 'O1', 'O3', 'O2', 'FE1']

lig_connects = [['FE1', 'O3','O1']]
geo_sel = 'chid X and name FE1 O2 O3'
clash_dist = 3.0
lig_name = 'TTS'

### Search ligands paramters.
#predefined_resnums = [11, 12, 16, 24, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119] 
predefined_resnums = [11, 12, 15, 16, 19, 24, 27, 30, 31, 35, 37, 39, 43, 46, 52, 55, 56, 59, 60, 65, 67, 69, 83, 89, 98, 102, 104, 106, 112, 115, 117, 119]


load_cg_aa_vdm_dict = {
    'coo': [('AGKNQRSTY', True, False)], # (aas, filter_hb, filter_cc)
    'bb_cco': [('AGKNQRSTY', True, False)],
    'phenol': [('AGDEKNQRST', True, False), ('FWY', False, True)],
    'ph': [('FWY', False, True)]
}

rmsd = 0.6

vdm_cg_aa_cc_dict = {
    ('coo_0'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD1', 'OD2']
    },
    ('coo_1'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD2', 'OD1']
    },    
    ('coo_2'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE1', 'OE2']
    },
    ('coo_3'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O2'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE2', 'OE1']
    },
    ('phenol_0'):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O4'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE1', 'CZ', 'CE2', 'OH']
    },
    ('phenol_1'):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O4'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE2', 'CZ', 'CE1', 'OH']
    },
    ('bb_cco_0'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_1'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_2'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('ph_0'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ']
    },
    ('ph_1'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ']
    },
    ('ph_2'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2']
    },
    ('ph_3'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2']
    },
    ('ph_4'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2']
    },
    ('ph_5'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2']
    },
    ('ph_6'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG']
    },
    ('ph_7'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG']
    },
    ('ph_8'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1']
    },
    ('ph_9'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1']
    },
    ('ph_10'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1']
    },
    ('ph_11'):{
        'cg' : 'ph',
        'lgd_sel' : ['C1', 'C2', 'C6', 'C4'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1']
    }
}