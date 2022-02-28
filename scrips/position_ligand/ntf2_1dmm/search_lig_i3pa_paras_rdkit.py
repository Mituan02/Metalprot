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

        lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/i3pa_rdkit_fe.pdb'


    else:
        with open('/mnt/e/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
            ideal_alanine_bb_only = pickle.load(f)
        ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])

        path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'

        workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'

        lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/i3pa_rdkit_fe.pdb'

    return workdir, path_to_database, ideal_ala_coords, lig_path

task_type = 'search_unknow'
### Target strcture 
predefined_win_filters = [15, 19, 27]
lig_cgs = [None, None, None]

### Database para
use_enriched = True
use_abple=True

### Ligand 
find_unknow = True
ro1 = ['C2', 'C3']
rest1 = ['H1', 'C1', 'O1', 'O2', 'O3', 'FE1']
ro2 = ['C3', 'C4']
rest2 = ['H2', 'H3', 'H1', 'C2', 'C1', 'O1', 'O2', 'O3', 'FE1']

lig_connects = [['FE1', 'O2','O3']]  #On ligand
geo_sel = 'chid X and name FE1 O2 O3' #On protein
clash_dist = 3.0
lig_name = 'UNL'

### Search ligands paramters.
predefined_resnums = [11, 12, 16, 24, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119]

load_cg_aa_vdm_dict = {
    #'coo': [('AGKNQRSTY', True, False)], # (aas, filter_hb, filter_cc)
    #'bb_cco': [('AGKNQRSTY', True, False)],
    'indole': [('AGDEKNQRST', True, False), ('FWY', False, True)],
    'hid': [('AGDEKNQRST', True, False), ('FWY', False, True)],
    'hie': [('AGDEKNQRST', True, False), ('FWY', False, True)],
    'ph': [('FWY', False, True)]
}
rmsd = 0.6

vdm_cg_aa_cc_dict = {
    ('coo_0'):{
        'cg' : 'coo',
        'lgd_sel' : ['C2', 'C1', 'O2', 'O1'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD1', 'OD2']
    },
    ('coo_1'):{
        'cg' : 'coo',
        'lgd_sel' : ['C2', 'C1', 'O2', 'O1'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD2', 'OD1']
    },    
    ('coo_2'):{
        'cg' : 'coo',
        'lgd_sel' : ['C2', 'C1', 'O2', 'O1'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE1', 'OE2']
    },
    ('coo_3'):{
        'cg' : 'coo',
        'lgd_sel' : ['C2', 'C1', 'O2', 'O1'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE2', 'OE1']
    },
    ('indole_0'):{
        'cg' : 'indole',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'NE1',
        'correspond_resname' : 'TRP',
        'correspond_names' : ['NE1', 'CE2', 'CD1', 'CG']
    },
    ('indole_1'):{
        'cg' : 'indole',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'NE1',
        'correspond_resname' : 'TRP',
        'correspond_names' : ['NE1', 'CD1', 'CE2', 'CD2']
    },
    ('hid_0'):{
        'cg' : 'hid',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['ND1', 'CE1', 'CG', 'CD2']
    },
    ('hid_1'):{
        'cg' : 'hid',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['ND1', 'CG', 'CE1', 'NE2']
    },
    ('hie_0'):{
        'cg' : 'hie',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['NE2', 'CE1', 'CD2', 'CG']
    },
    ('hie_1'):{
        'cg' : 'hie',
        'lgd_sel' : ['N1', 'C6', 'C5', 'C4'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['NE2', 'CD2', 'CE1', 'ND1']
    },
    ('bb_cco_0'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_1'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_2'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C2', 'C1', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('ph_0'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CD1', 'CD2', 'CZ']
    },
    ('ph_1'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CG', 'CD2', 'CD1', 'CZ']
    },
    ('ph_2'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CG', 'CE1', 'CE2']
    },
    ('ph_3'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CE1', 'CG', 'CE2']
    },
    ('ph_4'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CD1', 'CZ', 'CD2']
    },
    ('ph_5'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE1', 'CZ', 'CD1', 'CD2']
    },
    ('ph_6'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CE1', 'CE2', 'CG']
    },
    ('ph_7'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CE2', 'CE1', 'CG']
    },
    ('ph_8'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CZ', 'CD2', 'CD1']
    },
    ('ph_9'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CE2', 'CD2', 'CZ', 'CD1']
    },
    ('ph_10'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CE2', 'CG', 'CE1']
    },
    ('ph_11'):{
        'cg' : 'ph',
        'lgd_sel' : ['C7', 'C8', 'C6', 'C10'],
        'represent_name' : 'CG',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CG', 'CE2', 'CE1']
    }
}