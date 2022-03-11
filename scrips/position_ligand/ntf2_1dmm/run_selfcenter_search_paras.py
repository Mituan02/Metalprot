import os
import pickle
import numpy as np

def get_file_path(on_wynton):
    if on_wynton:
        #query_dir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/database/pickle_all_fe_220119/'
        query_dir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/database/pickle_noCYS_mn_fe_co_220119/'
        workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/data/ntf2_fe/family_3vsy/'

    else:
        #query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/pickle_noCYS/'
        query_dir = '/mnt/e/DesignData/ligands/all/20220116_FE_MN_CO/20220116_selfcenter/pickle_noCYS/'
        workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_3rdRound/output_55F_newlop/output_55F_newlop_sel/eval/'

    return query_dir, workdir


win_filters = [[15, 19, 27] for i in range(4)]
#win_filters = [[16, 20, 28] for i in range(3)]
#win_filters = [[] for i in range(len(pdb_files))]
#win_filters = predefined_win_filters

geometry_path = None
geometry_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/fe_geo.pdb'


metal_metal_dist = 0.75

num_contact_vdms = [3]

allowed_aa_combinations = [['H', 'H', 'D'], ['H', 'H', 'E']] 
#allowed_aa_combinations = []