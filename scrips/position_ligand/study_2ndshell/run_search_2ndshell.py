'''
The script here is for the purposes of evaluate/search 2nd shell for metalloprotein.
'''

import os
import sys
from metalprot.combs import search_lig_indep_wrap 
import datetime

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/study_2ndshell/run_search_2ndshell.py local_single 1cbx.pdb
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/study_2ndshell/run_search_2ndshell.py
''' 

class Para:
    #>>> On wynton path
    #path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'
    #workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_2rd_sel/'

    #>>> On local path
    #path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'
    path_to_database='/mnt/e/DesignData/Database/vdMs_all/'

    #workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/'
    workdir = '/mnt/e/DesignData/LigandBB/study2ndshell/'

    #>>> There are three type of tasks: 'search_unknow', 'search_eval', or 'search_2ndshell'
    task_type = 'search_2ndshell'

    ### Target strctures
    #>>> Search ligands only with vdMs on the predefined chidres. If task_type is 'search_2ndshell' and predefined_chidres is None, the program will calc nearby residues.

    # predefined_win_filters = [('B', 8), ('B', 10), ('B', 178)]
    # lig_cg_2ndshell = [['hid', 'hie'],['hid', 'hie'], ['hid', 'hie']]

    #predefined_win_filters = [('A', 37), ('A', 66), ('A', 164)]
    #lig_cg_2ndshell = [['coo'], ['hid', 'hie'], ['hid', 'hie']]

    #>>> 5od1
    predefined_win_filters = [('A', 61)]
    lig_cg_2ndshell = [['hid', 'hie']]
    predefined_chidres = [('A', 58)]

    #>>> 6dwv
    predefined_win_filters = [('B', 8), ('B', 10), ('B', 178)]
    lig_cg_2ndshell = [['hid', 'hie'], ['hid', 'hie'], ['hid', 'hie']]    
    predefined_chidres = [('B', 69), ('B', 258), ('B', 179)]

    #>>> 1cbx
    predefined_win_filters = [('A', 69)]
    lig_cg_2ndshell = [['hid', 'hie']]
    predefined_chidres = [('A', 142)]

    # #>>> 1os0
    # predefined_win_filters = [('A', 142), ('A', 146)]
    # lig_cg_2ndshell = [['hid', 'hie'], ['hid', 'hie']]    
    # predefined_chidres = [('A', 165), ('A', 170)]

    # #>>> 2afw
    # predefined_win_filters = [('A', 159)]
    # lig_cg_2ndshell = [['coo']]
    # predefined_chidres = [('A', 140)]   

    # #>>> 3c56
    # predefined_win_filters = [('A', 83), ('A', 210)]
    # lig_cg_2ndshell = [['hid', 'hie'], ['hid', 'hie']]    
    # predefined_chidres = [('A', 104), ('A', 132)]

    # #>>> 3d4y
    # predefined_win_filters = [('A', 204), ('A', 471)]
    # lig_cg_2ndshell = [['coo'], ['hid', 'hie']]    
    # predefined_chidres = [('A', 167), ('A', 228)]

    # #>>> 3ebh
    # predefined_win_filters = [('A', 500), ('A', 519)]
    # lig_cg_2ndshell = [['hid', 'hie'], ['coo']]    
    # predefined_chidres = [('A', 503), ('A', 518)]

    # #>>> 4hwo
    # predefined_win_filters = [('A', 385), ('A', 511)]
    # lig_cg_2ndshell = [['hid', 'hie'], ['hid', 'hie']]    
    # predefined_chidres = [('A', 337), ('A', 486)]

    #>>> vdM database filter.
    use_enriched = True
    use_abple = True
    rmsd = 1.0

    #TO DO: rename the atoms
    vdm_cg_aa_atommap_dict = {
        ('coo_0'):{
            'cg' : 'coo',
            'lgd_sel' : ['CB', 'CG', 'OD1', 'OD2'],
            'represent_name' : 'OD2',
            'correspond_resname' : 'ASP',
            'correspond_names' : ['CB', 'CG', 'OD1', 'OD2'],
            'aas' : 'AGKRNQSTYH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('coo_1'):{
            'cg' : 'coo',
            'lgd_sel' : ['CG', 'CD', 'OE1', 'OE2'],
            'represent_name' : 'OE2',
            'correspond_resname' : 'GLU',
            'correspond_names' : ['CG', 'CD', 'OE1', 'OE2'],
            'aas' : 'AGKRNQSTYH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('hid_0'):{
            'cg' : 'hid',
            'lgd_sel' : ['ND1', 'CE1', 'CG', 'CD2'],
            'represent_name' : 'ND1',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['ND1', 'CE1', 'CG', 'CD2'],
            'aas' : 'AGDENQSTYH',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('hie_0'):{
            'cg' : 'hie',
            'lgd_sel' : ['NE2', 'CE1', 'CD2', 'CG'],
            'represent_name' : 'NE2',
            'correspond_resname' : 'HIS',
            'correspond_names' : ['NE2', 'CE1', 'CD2', 'CG'],
            'aas' : 'AGDENQSTYH',
            'filter_hb' : True,
            'filter_cc' : False
        }
    }


def run_local():
    para = Para()
    print('Task: ' + para.task_type)
    workdir, path_to_database, lig_path = para.workdir, para.path_to_database, None

    pdb_file = sys.argv[2]

    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = workdir + para.task_type  + '_result/output_' + '_'+ pdb_file + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)
    
    search_lig_indep_wrap.run_search_2ndshell(pdb_file, workdir, outdir, path_to_database, lig_path, para)
    return

def run_local_all():
    para = Para()
    print('Task: ' + para.task_type)
    workdir, path_to_database, lig_path = para.workdir, para.path_to_database, None

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    for pdb_file in pdb_files:
        outdir = workdir + para.task_type  + '_result/output_' + '_'+ pdb_file + '_' + time_tag + '/'
        os.makedirs(outdir, exist_ok=True)    
        search_lig_indep_wrap.run_search_2ndshell(pdb_file, workdir, outdir, path_to_database, lig_path, para)
    return

if __name__=='__main__':
    if sys.argv[1] == 'local_single':
        run_local()
    elif sys.argv[1] == 'local_all':
        run_local_all()