'''
The script here is for the development of the pose_lig_by_pair_vdm.py.
Use ABLE (6w70.pdb) as an example.
Use the run_verify_vdMs() function to find the solution vdMs to prove the existance of the vdMs.
'''
import os
import sys
import prody as pr
import datetime

#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/')

from metalprot.basic import utils
from metalprot.combs import search_lig_indep_wrap, search_lig_indep_inpair

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/pose_lig_by_pair_vdm_dev.py
'''

class Para:

    resnums = [49, 6]
    predefined_resnums = [('A', r) for r in resnums]

    use_enriched = True
    use_abple=True

    rmsd = 0.75

    vdm_cg_aa_atommap_dict_a = {
        ('bb_cco_0'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C13', 'C8', 'O3'],
            'correspond_resname' : 'GLY',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'H',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_1'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C13', 'C8', 'O3'],
            'correspond_resname' : 'ALA',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'H',
            'filter_hb' : True,
            'filter_cc' : False
        },
        ('bb_cco_2'):{
            'cg' : 'bb_cco',
            'lgd_sel' : ['C13', 'C8', 'O3'],
            'correspond_resname' : 'PRO',
            'represent_name' : 'CA',
            'correspond_names' : ['CA', 'C', 'O'],
            'aas' : 'H',
            'filter_hb' : True,
            'filter_cc' : False
        },
    }

    vdm_cg_aa_atommap_dict_b = {
        ('conh2_0'):{
            'cg' : 'conh2',
            'lgd_sel' : ['O1', 'C11', 'N3'],
            'correspond_resname' : 'ASN',
            'represent_name' : 'OD1',
            'correspond_names' : ['OD1', 'CG', 'ND2'],
            'aas' : 'YQ', #'HDE' is also provide good hb but we want to get rid of confusing the metal binding.
            'filter_hb' : True,
            'filter_cc' : False
        },  
    }


def run_local():
    workdir = '/mnt/e/DesignData/Metalloenzyme/6w70_vdM/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'gg2.pdb')

    target = pr.parsePDB(workdir + '6w70_bb.pdb')

    para = Para()

    outdir = workdir + 'output_pair-search_' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    os.makedirs(outdir, exist_ok = True)

    select_chidres_keys = search_lig_indep_inpair._select_chidres_keys(target, lig, para, path_to_database)
    for key_a, key_b, chidres_a, chidres_b, abple_a, abple_b in select_chidres_keys:
        search_lig_indep_inpair.search_select_pair_vdm(outdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return

def run_verify_vdMs():
    '''
    For a protein bind to ligand, check the vdMs that satisfy the binding. 
    '''
    workdir = '/mnt/e/DesignData/Metalloenzyme/6w70_vdM/'

    outdir = workdir + 'output_verify-vdM_' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    os.makedirs(outdir, exist_ok = True)

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'gg2.pdb')

    target = pr.parsePDB(workdir + '6w70_bb.pdb')

    para = Para()

    abples, phipsi = utils.seq_get_ABPLE(target)

    search_lig_indep_wrap.verify_vdMs(outdir, target, para.predefined_resnums[0], para, abples, path_to_database, lig, para.vdm_cg_aa_atommap_dict_a)

    search_lig_indep_wrap.verify_vdMs(outdir, target, para.predefined_resnums[1], para, abples, path_to_database, lig, para.vdm_cg_aa_atommap_dict_b)

    return

if __name__=='__main__':
    run_verify_vdMs()
    run_local()