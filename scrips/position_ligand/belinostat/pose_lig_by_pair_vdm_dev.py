'''
The script here is for the development of the pose_lig_by_pair_vdm.py.
Use ABLE (6w70.pdb) as an example.
Use the run_verify_vdMs() function to find the solution vdMs to prove the existance of the vdMs.
'''

import sys
import prody as pr

sys.path.append(r'/mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/')
import pose_lig_by_pair_vdm as plby

from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep

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

    select_chidres_keys = plby._select_chidres_keys(workdir, target, lig, para, path_to_database)
    for key_a, key_b, chidres_a, chidres_b, abple_a, abple_b in select_chidres_keys:
        plby._search_select_pair_vdm(workdir, target, lig, para, path_to_database, key_a, key_b, chidres_a, chidres_b, abple_a, abple_b)
    return


def _verify_vdMs(workdir, target, chidres, para, abples, path_to_database, lig, vdm_cg_aa_atommap_dict):
    for cg_id in vdm_cg_aa_atommap_dict.keys():

        pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
        abple = abples[pos.getResindices()[0]]

        for a in vdm_cg_aa_atommap_dict[cg_id]['aas']:
            aa = constant.inv_one_letter_code[a]

            df_vdms = search_lig_indep.load_old_vdm(path_to_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
            df_vdm_filter = search_lig_indep.filter_db(df_vdms, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)   

            ligs = [lig]
            ideal_ala_coords = constant.ideal_ala_coords
            results = search_lig_indep.search_lig_at_cg_aa_resnum(target, chidres, pos, abple, ligs, vdm_cg_aa_atommap_dict, cg_id, df_vdm_filter, ideal_ala_coords, para.rmsd, True, False)
            print(len(results))

            if len(results) <= 0:
                continue

            cg = vdm_cg_aa_atommap_dict[cg_id]['cg']
            for cg_id, lig_id, _rmsd, vdm_ag, vdm_id, score, v, contact_hb, contact_cc in results:
                prefix = 'Lig-' + str(lig_id) + '_' + cg + '_' + chidres[0] + str(chidres[1]) + '_' + aa + '_rmsd_' + str(round(_rmsd, 2)) + '_v_' + str(round(score, 1)) + '_'       

                vdm_ag.setTitle(prefix + '_' + vdm_ag.getTitle())            
                pr.writePDB(workdir + vdm_ag.getTitle(), vdm_ag)  
    return

def run_verify_vdMs():
    '''
    For a protein bind to ligand, check the vdMs that satisfy the binding. 
    '''
    workdir = '/mnt/e/DesignData/Metalloenzyme/6w70_vdM/'

    path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

    lig = pr.parsePDB(workdir + 'gg2.pdb')

    target = pr.parsePDB(workdir + '6w70_bb.pdb')

    para = Para()

    abples, phipsi = utils.seq_get_ABPLE(target)

    _verify_vdMs(workdir, target, para.predefined_resnums[0], para, abples, path_to_database, lig, para.vdm_cg_aa_atommap_dict_a)

    _verify_vdMs(workdir, target, para.predefined_resnums[1], para, abples, path_to_database, lig, para.vdm_cg_aa_atommap_dict_b)

    return

if __name__=='__main__':
    run_verify_vdMs()
    run_local()