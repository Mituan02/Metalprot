'''
The script is to design ligand-metal binding enzyme using the vdM idea. 
Specially, here is to search the combs vdM library without using Combs2.
'''

import os
import sys
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant
from metalprot.basic import utils
from metalprot.combs import search_lig_indep
from metalprot.combs import position_ligand
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import gc
import shutil
import datetime

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_ligand_indep.py
'''

with open('/mnt/e/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
#with open('/wynton/home/degradolab/lonelu/GitHub_Design/Metalprot/metalprot/constants/ideal_alanine_bb_only.pkl', 'rb') as f:
    ideal_alanine_bb_only = pickle.load(f)
ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])


path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'
#path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/'


def run_ligand(outdir, target):
    '''
    Generate all potential ligands for each binding position. 
    '''
    lig = pr.parsePDB(lig_path)

    all_ligs = position_ligand.generate_rotated_ligs(lig, [ro1, ro2], [rest1, rest2], [5, 5])

    # points = np.array(position_ligand.fibonacci_sphere(10, scale=0.2))
    # point_sel = 'name FE1'

    filtered_ligs, _ = position_ligand.generate_ligands(all_ligs, target, lig_connects, geo_sel, clash_dist = clash_dist)

    position_ligand.write_ligands(outdir, filtered_ligs)

    return


def run_search(target, ligs):
    '''
    cg = 'phenol'
    resnum = 102
    '''
    # load_cg_aa_vdm_dict = {'phenol': ('FW', True, True)}
    # predefined_resnums = [102]
    abples, phipsi = utils.seq_get_ABPLE(target)

    #summaries = []
    lig_vdm_dict = {}
    for cg in load_cg_aa_vdm_dict.keys():
        for cg_aa_vdms in load_cg_aa_vdm_dict[cg]:
            for a in cg_aa_vdms[0]:
                aa = constant.inv_one_letter_code[a]
                #df_vdm, df_gr, df_score = load_new_vdm(path_to_database, cg, aa)
                df_vdm = search_lig_indep.load_old_vdm(path_to_database, cg, aa)

                for resnum in predefined_resnums:
                    print('Searching: ' + cg + ' ' + aa + ' ' + str(resnum))
                    pos = target.select('resnum ' + str(resnum))
                    abple = abples[pos.getResindices()[0]]
                    df_vdm_filter = search_lig_indep.filter_db(df_vdm, abple=abple)        

                    results = search_lig_indep.search_lig_at_cg_aa_resnum(target, resnum, pos, abple, ligs, input_dict, cg, df_vdm_filter, ideal_ala_coords, rmsd, cg_aa_vdms[1], cg_aa_vdms[2])


                    for cg_id, lig_id, _rmsd, vdm_ag, vdm_id, score, v, contact_hb, contact_cc in results:
                        prefix = 'Lig-' + str(lig_id) + '_' + cg  + '_' + aa + '_' + str(resnum) + '_rmsd_' + str(round(_rmsd, 2)) + '_v_' + str(round(score, 1)) + '_'       

                        #pr.writePDB(outdir_all + prefix + str(vdm_id) + '.pdb.gz', vdm_ag)
                        #summaries.append((cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score))

                        vdm_ag.setTitle(prefix + str(vdm_id))                       
                        if lig_id in lig_vdm_dict.keys():
                            lig_vdm_dict[lig_id].append((vdm_ag, (cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc)))
                        else:
                            lig_vdm_dict[lig_id] = [(vdm_ag, (cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc))]

                        
                del [[df_vdm]]
                gc.collect()
                df_vdm = pd.DataFrame()
        
    return lig_vdm_dict


########################################################################
### Position ligand paramters.

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/'
#workdir = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe_1dmm_rosetta_sel/'


predefined_win_filters = [15, 19, 27]

lig_path = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/tts_fe_adj.pdb'
#lig_path = '/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/ntf2_fe/tts_fe_adj.pdb'

ro1 = ['C8', 'C7']
rest1 = ['C9', 'O1', 'O3', 'O4', 'FE1']
ro2 = ['C7', 'C6']
rest2 = ['H6', 'H7', 'C8', 'C9', 'O1', 'O3', 'O4', 'FE1']

lig_connects = [['FE1', 'O3','O1'], ['FE1', 'O1','O3']]
geo_sel = 'chid X and name FE1 O2 O3'
clash_dist = 3.0

########################################################################
### Search ligands paramters.

load_cg_aa_vdm_dict = {
    'coo': [('AGKNQRSTY', True, False)], # (aas, filter_hb, filter_cc)
    'bb_cco': [('AGKNQRSTY', True, False)],
    'phenol': [('AGDEKNQRST', True, False), ('FWY', False, True)],
    'ph': [('FWY', False, True)]
}

predefined_resnums = [11, 12, 15, 16, 19, 24, 27, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119]

input_dict = {
    ('coo_0'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD1', 'OD2']
    },
    ('coo_1'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OD2',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['CB', 'CG', 'OD2', 'OD1']
    },    
    ('coo_2'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE1', 'OE2']
    },
    ('coo_3'):{
        'cg' : 'coo',
        'lgd_sel' : ['C8', 'C9', 'O3', 'O4'],
        'represent_name' : 'OE2',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['CG', 'CD', 'OE2', 'OE1']
    },
    ('phenol_0'):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE1', 'CZ', 'CE2', 'OH']
    },
    ('phenol_1'):{
        'cg' : 'phenol',
        'lgd_sel' : ['C2', 'C3', 'C4', 'O2'],
        'represent_name' : 'OH',
        'correspond_resname' : 'TYR',
        'correspond_names' : ['CE2', 'CZ', 'CE1', 'OH']
    },
    ('bb_cco_0'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O4'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_1'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O4'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O']
    },
    ('bb_cco_2'):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C8', 'C9', 'O4'],
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

rmsd = 0.6

########################################################################
### Run Search ligands.

def run_all(file):
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 

    target_path = workdir + file
    #target_path = workdir + 'o2_1dmm_16-20-28_H-H-D_a_842.pdb'
    outdir = workdir + 'output_' + file + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)

    target = search_lig_indep.prepare_rosetta_target(outdir, target_path, predefined_win_filters)
    #target = pr.parsePDB(workdir + 'o1_1dmm_16-20-28_H-H-D_a_820_gly.pdb.pdb')
    run_ligand(outdir, target)

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)

    outdir_ligs = outdir + 'ligs_inorder/'
    os.makedirs(outdir_ligs, exist_ok=True)
    ligs = []
    lig_id = 0
    for file in os.listdir(outdir + 'filtered_ligs/'):
        if not '.pdb' in file and not '.pdb.gz' in file:
            continue
        lig = pr.parsePDB(outdir + 'filtered_ligs/' + file)
        shutil.copy2(outdir + 'filtered_ligs/' + file, outdir_ligs + 'Lig-'+ str(lig_id) + '_'+ file)
        lig_id += 1
        ligs.append(lig)

    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    lig_vdm_dict = run_search(target, ligs)

    ### write filtered vdms.
    summaries = []
    for lig_id in lig_vdm_dict.keys():
        values = lig_vdm_dict[lig_id]
        exist_required_cg_contact = False
        for v in values:
            if v[1][0] == 'phenol' and v[1][-2] == True:
                exist_required_cg_contact = True

        if not exist_required_cg_contact:
            continue

        for vdm_ag, _info in values:
            pr.writePDB(outdir_all + vdm_ag.getTitle() + '.pdb.gz', vdm_ag)
            summaries.append(_info)

    if len(summaries) <= 0:
        print('Filtered vdms is 0.')
        return

    with open(outdir + target.getTitle() + '_summary.tsv', 'w') as f:
        f.write('file\tcg\taa\tpos\tcg_tag\tlig\trmsd\tvdm_id\tscore\tvdm\tContact_hb\tContact_cc\n')
        for s in summaries:
            #print(s)
            f.write(target.getTitle() + '\t' + '\t'.join([str(x) for x in s]) + '\n')


    ### Write unique vdms.
    unique_vdms = {}
    for lig_id in lig_vdm_dict.keys():
        values = lig_vdm_dict[lig_id]
        for v in values:
            vdm_ag = v[0]
            (cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc) = v[1]
            if (cg, aa, vdm_id) in unique_vdms.keys():
                continue
            else:
                unique_vdms[(cg, aa, vdm_id)] = (vdm_ag, v)

    for (cg, aa, vdm_id) in unique_vdms.keys():
        prefix = cg  + '_' + aa + '_' + str(resnum) + '_' + str(vdm_id)                 
        pr.writePDB(outdir_uni + prefix  + '.pdb.gz', unique_vdms[(cg, aa, vdm_id)][0])


    ### Write summary file.
    df = pd.read_csv(outdir + target.getTitle() + '_summary.tsv', sep = '\t')
    df_group = df.groupby(['file', 'lig'])

    df_score = df_group[['score']].sum()

    df_score.to_csv(outdir + target.getTitle() + '_sum_score.tsv', sep = '\t')

    scores = []
    for g_name, g in df_group:
        df_gg_group = g.groupby(['cg','aa','pos'])
        score = 0
        for gg_name, gg in df_gg_group:
            score += gg[['score']].max()
        scores.append((g_name, score.sum()))

    with open(outdir + + target.getTitle() + '_sum_score_rmdu.tsv', 'w') as f:
        f.write('file\tlig\tscore\n')
        for s in scores:
            f.write('\t'.join([str(x) for x in s[0]]) + '\t' +  str(s[1]) + '\n')

#run_all('o2_1dmm_16-20-28_H-H-D_a_842.pdb')

def main():
    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[1]) -1
    if ind > len(pdb_files) -1:
        return

    run_all(pdb_files[ind])
    return

if __name__=='__main__':
    main()

