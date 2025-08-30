'''
The direct run and write functions to search/eval ligs or 2ndshell.
'''
import os
import prody as pr
import pandas as pd

from ..basic import constant, utils
from . import gvdm_helper, search_lig_indep, position_ligand, search_lig_indep_2ndshell 

import gc
import shutil


def verify_vdMs(outdir, target, chidres, para, abples, path_to_database, lig, vdm_cg_aa_atommap_dict):
    '''
    For a ligand-binding protein, verify the vdMs used by the protein in a quick way.
    For usage, please check '/mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/pose_lig_by_pair_vdm_dev.py'
    '''
    for cg_id in vdm_cg_aa_atommap_dict.keys():

        pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
        abple = abples[pos.getResindices()[0]]

        for a in vdm_cg_aa_atommap_dict[cg_id]['aas']:
            aa = constant.inv_one_letter_code[a]

            df_vdms = gvdm_helper.load_old_vdm(path_to_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
            df_vdm_filter = gvdm_helper.filter_db(df_vdms, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)   

            ligs = [lig]
            ideal_ala_coords = constant.ideal_ala_coords
            results = search_lig_indep.search_lig_at_cg_aa_resnum(target, chidres, pos, abple, ligs, vdm_cg_aa_atommap_dict, cg_id, df_vdm_filter, ideal_ala_coords, para.rmsd, True, False)
            print(len(results))

            if len(results) <= 0:
                continue

            cg = vdm_cg_aa_atommap_dict[cg_id]['cg']
            for cg_id, lig_id, _rmsd, vdm_ag, vdm_id, score, v, contact_hb, contact_cc in results:
                prefix = 'vdM_Lig-' + str(lig_id) + '_' + cg + '_' + chidres[0] + str(chidres[1]) + '_' + aa + '_rmsd_' + str(round(_rmsd, 2)) + '_v_' + str(round(score, 1)) + '_'       

                vdm_ag.setTitle(prefix + '_' + vdm_ag.getTitle())            
                pr.writePDB(outdir + vdm_ag.getTitle(), vdm_ag)  
    return


def prepare_ligs(outdir, target, task_type = 'search_eval', lig_path = None, para_lig = None):
    '''
    Prepare ligands for vdm search.
    'search_eval': We can extract ligand from the target for the evaluation purpose.
    'search_unknow': We can generate ligand by rotating the ligands.
    'search_exists': We can use pre-generated ligands.
    '''
    if task_type == 'search_unknow':
        position_ligand.run_ligand(outdir, target, lig_path, para_lig.ro1, para_lig.ro2, para_lig.rest1, para_lig.rest2, para_lig.lig_connects, para_lig.geo_sel, para_lig.rot_degree, para_lig.interMolClashSets, clash_dist = para_lig.clash_dist, write_all_ligands=False)
    elif task_type == 'search_eval':
        position_ligand.extract_ligand(outdir, target, para_lig.lig_name)
    elif task_type == 'search_exists':
        for file in os.listdir(lig_path):
            os.makedirs(outdir + 'filtered_ligs/', exist_ok=True)
            if not '.pdb' in file and not '.pdb.gz' in file:
                continue
            shutil.copy2(lig_path + file, outdir + 'filtered_ligs/')

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
    return ligs


def write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict, para):
    ### write filtered vdms.
    summaries = []
    for lig_id in lig_vdm_dict.keys():
        values = lig_vdm_dict[lig_id]
        exist_required_cg_contact = False
        for v in values:
            if v[1][0] == 'phenol' and v[1][-2] == True:
                exist_required_cg_contact = True

        if para.task_type == 'search_unknow' and not exist_required_cg_contact:
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
            (cg, aa, chidres, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc) = v[1]
            if (cg, chidres, aa, vdm_id) in unique_vdms.keys():
                continue
            else:
                unique_vdms[(cg, chidres, aa, vdm_id)] = (vdm_ag, v)

    for (cg, chidres, aa, vdm_id) in unique_vdms.keys():
        prefix = cg  + '_' + chidres[0] + str(chidres[1]) +  '_' + aa + '_' + str(vdm_id)                 
        pr.writePDB(outdir_uni + prefix  + '.pdb.gz', unique_vdms[(cg, chidres, aa, vdm_id)][0])

    ### Write summary file.
    df = pd.read_csv(outdir + target.getTitle() + '_summary.tsv', sep = '\t')
    df_group = df.groupby(['file', 'lig'])

    df_score = df_group[['score']].sum()

    df_score.to_csv(outdir + target.getTitle() +  '_sum_score.tsv', sep = '\t')

    scores = []
    for g_name, g in df_group:
        df_gg_group = g.groupby(['cg','aa','pos'])
        score = 0
        for gg_name, gg in df_gg_group:
            score += gg[['score']].max()
        scores.append((g_name, score.sum()))

    with open(outdir + target.getTitle() + '_' + '_sum_score_rmdu.tsv', 'w') as f:
        f.write('file\tlig\tscore\n')
        for s in scores:
            f.write('\t'.join([str(x) for x in s[0]]) + '\t' +  str(s[1]) + '\n')
    return 


def run_search(target, ligs, path_to_database, para, predefined_chidres, lig_cg_2ndshell = None):
    '''
    cg = 'phenol'
    resnum = 102
    '''
    # load_cg_aa_vdm_dict = {'phenol': ('FW', True, True)}
    # predefined_resnums = [102]
    abples, phipsi = utils.seq_get_ABPLE(target)

    #summaries = []
    lig_vdm_dict = {}
    for cg_id in para.vdm_cg_aa_atommap_dict.keys():
        cg = para.vdm_cg_aa_atommap_dict[cg_id]['cg']
        correspond_resname = para.vdm_cg_aa_atommap_dict[cg_id]['correspond_resname']

        if para.task_type == 'search_2ndshell' and (para.vdm_cg_aa_atommap_dict[cg_id]['correspond_resname'] != lig_cg_2ndshell[0] or cg not in lig_cg_2ndshell[1]):
            print('cg {} is not for the 2ndshell of aa {}.'.format(cg, lig_cg_2ndshell))
            continue

        for a in para.vdm_cg_aa_atommap_dict[cg_id]['aas']:
            aa = constant.inv_one_letter_code[a]
            #df_vdm, df_gr, df_score = load_new_vdm(path_to_database, cg, aa)
            df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)

            for chidres in predefined_chidres:
                print('Searching: ' + cg_id + ' ' + aa + ' ' + chidres[0] + str(chidres[1]))
                pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
                abple = abples[pos.getResindices()[0]]

                if abple == 'n': #>>> The abple can be 'n' at the end of a chid.
                    print('abple is not ABPLE.')
                    continue

                df_vdm_filter = gvdm_helper.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        

                results = search_lig_indep.search_lig_at_cg_aa_resnum(target, chidres, pos, abple, ligs, para.vdm_cg_aa_atommap_dict, 
                    cg_id, df_vdm_filter, constant.ideal_ala_coords, para.rmsd, 
                    para.vdm_cg_aa_atommap_dict[cg_id]['filter_hb'], para.vdm_cg_aa_atommap_dict[cg_id]['filter_cc'])

                for cg_id, lig_id, _rmsd, vdm_ag, vdm_id, score, v, contact_hb, contact_cc in results:
                    prefix = 'Lig-' + str(lig_id) + '_' + cg + '_' + chidres[0] + str(chidres[1]) + '_' + aa + '_rmsd_' + str(round(_rmsd, 2)) + '_v_' + str(round(score, 1)) + '_'       

                    #pr.writePDB(outdir_all + prefix + str(vdm_id) + '.pdb.gz', vdm_ag)
                    #summaries.append((cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score))

                    vdm_ag.setTitle(prefix + str(vdm_id))                       
                    if lig_id in lig_vdm_dict.keys():
                        lig_vdm_dict[lig_id].append((vdm_ag, (cg, aa, chidres, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc)))
                    else:
                        lig_vdm_dict[lig_id] = [(vdm_ag, (cg, aa, chidres, cg_id, lig_id, _rmsd, vdm_id, score, v, contact_hb, contact_cc))]

                    
            del [[df_vdm]]
            gc.collect()
            df_vdm = pd.DataFrame()
        
    return lig_vdm_dict


def run_search_2ndshell(file, workdir, outdir, path_to_database, lig_path, para):

    target_path = workdir + file

    #target, chidres2ind = search_lig_indep.prepare_rosetta_target(outdir, target_path, para.predefined_win_filters)
    target = pr.parsePDB(target_path)
    
    ligs = search_lig_indep_2ndshell.extract_2ndshell_cgs(outdir, target, para.predefined_win_filters)

    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)

    for i in range(len(ligs)):
        lig = ligs[i]
        predefined_chidres = para.predefined_chidres
        if not predefined_chidres and para.task_type == 'search_2ndshell':
            predefined_chidres = search_lig_indep_2ndshell.calc_chidres_around_pose(target, para.predefined_win_filters[i], para.predefined_win_filters, dist = 12)

        lig_vdm_dict = run_search(target, [lig], path_to_database, para, predefined_chidres, (lig.getResnames()[0], para.lig_cg_2ndshell[i]))
        print('lig_vdm_dict size {}'.format(len(list(lig_vdm_dict.keys()))))
        write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict, para)

    return