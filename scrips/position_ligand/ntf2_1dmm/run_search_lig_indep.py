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
from metalprot.basic import constant, utils
from metalprot.combs import search_lig_indep
from metalprot.combs import position_ligand
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import gc
import shutil
import datetime
import importlib.machinery

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_lig_indep.py 0 /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_lig_paras_rdkit.py 7

python run_search_lig_indep.py 0 search_lig_paras_eval.py 7
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_search_lig_indep.py 0 /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_2ndshell_paras.py 5
'''

def extract_2ndshell_cgs(outdir, target, win_filters):
    cgs = []
    for w in win_filters:
        cg = target.select('protein and heavy and sc and chid ' + w[0] + ' and resnum ' + str(w[1])).toAtomGroup()
        cg.setTitle('w_' + w[0] + str(w[1]))
        cgs.append([cg])
        pr.writePDB(outdir + cg.getTitle() + '.pdb', cg)
    return cgs


def run_search(target, ligs, path_to_database, para, lig_cg):
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
        if para.task_type == 'search_2ndshell' and cg not in lig_cg:
            print('cg {} is not for the 2ndshell of aa {}.'.format(cg, lig_cg))
            continue

        for a in para.vdm_cg_aa_atommap_dict[cg_id]['aas']:
            aa = constant.inv_one_letter_code[a]
            #df_vdm, df_gr, df_score = load_new_vdm(path_to_database, cg, aa)
            df_vdm = search_lig_indep.load_old_vdm(path_to_database, cg, aa)

            for chidres in para.predefined_resnums:
                print('Searching: ' + cg + ' ' + aa + ' ' + chidres[0] + str(chidres[1]))
                pos = target.select('chid ' + chidres[0] + ' and resnum ' + str(chidres[1]))
                abple = abples[pos.getResindices()[0]]
                df_vdm_filter = search_lig_indep.filter_db(df_vdm, use_enriched = para.use_enriched, use_abple=para.use_abple, abple=abple)        

                results = search_lig_indep.search_lig_at_cg_aa_resnum(target, chidres, pos, abple, ligs, para.vdm_cg_aa_atommap_dict, cg_id, df_vdm_filter, constant.ideal_ala_coords, para.rmsd, cg_aa_vdms[1], cg_aa_vdms[2])

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


### Run Search ligands.

def prepare_ligs(outdir, target, lig_path, para):
    if para.task_type == 'search_unknow':
        position_ligand.run_ligand(outdir, target, lig_path, para.ro1, para.ro2, para.rest1, para.rest2, para.lig_connects, para.geo_sel, para.rot_degree, para.interMolClashSets, clash_dist = 2.7, write_all_ligands=False)
    else:
        position_ligand.extract_ligand(outdir, target, para.lig_name)

    outdir_ligs = outdir + 'ligs_inorder/'
    os.makedirs(outdir_ligs, exist_ok=True)
    ligss = []
    ligs = []
    lig_id = 0
    for file in os.listdir(outdir + 'filtered_ligs/'):
        if not '.pdb' in file and not '.pdb.gz' in file:
            continue
        lig = pr.parsePDB(outdir + 'filtered_ligs/' + file)
        shutil.copy2(outdir + 'filtered_ligs/' + file, outdir_ligs + 'Lig-'+ str(lig_id) + '_'+ file)
        lig_id += 1
        ligs.append(lig)
    ligss.append(ligs)
    return ligss


def write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict, i, para):
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

    with open(outdir + target.getTitle() + '_' + str(i) + '_summary.tsv', 'w') as f:
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
    df = pd.read_csv(outdir + target.getTitle() + '_' + str(i) + '_summary.tsv', sep = '\t')
    df_group = df.groupby(['file', 'lig'])

    df_score = df_group[['score']].sum()

    df_score.to_csv(outdir + target.getTitle() + '_' + str(i) + '_sum_score.tsv', sep = '\t')

    scores = []
    for g_name, g in df_group:
        df_gg_group = g.groupby(['cg','aa','pos'])
        score = 0
        for gg_name, gg in df_gg_group:
            score += gg[['score']].max()
        scores.append((g_name, score.sum()))

    with open(outdir + target.getTitle() + '_' + str(i) +  '_sum_score_rmdu.tsv', 'w') as f:
        f.write('file\tlig\tscore\n')
        for s in scores:
            f.write('\t'.join([str(x) for x in s[0]]) + '\t' +  str(s[1]) + '\n')
    return 



def run_all(file, workdir, path_to_database, lig_path, para):
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 

    target_path = workdir + file
    outdir = workdir + para.task_type  + '_result/output_' + '_'+ file + '_' + time_tag + '/'
    os.makedirs(outdir, exist_ok=True)

    target, chidres2ind = search_lig_indep.prepare_rosetta_target(outdir, target_path, para.predefined_win_filters)
    if para.task_type == 'search_2ndshell':
        ligs = extract_2ndshell_cgs(outdir, target, para.predefined_win_filters)
    else:
        ligs = prepare_ligs(outdir, target, lig_path, para)

    print('Filtered Ligs: ' + str(len(ligs)))
    if len(ligs) <= 0:
        return

    outdir_uni = outdir + 'vdms_output_uni/'
    os.makedirs(outdir_uni, exist_ok=True)

    outdir_all = outdir + 'vdms_output_all/'
    os.makedirs(outdir_all, exist_ok=True)
    print('number of ligs: {}'.format(len(ligs)))
    for i in range(len(ligs)):
        lig = ligs[i]
        lig_vdm_dict = run_search(target, lig, path_to_database, constant.ideal_ala_coords, para, para.lig_cgs[i])
        print('lig_vdm_dict size {}'.format(len(list(lig_vdm_dict.keys()))))
        write_vdm(outdir, outdir_all, outdir_uni, target, lig_vdm_dict, i, para)

    return


def main():
    #path = '/mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/search_lig_paras.py'
    on_wynton = bool(int(sys.argv[1]))
    path = sys.argv[2]
    para = importlib.machinery.SourceFileLoader('para', path).load_module()
    print(path)
    print('Task: ' + para.task_type)
    workdir, path_to_database, lig_path = para.get_file_path(on_wynton)
    print('on_wynton: ' + str(on_wynton))

    pdb_files = sorted([fp for fp in os.listdir(workdir) if fp[0] != '.' and '.pdb' in fp])

    ind = int(sys.argv[3]) -1
    if ind > len(pdb_files) -1:
        return
    print(pdb_files[ind])
    
    run_all(pdb_files[ind], workdir,  path_to_database, lig_path, para)
    return


if __name__=='__main__':
    main()

