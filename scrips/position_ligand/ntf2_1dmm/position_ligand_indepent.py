'''
The script is to design ligand-metal binding enzyme using the vdM idea. 
Specially, here is to search the combs vdM library without using Combs2.
'''

import os
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant
from metalprot.basic import utils
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import gc
import shutil

with open('/mnt/e/GitHub_Design/Combs2/combs2/files/ideal_alanine_bb_only.pkl', 'rb') as f:
    ideal_alanine_bb_only = pickle.load(f)
ideal_ala_coords = np.array(ideal_alanine_bb_only[['c_x', 'c_y', 'c_z']])


def load_new_vdm(path_to_database, cg, aa):
    '''
    The function is used to load vdM database after 2022.02.
    Pelase check combs2_da.ipynb to learn the loaded database.
    '''
    cg_aa = cg + '/' + aa 
    df_vdm = pd.read_parquet(path_to_database + 'vdMs/' + cg_aa + + '.parquet.gzip')

    with open(path_to_database + 'vdMs_gr_indices/' + cg_aa + '.pkl', 'rb') as f:
        df_gr = pickle.load(f)

    df_score = pd.read_parquet(path_to_database + 'nbrs/vdMs_cg_nbrs_scores/' + cg_aa + '.parquet.gzip')

    return df_vdm, df_gr, df_score


def load_old_vdm(path_to_database, cg, aa):
    '''
    The function is used to load vdM database before 2022.02.
    Pelase check combs2_da.ipynb to learn the loaded database.
    '''
    cg_aa = cg + '/' + aa 
    df_vdm = pd.read_parquet(path_to_database + 'vdMs/' + cg_aa + '.parquet.gzip')

    return df_vdm


def ligscoords_2_ideal_ala(pos, ligs, ideal_ala_coords):
    '''
    For all the pre generated ligs, combine them to each position in bb and then transform the pos_lig to 'ideal alanine'.
    The 'ideal alanine' is the alanine used for the vdM database.
    TO DO: this function can be more efficient to calc the matrix instead of calc each ligs.
    '''
    pos_coords = pos.select('name N CA C').getCoords()
    tf = pr.calcTransformation(pos_coords, ideal_ala_coords)
    _ligs = []
    for lig in ligs:
        _lig = lig.copy()
        tf.apply(_lig)
        _ligs.append(_lig)
    
    tf_rev = pr.calcTransformation(ideal_ala_coords, pos_coords)
    return _ligs, tf_rev



def get_ligand_coords(filtered_ligands, lgd_sel):
    # Get all ligand coords.
    ligand_coords = []
    for i in range(len(filtered_ligands)):
        lgd = filtered_ligands[i]
        coords = []
        for lsa in lgd_sel:
            coords.extend(lgd.select('name ' + lsa).getCoords().flatten())
        ligand_coords.append(coords)
    ligand_coords = np.array(ligand_coords)

    return ligand_coords
    

def get_vdm_labels_coords_4old_vdm_db(dfa, represent_name, correspond_resname, correspond_names):
    '''
    Example:
    correspond_resname = 'ASP'
    represent_name = 'OD2'
    correspond_names = ['CG', 'OD1', 'OD2']
    '''
    labels = dfa[(dfa['chain'] == 'Y') & (dfa['resname'] == correspond_resname) & (dfa['name'] == represent_name)][['CG', 'rota', 'probe_name']]

    vdm_coords =[]
    for k in correspond_names:
        df_contacts = dfa[
            (dfa['resname'] == correspond_resname) 
            & (dfa['chain'] == 'Y') 
            & (dfa['name'] == k)
        ]
        vdm_coords.append(df_contacts[['c_x', 'c_y', 'c_z']].values.T)
    try:
        vdm_coords = np.concatenate(vdm_coords).T
    except:
        vdm_coords = np.array([])
        print(dfa.head(50))
        print(represent_name)
        print(correspond_resname)
        print(correspond_names)
        print('get_vdm_labels_coords. labels length is ' + str(len(labels)))

    return labels, vdm_coords


def get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = 1):
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    dists, inds = nbr.radius_neighbors(ligand_coords)

    return dists, inds


def filter_db(df_vdm, use_enriched = True, use_abple = True, abple = 'A'):
    '''
    Filter database based on ABPLE, enriched. Check combs2._sample._load_res().
    '''
    if use_enriched and use_abple:
        return df_vdm[df_vdm['C_score_' + 'ABPLE_' + abple] > 0]
    if use_enriched and not use_abple:
        return df_vdm[df_vdm['C_score_bb_ind'] > 0]
    return df_vdm


def search_lig_at_cg_aa_resnum(target, resnum, pos, abple, ligs, input_dict, cg, df_vdm, rmsd = 0.7):
    
    results = []

    _ligs, tf_rev = ligscoords_2_ideal_ala(pos, ligs, ideal_ala_coords)

    for cg_id in input_dict.keys():
        if input_dict[cg_id]['cg'] != cg:
            continue
        ligand_coords = get_ligand_coords(_ligs, input_dict[cg_id]['lgd_sel'])
        
        labels, vdm_coords = get_vdm_labels_coords_4old_vdm_db(df_vdm, input_dict[cg_id]['represent_name'], input_dict[cg_id]['correspond_resname'], input_dict[cg_id]['correspond_names'])

        if vdm_coords.shape[0] <=0 or vdm_coords.shape[0] != labels.shape[0]:
            print('cg_id not working: {}'.format(cg_id))
            continue
        
        num_cg_atoms = len(input_dict[cg_id]['lgd_sel'])
        radius = np.sqrt(num_cg_atoms) * rmsd

        dists, inds = get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = radius)

        #u is lig id, v is vdM id.
        for u in range(len(inds)): 
            for v in range(len(inds[u])):                
                ind = inds[u][v]
                rmsd = dists[u][v]/np.sqrt(num_cg_atoms)
                try:
                    x = labels.iloc[ind]
                    v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
                    ag = df2ag(v)
                    tf_rev.apply(ag)

                    if clash_filter_protein_single(target, resnum, ag) or clash_filter_lig(ligs[u], ag):
                        #print('filtered')
                        continue
                    results.append((cg_id, u, rmsd, ag, ind, v['C_score_ABPLE_' + abple].values[0]))
                except:
                    print(cg_id)
                    print(x)

    return results


def df2ag(df, b_factor_column=None):
    df = df.copy()
    if 'chain' not in df.columns:
        df['chain'] = 'A'
    ag = pr.AtomGroup()
    ag.setCoords(df[['c_x','c_y','c_z']].values)
    ag.setResnums(df['resnum'].values)
    ag.setResnames(df['resname'].values)
    ag.setNames(df['name'].values)
    ag.setChids(df['chain'].values)
    ag.setSegnames(df['chain'].values)
    heteroflags = ag.getSegnames() == 'L'
    ag.setFlags('hetatm', heteroflags)
    if 'beta' in df.columns and b_factor_column is None:
        ag.setBetas(df['beta'].values)
    elif b_factor_column is not None:
        ag.setBetas(df[b_factor_column].values)
    if 'occ' not in df.columns:   
        df['occ'] = 1
    #pr.writePDB(outpath + prefix + filename + tag + '.pdb.gz', ag, occupancy=df['occ'].values)
    return ag


def clash_filter_protein(target, all_coords, dist_cut = 2.5):
    '''
    # clashing: any atom close to the bb let's say dist <= dist_cut

    The bb here should contain 'N CA C O of any aa, the sc of binding aa, Metal'
    '''
    
    bb_coords = target.select('heavy').getCoords()
    nn = NearestNeighbors(radius= dist_cut).fit(bb_coords)
    x_in_y = nn.radius_neighbors_graph(all_coords).toarray().astype(bool)
    clashing = np.any(x_in_y, axis=1)

    #np.sum(clashing) #The clashing is too less to be considered.
    return clashing


def clash_filter_protein_single(target, resnum, vdm, dist_cut = 2.7):
    '''
    # clashing: any atom close to the bb let's say dist <= dist_cut

    The bb here should contain 'N CA C O of any aa, the sc of binding aa, Metal'
    '''
    target_coords = target.select('heavy and not resnum ' + str(resnum)).getCoords()
    vdm_sel = vdm.select('heavy and sc and chain X and resnum 10') #For Gly, it can be None.
    if not vdm_sel:
        return True
    vdm_coords = vdm_sel.getCoords()
    nbrs = NearestNeighbors(radius= dist_cut).fit(target_coords)
    adj_matrix = nbrs.radius_neighbors_graph(vdm_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False


def clash_filter_lig(lig, vdm, dist_cut = 2.5):
    '''
    Check if the vdm clashing the ligand
    '''
    lig_coords = lig.select('heavy').getCoords()
    vdm_sel = vdm.select('heavy and sc and chain X and resnum 10') #For Gly, it can be None.
    if not vdm_sel:
        return True
    vdm_coords = vdm_sel.getCoords()
    nbrs = NearestNeighbors(radius= dist_cut).fit(lig_coords)
    adj_matrix = nbrs.radius_neighbors_graph(vdm_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False

def run():
    '''
    cg = 'phenol'
    resnum = 102
    '''
    # load_cg_aa_vdm_dict = {'phenol': 'F'}
    # predefined_resnums = [102]

    summaries = []
    for cg in load_cg_aa_vdm_dict.keys():

        for a in load_cg_aa_vdm_dict[cg]:
            aa = constant.inv_one_letter_code[a]
            #df_vdm, df_gr, df_score = load_new_vdm(path_to_database, cg, aa)
            df_vdm = load_old_vdm(path_to_database, cg, aa)

            for resnum in predefined_resnums:
                print('Searching: ' + cg + ' ' + aa + ' ' + str(resnum))
                pos = target.select('resnum ' + str(resnum))
                abple = abples[pos.getResindices()[0]]
                df_vdm_filter = filter_db(df_vdm, abple=abple)        

                results = search_lig_at_cg_aa_resnum(target, resnum, pos, abple, ligs, input_dict, cg, df_vdm_filter, rmsd)

                pdbs = {}
                for cg_id, lig_id, _rmsd, vdm_ag, vdm_id, score in results:
                    prefix = 'Lig-' + str(lig_id) + '_' + cg  + '_' + aa + '_' + str(resnum) + '_rmsd_' + str(round(_rmsd, 2)) + '_v_' + str(round(score, 1)) + '_'       
                    pr.writePDB(outdir_all + prefix + str(vdm_id) + '.pdb.gz', vdm_ag)

                    summaries.append((cg, aa, resnum, cg_id, lig_id, _rmsd, vdm_id, score))
                    if (cg, aa, vdm_id) in pdbs.keys():
                        continue
                    else:
                        pdbs[(cg, aa, vdm_id)] = vdm_ag

                for (cg, aa, vdm_id) in pdbs.keys():
                    prefix = cg  + '_' + aa + '_' + str(resnum) + '_'                   
                    #print_dataframe(vdm, filename = str(vdm_id), outpath=outdir, prefix = prefix)
                    pr.writePDB(outdir + prefix + str(vdm_id) + '.pdb.gz', pdbs[(cg, aa, vdm_id)])
                    
            del [[df_vdm]]
            gc.collect()
            df_vdm = pd.DataFrame()

    
    with open(outdir + '_summary.tsv', 'w') as f:
        for s in summaries:
            #print(s)
            f.write('\t'.join([str(x) for x in s]) + '\n')
        
    return results

#path_to_database='/mnt/e/DesignData/Combs_update_20220210/'
path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'

load_cg_aa_vdm_dict = {
    'coo': 'AGKNQRSTWY',
    'bb_cco': 'AGKNQRSTWY',
    'phenol': 'ADEFGKNQRSTWY'
}

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/output_selfcenter_o1_1dmm_16-20-28_H-H-D_a_820__20220210-101859/represents_combs/'

outdir = workdir + 'W_15-19-27_H-H-D_1000-404-467_/vdms_output/'
os.makedirs(outdir, exist_ok=True)

outdir_all = workdir + 'W_15-19-27_H-H-D_1000-404-467_/vdms_output_all/'
os.makedirs(outdir_all, exist_ok=True)

target = pr.parsePDB(workdir + 'W_15-19-27_H-H-D_1000-404-467_allgly.pdb')
abples, phipsi = utils.seq_get_ABPLE(target)

predefined_resnums = [11, 12, 15, 16, 19, 24, 27, 30, 31, 35, 37, 39, 46, 52, 55, 56, 60, 65, 67, 69, 83, 85, 87, 98, 100, 102, 104, 106, 112, 115, 117, 119]

# outdir_ligs = workdir + 'W_15-19-27_H-H-D_1000-404-467_/ligs_inorder/'
# os.makedirs(outdir_ligs, exist_ok=True)
ligs = []
lig_id = 0
for file in os.listdir(workdir + 'W_15-19-27_H-H-D_1000-404-467_/filtered_ligs/'):
    if '.pdb' in file:
        lig = pr.parsePDB(workdir + 'W_15-19-27_H-H-D_1000-404-467_/filtered_ligs/' + file)
        # shutil.copy2(workdir + 'W_15-19-27_H-H-D_1000-404-467_/filtered_ligs/' + file, outdir_ligs + 'Lig-'+ str(lig_id) + '_'+ file)
        # lig_id += 1
        ligs.append(lig)

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
    }
}

rmsd = 0.6

run()

'''
# Testing clashing. loading files.
import prody as pr
from sklearn.neighbors import NearestNeighbors

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/output_selfcenter_o1_1dmm_16-20-28_H-H-D_a_820__20220210-101859/represents_combs/'

target = pr.parsePDB(workdir + 'W_15-19-27_H-H-D_1000-404-467_allgly.pdb')

outdir = workdir + 'W_15-19-27_H-H-D_1000-404-467_/vdms_output/'

vdm = pr.parsePDB(outdir + 'phenol_PHE_102_2763.pdb')

outdir_ligs = workdir + 'W_15-19-27_H-H-D_1000-404-467_/ligs_inorder/'

lig = pr.parsePDB(outdir_ligs + 'Lig-150_Geo_5_tts_fe_C8-C7_255_C7-C6_355.pdb')


'''