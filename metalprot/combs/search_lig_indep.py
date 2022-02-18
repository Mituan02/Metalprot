import os
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant
from metalprot.basic import utils
from metalprot.basic import prody_ext
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix
import gc
import shutil


def prepare_rosetta_target(outdir, target_path, designed_site = [15, 19, 27]):
    '''
    Or searching rosetta 1st round designed target backbone. 
    All aas are mutated into gly except the metal binding ones.
    '''
    target = pr.parsePDB(target_path)

    resindices, target_index_dict = prody_ext.transfer2resindices(target, designed_site)

    pdb_gly = prody_ext.target_to_all_gly_ala(target, target.getTitle() + '_gly.pdb', resindices, aa= 'GLY', keep_no_protein = True )
    
    pr.writePDB(outdir + pdb_gly.getTitle() + '.pdb', pdb_gly)

    return pdb_gly


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


def search_lig_at_cg_aa_resnum(target, resnum, pos, abple, ligs, input_dict, cg, df_vdm, ideal_ala_coords, rmsd = 0.7, filter_hb = False, filter_cc = False):
    
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

                    if (not filter_hb and not filter_cc) or (filter_hb and v[['contact_hb']].any()[0]) or (filter_cc and v[['contact_cc']].any()[0]):
                        
                        ag = df2ag(v)
                        tf_rev.apply(ag)

                        if clash_filter_protein_single(target, resnum, ag) or clash_filter_lig(ligs[u], ag):
                            #print('filtered')
                            continue
                        results.append((cg_id, u, rmsd, ag, ind, v['C_score_ABPLE_' + abple].values[0], None, v[['contact_hb']].any()[0], v[['contact_cc']].any()[0]))
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