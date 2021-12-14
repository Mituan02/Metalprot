import sys
sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
import combs2
import prody as pr
import pandas as pd
import os
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import lil_matrix

class CgCombInfo:
    def __init__(self):

        ligand_id = None

        vdm_cgs = {}


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


def get_vdm_labels_coords(dfa, represent_name, correspond_resname, correspond_names):
    '''
    Example:
    correspond_resname = 'ASP'
    represent_name = 'OD2'
    correspond_names = ['CG', 'OD1', 'OD2']
    '''
    labels = dfa[(dfa['chain'] == 'Y') & (dfa['resname'] == correspond_resname) & (dfa['name'] == represent_name)][['CG', 'rota', 'probe_name', 'seg_chain_resnum']]

    vdm_coords =[]
    for k in correspond_names:
        df_contacts = dfa[
            (dfa['resname'] == correspond_resname) 
            & (dfa['chain'] == 'Y') 
            & (dfa['name'] == k)
        ]
        vdm_coords.append(df_contacts[['c_x', 'c_y', 'c_z']].values.T)
    vdm_coords = np.concatenate(vdm_coords).T

    return labels, vdm_coords


def get_nearest_vdms(vdm_coords, ligand_coords, radius = 2):
    '''
    #to be predecated. 
    '''
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)
    print(adj_matrix.sum())
    print(adj_matrix.shape)

    m_adj_matrix = adj_matrix.tolil()

    all_inds = {}
    for r in range(m_adj_matrix.shape[0]):
        inds = m_adj_matrix.rows[r]
        if len(inds) <= 0:
            continue
        for ind in inds:
            all_inds.add(ind)

    return all_inds


def write_vdms(outdir, all_inds, labels, dfa, prefix):

    for ind in all_inds:
        x = labels.iloc[ind]
        # if x['seg_chain_resnum'][2] != 45:
        #     #print(x['seg_chain_resnum'][2])
        #     continue
        v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
        # if v['resname'].iloc[-1] != 'ARG':
        #     continue
        print(ind)
        combs2.design.functions.print_dataframe(v, outpath=outdir, tag = '_' + str(ind), prefix = prefix)
    return


def get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = 2):
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    dists, inds = nbr.radius_neighbors(ligand_coords)

    return dists, inds


def vdm_ligand_clash(vdm, ligand, clash_radius = 3):
    '''
    The vdm X sidechain may clash with the ligand.
    return True if clash.
    '''
    vdm_coords =[]
    for k in vdm[(vdm['chain'] == 'X')]:
        if k['name'][0] == 'H':
            continue
        vdm_coords.append(k[['c_x', 'c_y', 'c_z']].values)

    ligand_coords = ligand.select('heavy').getCoords()

    nbr = NearestNeighbors(radius=clash_radius).fit(vdm_coords)
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False


def search_vdm(s, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs):
    '''
    s: combs2.design._sample.Sample

    cg_id: (0, 0)  #The first is to record the cg in ligand, the second is to record the cg in vdms. 

    input_dict[cg_id] = {

        cg : 'coo'  #dfa = s.cg_dict['coo']

        lgd_sel : ['C9', 'O3', 'O4']  # correpond to coo: ASP ['CG', 'OD1', 'OD2'] 

        represent_name : 'OD2'

        correspond_resname : 'ASP'

        correspond_names : ['CG', 'OD1', 'OD2']
    }
    '''

    dfa = s.cg_dict[input_dict[cg]]

    ligand_coords = get_ligand_coords(ligands, input_dict[lgd_sel])

    labels, vdm_coords = get_vdm_labels_coords(dfa, input_dict[represent_name], input_dict[correspond_resname], input_dict[correspond_names])

    dists, inds = get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = 2)

    labels_cgs[cg_id] = labels
    df_cgs[cg_id] = dfa
    dist_ind_cgs[cg_id] = (dists, inds)

    return 




def construct_vdm_write(outdir, ligands, labels_cgs, df_cgs, dist_ind_cgs):
    '''
    dist_ind_cgs: dict. {cg: (dists, inds)}, where dists in shape: (len(ligands), )
    df_cgs: {cg: df}
    '''
    CgCombInfoDict = {}

    for i in range(len(ligands)):
        cgCombInfo = CgCombInfo()
        cgCombInfo.ligand_id = i

        pr.writePDB(outdir + str(i) + '_ligand_' + ligand.getTitle(), ligands[i])

        for cg in dist_ind_cgs.keys():
            info = []

            labels = labels_cgs[cg]
            dfa = df_cgs[cg]

            dists = dist_ind_cgs[cg][0][i]
            inds = dist_ind_cgs[cg][1][i]

            for j in range(len(inds)):
                ind = inds[j]
                dist = dists[j]
                x = labels.iloc[ind]
                prefix = str(i) + '_' + '-'.join(str(z) for z in cg) + '_' + str(dist) + '_' + str(x['score']) + '_'              
                v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
                if vdm_ligand_clash(v, ligand):
                    continue
                info.append((dist, v)) 
                combs2.design.functions.print_dataframe(v, outpath=outdir, tag = '_' + str(ind), prefix = prefix)
            cgCombInfo[cg] = info
        CgCombInfoDict[i] = cgCombInfo
    return CgCombInfoDict

            
def write_summary(outdir, CgCombInfoDict, name = '_summary.tsv'):
    with open(outdir + name, 'w') as f:
        f.write('')
        for k in CgCombInfoDict.keys():          
            cgInfo = CgCombInfoDict[k]
            for cg in cgInfo.keys():
                info = cgInfo[cg]
                for o in info:
                    f.write(str(cgInfo.ligand_id) + '\t')
                    f.write(str(cg) + '\t')
                    f.write(str(o[0]) + '\n')
    return