import os
import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors

from ..basic.vdmer import get_contact_atom
from ..basic import cluster
from ..basic import transformation
from ..basic import constant

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


def construct_pseudo_2ndshellVdm(target, vdM, w, wins):
    '''
    Find the resind of the vdm on target. Then extract resinds of atoms within a distance. 
    Followed by extracting the vdm resind and the atoms resind pairs together with the metal. 
    '''
        
    nearby_aas = target.select('protein and not carbon and not hydrogen and within 10 of resindex ' + str(w))
    nearby_aa_resinds = [x for x in np.unique(nearby_aas.getResindices()) if x not in wins]

    ags = []
    count = 0
    for resind in nearby_aa_resinds:
        if vdM.win and resind in vdM.win:
            continue
        neary_aas_coords = []
        neary_aas_coords.extend(target.select('name N C CA O and resindex ' + str(resind)).getCoords())
        neary_aas_coords.extend(vdM.query.select('name N C CA O and resindex ' + str(vdM.contact_resind)).getCoords())
        neary_aas_coords.extend(vdM.query.select(metal_sel).getCoords())
        coords = np.array(neary_aas_coords)

        names = []
        names.extend(target.select('name N C CA O and resindex ' + str(resind)).getNames())
        names.extend(vdM.query.select('name N C CA O and resindex ' + str(vdM.contact_resind)).getNames())
        names.extend(vdM.query.select(metal_sel).getNames())
        
        resnums = [0, 0, 0, 0, 1, 1, 1, 1, 2]    

        atom_contact_pdb = pr.AtomGroup('nearby_bb' + str(count))
        atom_contact_pdb.setCoords(coords)
        atom_contact_pdb.setNames(names)
        atom_contact_pdb.setResnums(resnums)
        ags.append(atom_contact_pdb)
        count +=1

    return ags


def get_2ndvdmCoords_rot(ag, secondshell_vdm_coords):
    '''
    The 'transformation' way. Note the method is not used or tested. It is a backup method. 
    '''
    secondshell_vdm_coords_t = []
    agCoords = ag.select('resindex 1 and name N C CA').getCoords()
    R, m_com, t_com = transformation.get_rot_trans(secondshell_vdm_coords[0], agCoords)
    for i in range(len(secondshell_vdm_coords)):
        secondshell_vdm_coords_t.append(np.dot((secondshell_vdm_coords - m_com), R) + t_com)
    return secondshell_vdm_coords_t


def get_2ndvdmCoords_prody(ag, root_2ndvdm_query, allInOne_2ndvdm, secondshell_vdms_count):
    '''
    The prody way.
    '''
    _allInOne_2ndvdm = allInOne_2ndvdm.copy()
    transform = pr.calcTransformation(root_2ndvdm_query.select('resindex 1 and name N C CA') ,ag.select('resindex 1 and name N C CA'))
    transform.apply(_allInOne_2ndvdm)
    all2ndCoords = _allInOne_2ndvdm.getCoords().reshape(secondshell_vdms_count, 27)
    return all2ndCoords, transform


def extract_candidates(ags, secondshell_vdms, candidateIds, transform):
    '''
    
    '''
    candidates = []
    if len(candidateIds) <= 0:
        return candidates
    for x, y in candidateIds:
        ag = ags[x]
        vdm = secondshell_vdms[y].copy()
        transform.apply(vdm.query)
        candidates.append(vdm)
    return candidates


def search_2ndshell(comb_dict, key, target, secondshell_vdms, secondshell_vdm_aatype, allInOne_2ndvdm, rmsd_2ndshell):
    '''
    Construct all possible 2nd shell from target and the 1st shell vdm candidate. 
    Then transform the secondashell_coords to the 1st shell vdm's 'N CA C' atoms. 
    Following using NearestNeighbor to get radius_neighbors_graph. 
    '''
    for w in key[0]:
        vdm = comb_dict[key].centroid_dict[w]

        ags = construct_pseudo_2ndshellVdm(target, vdm, w, key[0])
        ags_coords = [ag.getCoords().flatten() for ag in ags]         
        
        _2ndvdm_coords, transform = get_2ndvdmCoords_prody(ags[0], secondshell_vdms[0].query, allInOne_2ndvdm, len(secondshell_vdms))

        radius = np.sqrt(len(ags[0].getCoords())) * rmsd_2ndshell

        nbr = NearestNeighbors(radius=radius).fit(_2ndvdm_coords)
        #dists, inds = nbr.radius_neighbors(_2ndvdm_coords)
        adj_matrix = nbr.radius_neighbors_graph(ags_coords).astype(bool)

        mask_aa_type = secondshell_vdm_aatype == vdm.aa_type

        candidateIds = []
        for r in range(adj_matrix.shape[0]):
            inds = np.where(adj_matrix.getrow(r).toarray())[1]
            for c in inds:
                if mask_aa_type[c]:
                    candidateIds.append((r, c))

        candidates = extract_candidates(ags, secondshell_vdms, candidateIds, transform)
        comb_dict[key].secondshell_dict[w] = candidates
    return

'''
Inheritated from Search_vdM. 
Search the 2nd shell h-hond. 
'''

def run_search_2ndshell(comb_dict, target, secondshell_vdms, allInOne_2ndvdm, rmsd_2ndshell):
    '''
    
    '''
    print('run search 2nd-shell vdM.')
    if not comb_dict:
        print('No 1st shell metal-binding vdM found. No need to search 2ndshell.')

    secondshell_vdm_aatype = np.array([constant.one_letter_code[v.query.select('name C and resindex 1').getResnames()[0]] for v in secondshell_vdms])

    for key in comb_dict.keys():
        search_2ndshell(comb_dict, key, target, secondshell_vdms, secondshell_vdm_aatype, allInOne_2ndvdm, rmsd_2ndshell)

    return 


def write_2ndshell(ss, workdir, comb_dict):
    '''
    #Could be combined in write_comb_info.
    '''
    for key in comb_dict.keys():
        outdir = workdir + 'win_' + '-'.join([ss.target_index_dict[w] for w in key[0]]) + '/'
        
        tag = 'W_' + '-'.join([ss.target_index_dict[w] for w in key[0]]) + '_X_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1])
        
        outdir += tag + '_hb/'
        os.makedirs(outdir, exist_ok=True)

        for w in key[0]:
            for c in comb_dict[key].secondshell_dict[w]:
                c_names = c.query.getTitle().split('_')
                out_file = outdir + tag + '_w_' + str(w) + '_hb_' + '_'.join([c_names[i] for i in range(5)])
                pr.writePDB(out_file, c.query)

    with open(outdir + tag + '_2ndshell_summary.tsv', 'w') as f:
        f.write('name\twin\trmsd\n')
        for key in comb_dict.keys():
            for w in key[0]:
                for c in comb_dict[key].secondshell_dict[w]:
                    name = tag + '_w_' + str(w) + '_hb_' + c[0].query.getTitle()
                    f.write(name + '\t' + str(w) + '\t' + 'None' + '\t' + str(c.score) + '\n')
    
    return 