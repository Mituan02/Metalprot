import prody as pr
import numpy as np
from prody.utilities.catchall import getCoords
from scipy.spatial.transform import Rotation
from sklearn.neighbors import NearestNeighbors

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToPDBFile

from ..basic import prody_ext 

def rotate_ligs(orign_lig, rot, rest, rotation_degree = 10, dist = 3):
    '''
    Note that the ZN must be the last element of the ligand.
    '''
    lig = orign_lig.copy()

    # transform the lig to z axis.
    rot_coords = lig.select('name ' + ' '.join(rot)).getCoords()
    rot_dist = pr.calcDistance(rot_coords[1], rot_coords[0])

    z_coords = np.zeros((2, 3), dtype=float)
    z_coords[1, -1] = rot_dist

    pr.calcTransformation(rot_coords, z_coords).apply(lig)

    all_ligs = []
    for i in range(0, 360, rotation_degree):
        _lig = lig.copy()
        rotation = Rotation.from_rotvec(np.radians(i)*np.array([0, 0, 1]))        
        _coords = rotation.apply(_lig.getCoords())
        _lig.setCoords(_coords)
        if len(rest) > 0:
            _lig.select('name ' + ' '.join(rest)).setCoords(lig.select('name ' + ' '.join(rest)).getCoords())

        pr.calcTransformation(_lig.select('name ' + ' '.join(rest)).getCoords(), orign_lig.select('name ' + ' '.join(rest)).getCoords()).apply(_lig)
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)

        if ligand_rot_is_clash(_lig, rot, rest, dist):
            continue
        all_ligs.append(_lig)

    return all_ligs


def generate_rotated_ligs(lig, rots, rests, rotation_degrees, clash_dist = 3):
    '''
    The method only works for 2 rots.
    TO DO: use recursive algorithm to allow more than 2 rots.
    '''
    all_ligs = []
    i = 0
    rot = rots[i]
    rest = rests[i]
    degree = rotation_degrees[i]
    ligs = rotate_ligs(lig, rot, rest, degree, clash_dist)
    for _lig in ligs:
        i = 1
        rot = rots[i]
        rest = rests[i]
        degree = rotation_degrees[i]
        ligs2 = rotate_ligs(_lig, rot, rest, degree, clash_dist)
        all_ligs.extend(ligs2)
    return all_ligs


def generate_rotated_ligs_rdkit(lig_smiles, total, outdir):
    m = Chem.MolFromSmiles(lig_smiles)
    m2=Chem.AddHs(m)
    AllChem.EmbedMolecule(m2)

    cids = AllChem.EmbedMultipleConfs(m2, numConfs=total)
    for i in range(total):
        MolToPDBFile(m2, outdir + 'rdkit_' + str(i) + '.pdb', confId = i)

    m2s = []
    for p in os.listdir(outdir):
        if '.pdb' not in p:
            continue
        m2s.append(pr.parsePDB(outdir + p))
    return m2s


def add_metal2lig(lig, rig, lig_sel, rig_sel, metal):
    '''
    Add metal to ligs that do not contain metal for superimpose on binding geometry.
    The rdkit generated ligands generally do not contain metal.

    rig: the ligand-metal prody obejct extracted from pdb.
    '''
    prody_ext.ordered_sel_transformation(rig, lig, rig_sel, lig_sel)
    
    metal_point = rig.select('name ' + metal)[0].getCoords()
    points = [x.getCoords() for x in lig.select('heavy')]
    points.append(metal_point)

    names = [x.getName() for x in lig.select('heavy')]
    names.append(metal)
    resnames = [x.getResname() for x in lig.select('heavy')]
    resnames.append(metal)
    resnums = [0 for x in lig.select('heavy')]
    resnums.append(1)
    mm = prody_ext.transfer2pdb(points, names, resnums, resnames, title= lig.getTitle())
    return mm


def add_metal2ligs(ligs, rig, lig_sel, rig_sel, metal):
    lig_metals = []
    for lig in ligs:
        lig_metal = add_metal2lig(lig, rig, lig_sel, rig_sel, metal)
        lig_metals.append(lig_metal)
    return lig_metals

def ligand_rot_is_clash(lig, rot, rest, dist = 3):
    '''
    The ligand it self could clash after rotation.
    The idea is to sep the ligand by rot into 2 groups anc check their dists.

    return True if clash
    '''
    atoms_rot = lig.select('heavy and not name ' + ' '.join(rot) + ' ' + ' '.join(rest)).getCoords()
    rest_coords = lig.select('heavy and name ' +  ' '.join(rest)).getCoords()

    nbrs = NearestNeighbors(radius= dist).fit(atoms_rot)
    adj_matrix = nbrs.radius_neighbors_graph(rest_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False


def lig_2_ideageo(ligs, lig_connect_sel, ideal_geo_o = None, geo_sel = 'OE2 ZN'):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''
    _lig = ligs[0]

    mobile_sel_coords = []
    for s in lig_connect_sel:
        mobile_sel_coords.append(_lig.select('name ' + s).getCoords()[0])

    ideal_geo_sel_coords = []
    for s in geo_sel.split(' '):
        ideal_geo_sel_coords.append(ideal_geo_o.select('name ' + s).getCoords()[0])

    transformation = pr.calcTransformation(np.array(mobile_sel_coords), np.array(ideal_geo_sel_coords))

    for lg in ligs:
        transformation.apply(lg)

    return


def ligand_clashing_filter(ligs, target, dist = 3):
    '''
    The ligand clashing: the ligs cannot have any heavy atom within 3 A of a target bb.
    Nearest neighbor is used to calc the distances. 
    '''
    all_coords = []
    labels = []

    for i in range(len(ligs)):
        coords = ligs[i].select('heavy').getCoords()
        all_coords.extend(coords)
        labels.extend([i for j in range(len(coords))])

    
    target_coords = target.select('name N C CA O CB').getCoords()

    nbrs = NearestNeighbors(radius= dist).fit(target_coords)
    adj_matrix = nbrs.radius_neighbors_graph(all_coords).astype(bool)

    failed = set()
    for i in range(adj_matrix.shape[0]):
        if adj_matrix.getrow(i).toarray().any():
            failed.add(labels[i])
    
    filtered_ligs = []
    for i in range(len(ligs)):
        if i in failed:
            continue
        filtered_ligs.append(ligs[i])

    return filtered_ligs
    


