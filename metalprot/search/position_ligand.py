import prody as pr
import numpy as np
from scipy.spatial.transform import Rotation
from sklearn.neighbors import NearestNeighbors


def rotate_ligs_first(lig, rot, rotation_degree = 5, metal = 'ZN'):
    '''
    Note that the ZN must be the last element of the ligand.
    '''
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
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)
        all_ligs.append(_lig)

    return all_ligs


def rotate_ligs_second(lig, rot, base_rot, rotation_degree = 5, metal = 'ZN'):
    '''
    Note that the ZN must be the last element of the ligand.
    '''
    base_coords = lig.select('name ' + ' '.join(base_rot) + ' ' + ' '.join(rot)).getCoords()

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
        _lig.select('name ' + metal).setCoords(lig.select('name ' + metal).getCoords())
        base_rot_coords = _lig.select('name ' + ' '.join(base_rot) + ' ' + ' '.join(rot)).getCoords()
        pr.calcTransformation(base_rot_coords, base_coords).apply(_lig)
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)
        all_ligs.append(_lig)

    return all_ligs


def generate_rotated_ligs(lig, ro1, ro2, rotation_degree = 10, metal = 'ZN'):

    all_ligs = []
    
    rot1_ligs = rotate_ligs_first(lig, ro1, rotation_degree=rotation_degree, metal= metal)
    for lig1 in rot1_ligs:
        _ligs = rotate_ligs_second(lig1, ro2, ro1, rotation_degree=rotation_degree, metal= metal)
        all_ligs.extend(_ligs)

    return all_ligs


def lig_2_ideageo(ligs, lig_connects, ideal_geo_o = None, geo_sel = 'name OE2 ZN', metal_sel = 'name NI MN ZN CO CU MG FE'):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''
    _lig = ligs[0]

    lig_sel = _lig.select('name ' + ' '.join(lig_connects))

    ideal_geo_sel = ideal_geo_o.select(geo_sel)

    transformation = pr.calcTransformation(lig_sel, ideal_geo_sel)

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
    


