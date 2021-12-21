import prody as pr
import numpy as np
from scipy.spatial.transform import Rotation
from sklearn.neighbors import NearestNeighbors


def rotate_porphyrin_first(lig, rot, rotation_degree = 1, metal = 'FE'):
    '''
    Note that the FE must be the last element of the ligand.
    '''
    # transform the lig to z axis.
    rot_coords = lig.select('name ' + ' '.join(rot)).getCoords()
    rot_dist = pr.calcDistance(rot_coords[1], rot_coords[0])

    z_coords = np.zeros((2, 3), dtype=float)
    z_coords[0, 0] = rot_dist

    pr.calcTransformation(rot_coords, z_coords).apply(lig)

    all_ligs = []
    for i in range(0, 12, rotation_degree):
        _lig = lig.copy()
        rotation = Rotation.from_rotvec(np.radians(i)*np.array([1, 0, 0]))        
        _coords = rotation.apply(_lig.getCoords())
        _coords[0] = _lig.select('name ' + 'NE2')[0].getCoords()
        _lig.setCoords(_coords)

        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)
        all_ligs.append(_lig)

    return all_ligs


def rotate_porphyrin_second(lig, rot, rotation_degree = 5, metal = 'FE'):
    '''
    Note that the FE must be the last element of the ligand.
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


def generate_rotated_porphyrins(lig, ro1, ro2, rotation_degree1 = 1, rotation_degree2 = 5, metal = 'ZN'):

    all_ligs = []
    
    rot1_ligs = rotate_porphyrin_first(lig, ro1, rotation_degree=rotation_degree1, metal= metal)
    for lig1 in rot1_ligs:
        _ligs = rotate_porphyrin_second(lig1, ro2, rotation_degree=rotation_degree2, metal= metal)
        all_ligs.extend(_ligs)

    return all_ligs


def porphyrin2vdm(ligs, lig_connects, vdm, metal_sel = 'name NI MN ZN CO CU MG FE'):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''
    _lig = ligs[0]

    lig_sel = _lig.select('name ' + ' '.join(lig_connects)).getCoords()

    vdm_coords = []
    vdm_coords.append(vdm.get_contact_coord())
    vdm_coords.extend(vdm.query.select(metal_sel).getCoords())


    transformation = pr.calcTransformation(lig_sel, np.array(vdm_coords))

    for lg in ligs:
        transformation.apply(lg)

    return 