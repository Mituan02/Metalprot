import prody as pr
import numpy as np
from scipy.spatial.transform import Rotation



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


def generate_rotated_ligs(lig, ro1, ro2, rotation_degree = 5, metal = 'ZN'):

    all_ligs = []
    
    rot1_ligs = rotate_ligs_first(lig, ro1, rotation_degree=rotation_degree, metal= metal)
    for lig1 in rot1_ligs:
        _ligs = rotate_ligs_second(lig1, ro2, ro1, rotation_degree=rotation_degree, metal= metal)
        all_ligs.extend(_ligs)

    return all_ligs


def calc_lig2ideageo_transformation(_lig, lig_connects, ideal_geo_o = None, geo_sel = 'name OE2 ZN'):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''
    lig_sel = _lig.select('name ' + ' '.join(lig_connects))

    ideal_geo_sel = ideal_geo_o.select(geo_sel)

    transformation = pr.calcTransformation(lig_sel, ideal_geo_sel)

    return transformation




