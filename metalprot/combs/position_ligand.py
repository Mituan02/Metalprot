import prody as pr
import numpy as np
from scipy.spatial.transform import Rotation
from sklearn.neighbors import NearestNeighbors
import math
import os
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem.rdmolfiles import MolToPDBFile

from ..basic import prody_ext 


def rotate_ligs(orign_lig, rot, rest, rotation_degree, interMolClashSets, interclash_dist = 3.0):
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

        if ligand_rot_is_clash(_lig, interMolClashSets, interclash_dist):
            continue
        all_ligs.append(_lig)

    return all_ligs


def generate_rotated_ligs(lig, rots, rests, rotation_degrees, interMolClashSets, interclash_dist = 3.0):
    '''
    The method only works for 2 rots.
    TO DO: use recursive algorithm to allow more than 2 rots.
    '''
    all_ligs = []
    i = 0
    rot = rots[i]
    rest = rests[i]
    degree = rotation_degrees[i]
    ligs = rotate_ligs(lig, rot, rest, degree, interMolClashSets, interclash_dist)
    for _lig in ligs:
        i = 1
        rot = rots[i]
        rest = rests[i]
        degree = rotation_degrees[i]
        ligs2 = rotate_ligs(_lig, rot, rest, degree, interMolClashSets, interclash_dist)
        all_ligs.extend(ligs2)
    return all_ligs


# def generate_rotated_ligs_rdkit(lig_smiles, total, outdir):
#     m = Chem.MolFromSmiles(lig_smiles)
#     m2=Chem.AddHs(m)
#     AllChem.EmbedMolecule(m2)

#     cids = AllChem.EmbedMultipleConfs(m2, numConfs=total)
#     for i in range(total):
#         MolToPDBFile(m2, outdir + 'rdkit_' + str(i) + '.pdb', confId = i)

#     m2s = []
#     for p in os.listdir(outdir):
#         if '.pdb' not in p:
#             continue
#         m2s.append(pr.parsePDB(outdir + p))
#     return m2s


def add_metal2lig(lig, rig, lig_sel, rig_sel, metal):
    '''
    Add metal to ligs that do not contain metal for superimpose on binding geometry.
    The rdkit generated ligands generally do not contain metal.

    rig: the ligand-metal prody obejct extracted from pdb.
    '''
    prody_ext.ordered_sel_transformation(lig, rig, lig_sel, rig_sel)
    
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


def ligand_rot_is_clash_depre(lig, rot, rest, interclash_dist = 3.0):
    '''
    The ligand it self could clash after rotation.
    The idea is to sep the ligand by rot into 2 groups anc check their dists.

    return True if clash
    '''
    atoms_rot = lig.select('heavy and not name ' + ' '.join(rot) + ' ' + ' '.join(rest)).getCoords()
    rest_coords = lig.select('heavy and name ' +  ' '.join(rest)).getCoords()

    nbrs = NearestNeighbors(radius= interclash_dist).fit(atoms_rot)
    adj_matrix = nbrs.radius_neighbors_graph(rest_coords).astype(bool)

    if np.sum(adj_matrix) >0:
        return True
    return False


def ligand_rot_is_clash(lig, interMolClashSets, interclash_dist = 3.0):
    '''
    The ligand it self could clash after rotation.
    The idea is to check the manually defined interMolClashSets 
    Such as interMolClashSets = [(['O1'], ['C1', 'C5']), (['O2'], ['C1', 'C5', 'C6'])], which means the rotation may induce clash between O1 and C1/C5.
    return True if clash
    '''
    for clashSet in interMolClashSets:
        atoms_rot = lig.select('name ' + ' '.join(clashSet[0])).getCoords()
        rest_coords = lig.select('name ' +  ' '.join(clashSet[1])).getCoords()

        nbrs = NearestNeighbors(radius= interclash_dist).fit(atoms_rot)
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
        try:
            mobile_sel_coords.append(_lig.select('name ' + s).getCoords()[0])
        except: 
            print('ERROR: (lig_2_ideageo) ' + ' '.join(lig_connect_sel))
            print(_lig.getNames())
            print(s)

    ideal_geo_sel_coords = []
    for s in geo_sel.split(' '):
        try:
            ideal_geo_sel_coords.append(ideal_geo_o.select('name ' + s).getCoords()[0])
        except:
            print('ERROR: (lig_2_ideageo lig) ' + ' '.join(lig_connect_sel))
            print('ERROR: (lig_2_ideageo) ' + geo_sel)
            print(s)

    transformation = pr.calcTransformation(np.array(mobile_sel_coords), np.array(ideal_geo_sel_coords))

    for lg in ligs:
        transformation.apply(lg)

    return


def lig_2_target(ligs, lig_connect_sel, target, geo_sel = 'chid X and name FE1 O2 O3'):
    '''
    supperimpose the ligand to the ideal metal binding geometry.
    '''
    _lig = ligs[0]

    mobile_sel_coords = []
    for s in lig_connect_sel:
        try:
            #print(_lig.getNames())
            mobile_sel_coords.append(_lig.select('name ' + s).getCoords()[0])
        except: 
            print('ERROR: (lig_2_ideageo) ' + ' '.join(lig_connect_sel))
            print(_lig.getNames())
            print(s)

    #print(target.select('chid X').getNames())
    try:
        target_sel_coords = target.select(geo_sel).getCoords()
    except:
        print('ERROR: (lig_2_ideageo lig) ' + ' '.join(lig_connect_sel))
        print('ERROR: (lig_2_ideageo) ' + geo_sel)


    transformation = pr.calcTransformation(np.array(mobile_sel_coords), np.array(target_sel_coords))

    for lg in ligs:
        transformation.apply(lg)

    return


def _ligand_clashing_filteredid(target_coords, lig_coords, lig_len, dist):

    nbrs = NearestNeighbors(radius= dist).fit(target_coords)
    adj_matrix = nbrs.radius_neighbors_graph(lig_coords).astype(bool)

    adj_matrix_reshape = adj_matrix.reshape((-1, adj_matrix.shape[1]*lig_len))
    successed = []
    #>>> TO DO: resahpe the adj_matrix (type csr_matrix) will improve the code.
    for i in range(adj_matrix_reshape.shape[0]):
        if adj_matrix_reshape.getrow(i).toarray().any():
            continue
        successed.append(i)

    return successed

def ligand_clashing_filter(ligs, target, dist = 3):
    '''
    The ligand clashing: the ligs cannot have any heavy atom within 3 A of a target bb.
    Nearest neighbor is used to calc the distances. 
    '''
    all_coords = []
    ids = []

    for i in range(len(ligs)):
        coords = ligs[i].select('heavy').getCoords()
        all_coords.extend(coords)
        ids.extend([i for j in range(len(coords))])

    lig_len = len(ligs[0].select('heavy'))
    target_coords = target.select('name N C CA O CB').getCoords()
    
    successed_id = _ligand_clashing_filteredid(target_coords, all_coords, lig_len, dist)
    
    filtered_ligs = []
    for i in successed_id:
        filtered_ligs.append(ligs[i])

    return filtered_ligs


def write_ligands(outdir, filtered_ligs, all_ligs = None, write_all_ligands = False):

    os.makedirs(outdir, exist_ok= True)

    os.makedirs(outdir + 'filtered_ligs/', exist_ok=True)

    for lg in filtered_ligs:
        pr.writePDB(outdir + 'filtered_ligs/' + lg.getTitle() + '.pdb.gz', lg)

    if write_all_ligands:
    # output all aligned ligands.
        os.makedirs(outdir + 'all_ligs/', exist_ok=True)
        for lg in all_ligs:
            pr.writePDB(outdir + 'all_ligs/' +  lg.getTitle() + '.pdb.gz', lg)

    '''
    # For benchmarking, the nature ligand exist. Try to get the minimum superimposed artificial ligand.
    nature_lig = pr.parsePDB(lig_path)
    min_RMSD = 100
    min_lg = None
    for lg in all_ligs:
        rmsd = pr.calcRMSD(lg, nature_lig)
        if rmsd < min_RMSD:
            min_RMSD = rmsd
            min_lg = lg
    print(min_RMSD)
    print(min_lg.getTitle())
    pr.writePDB(workdir + '_min_' + min_lg.getTitle(), min_lg)
    '''

    return 


def run_ligand(outdir, target, lig_path, ro1, ro2, rest1, rest2, lig_connects, geo_sel, rot_degree, interMolClashSets, clash_dist = 2.7, write_all_ligands = False):
    '''
    Generate all potential ligands for each binding position. 
    '''
    lig = pr.parsePDB(lig_path)

    all_ligs = generate_rotated_ligs(lig, [ro1, ro2], [rest1, rest2], rot_degree, interMolClashSets, interclash_dist= 3.0)

    # points = np.array(position_ligand.fibonacci_sphere(10, scale=0.2))
    # point_sel = 'name FE1'

    filtered_ligs, all_ligands = generate_ligands(all_ligs, target, lig_connects, geo_sel, clash_dist = clash_dist)

    write_ligands(outdir, filtered_ligs, all_ligands, write_all_ligands)

    return


def extract_ligand(outdir, target, lig_name = None):
    if not lig_name:
        lig = target.select('not protein')
    lig = target.select('resname ' + lig_name)
    os.makedirs(outdir + 'filtered_ligs/', exist_ok=True)
    pr.writePDB(outdir + 'filtered_ligs/' + lig_name + '.pdb', lig)


def fibonacci_sphere(samples=20, scale = 0.1):
    '''
    Generating points distributed evenly on a sphere. For metal position stimulation. 
    Based on: https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    '''
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((scale*x, scale*y, scale*z))

    return points



def generate_ligands(all_ligs, target, lig_connects, geo_sel, points = None, point_sel = None, clash_dist = 2.7):
    '''
    Generate all potential artificial ligand positions.

    points are from the stimulated potential metal positions. 
    '''
    all_ligands = []
    for i in range(len(lig_connects)):
        lig_connect = lig_connects[i]
        _ligs = [l.copy() for l in all_ligs]
        [l.setTitle('Geo_' + str(i) + '_P_0_' + l.getTitle() ) for l in _ligs]
        lig_2_target(_ligs, lig_connect, target, geo_sel = geo_sel)
        all_ligands.extend(_ligs)

        if points == None:
            continue
        tr = pr.calcTransformation(np.zeros((1, 3)), target.select(point_sel).getCoords())
        _points = tr.apply(points)

        for j in range(1, 1+len(_points)):
            p = _points[j-1]   
            _target = target.copy()
            _target.select(point_sel).setCoords(p)
            _ligs = [l.copy() for l in all_ligs]
            [l.setTitle('Geo_' + str(i) + '_P_' + str(j) + '_' + l.getTitle() ) for l in _ligs]
            lig_2_target(_ligs, lig_connect, _target, geo_sel = geo_sel)
            all_ligands.extend(_ligs)
    

    if len(all_ligands) <= 0:
        print('No ligands generated.')
    filtered_ligs = ligand_clashing_filter(all_ligands, target, dist = clash_dist)

    if len(filtered_ligs) <= 0:
        print('The position could not support the total ligand {}.'.format(len(all_ligands)))

    return filtered_ligs, all_ligands


