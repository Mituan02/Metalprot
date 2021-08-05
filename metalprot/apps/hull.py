from scipy.spatial import Delaunay
import prody as pr
import numpy as np

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

def calc_pairwise_hull(points1, points2):
    '''
    Use Delaunay data structure to calc hull.
    points1, points2 are two sets of atom coords from different clusters.
    return x: points1 in points2hull
    return y: points2 in points1hull

    return points_overlap: all overlap points
    '''
    points1hull = Delaunay(points1)
    points2hull = Delaunay(points2)
    x = points2hull.find_simplex(points1) > 0
    y = points1hull.find_simplex(points2) > 0

    points_1in2 = points1[x]
    points_2in1 = points2[y]
    ol_points = np.concatenate((points_1in2, points_2in1), axis=0)
    overlap = False
    if x.any() and y.any():
        overlap = True
    return overlap, x, y, ol_points

def transfer2pdb(points, names = None, title = 'MetalMol'):
    if not names:
        names = ['NI' for i in range(len(points))]
    resnums = [i for i in range(len(points))]
    chains = ['A' for i in range(len(points))]
    mm = pr.AtomGroup(title)
    mm.setCoords(points)
    mm.setResnums(resnums)
    mm.setNames(names)
    mm.setResnames(names)
    mm.setChids(chains)
    return mm

def write2pymol(points, outdir, filename, names = None):
    mm = transfer2pdb(points, names)
    pr.writePDB(outdir + filename + '.pdb', mm)