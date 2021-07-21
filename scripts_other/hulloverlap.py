import os
import sys
from metalprot.apps.quco import Comb
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot import ligand_database as ldb
from metalprot import extract_vdm
from metalprot.apps.quco import Query, Comb
import numpy as np
from prody.atomic import pointer
from scipy.spatial import Delaunay

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_rosetta/Rosetta_Output/eva_out_CUNI_M1-1/mems/'

outdir = workdir + 'hull_visualization/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

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
    points_overlap = np.concatenate((points_1in2, points_2in1), axis=0)
    return x, y, points_overlap

def calc_hull(pdbss, metal_sel = 'ion or name NI MN ZN CO CU MG FE' ):
    '''
    '''
    pointss = []
    for ps in pdbss:
        points = []
        for p in ps:
            points.append(p.select(metal_sel)[0].getCoords())
        pointss.append(np.array(points))
        
    hull_dict = {}
    points_rec = pointss[0]
    for i in range(len(pointss)):
        for j in range(i+1, len(pointss)):           
            x, y, points_overlap = calc_pairwise_hull(pointss[i], pointss[j])
            hull_dict[(i, j)] = (pointss[i], x, pointss[j], y, points_overlap)

            a,b,points_rec = calc_pairwise_hull(points_rec, pointss[j])

        
    return pointss, hull_dict, points_rec

def write2pymol(points, outdir, filename):
    names = ['NI' for i in range(len(points))]
    resnums = [i for i in range(len(points))]
    chains = ['A' for i in range(len(points))]
    mm = pr.AtomGroup('MetalMol')
    mm.setCoords(points)
    mm.setResnums(resnums)
    mm.setNames(names)
    mm.setResnames(names)
    mm.setChids(chains)

    pr.writePDB(outdir + filename + '.pdb', mm)
    
def hull2pymol(pointss, hull_dict, points_rec, outdir):
    '''
    visualization through creating pymol object.
    '''
    for id in range(len(pointss)):
        filename = 'chain' + str(id)
        write2pymol(pointss[id], outdir, filename)

    for key in hull_dict.keys():
        ps1, x, ps2, y, ps_vs = hull_dict[key]

        filename = str(key[0]) + '_' + str(key[1]) + '_psvs'
        write2pymol(ps_vs, outdir, filename)

    filename = 'all_psvs'
    write2pymol(points_rec, outdir, filename)


pdbs = extract_vdm.get_all_pbd_prody(workdir)

pdbss = []
pdbss.append([])
pdbss.append([])
pdbss.append([])
pdbss.append([])
for p in pdbs:
    if 'mem0' in p.getTitle():
        pdbss[0].append(p)
    elif 'mem1' in p.getTitle():
        pdbss[1].append(p)
    elif 'mem2' in p.getTitle():
        pdbss[2].append(p)
    elif 'mem3' in p.getTitle():
        pdbss[3].append(p)  

pointss, hull_dict, points_rec = calc_hull(pdbss)

hull2pymol(pointss, hull_dict, points_rec, outdir)

