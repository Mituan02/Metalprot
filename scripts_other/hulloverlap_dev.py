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

pdbs = extract_vdm.get_all_pbd_prody(workdir)

vdms = []
for pdb in pdbs:
    vdms.append(Query(pdb))

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

points1 = []
points2 = []
for p in pdbs:
    if 'mem0' in p.getTitle():
        points1.append(p.select(metal_sel)[0].getCoords())
    elif 'mem1' in p.getTitle():
        points2.append(p.select(metal_sel)[0].getCoords())
points1 = np.array(points1)
points2 = np.array(points2)

        
points1hull = Delaunay(points1)
points2hull = Delaunay(points2)

x = points1hull.find_simplex(points2) > 0

y = points2hull.find_simplex(points1) > 0

points_1in2 = np.array(points1[y])
points_2in1 = np.array(points2[x])


### Visualization

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def plot_tri(ax, points, tri, points2, tri2):
    for tr in tri.simplices:
        pts = points[tr, :]
        ax.plot3D(pts[[0,1],0], pts[[0,1],1], pts[[0,1],2], color='g', lw='0.1')
        ax.plot3D(pts[[0,2],0], pts[[0,2],1], pts[[0,2],2], color='g', lw='0.1')
        ax.plot3D(pts[[0,3],0], pts[[0,3],1], pts[[0,3],2], color='g', lw='0.1')
        ax.plot3D(pts[[1,2],0], pts[[1,2],1], pts[[1,2],2], color='g', lw='0.1')
        ax.plot3D(pts[[1,3],0], pts[[1,3],1], pts[[1,3],2], color='g', lw='0.1')
        ax.plot3D(pts[[2,3],0], pts[[2,3],1], pts[[2,3],2], color='g', lw='0.1')

    ax.scatter(points[:,0], points[:,1], points[:,2], color='b')

    for tr in tri2.simplices:
        pts = points2[tr, :]
        ax.plot3D(pts[[0,1],0], pts[[0,1],1], pts[[0,1],2], color='y', lw='0.1')
        ax.plot3D(pts[[0,2],0], pts[[0,2],1], pts[[0,2],2], color='y', lw='0.1')
        ax.plot3D(pts[[0,3],0], pts[[0,3],1], pts[[0,3],2], color='y', lw='0.1')
        ax.plot3D(pts[[1,2],0], pts[[1,2],1], pts[[1,2],2], color='y', lw='0.1')
        ax.plot3D(pts[[1,3],0], pts[[1,3],1], pts[[1,3],2], color='y', lw='0.1')
        ax.plot3D(pts[[2,3],0], pts[[2,3],1], pts[[2,3],2], color='y', lw='0.1')

    ax.scatter(points2[:,0], points2[:,1], points2[:,2], color='r')


fig = plt.figure()
ax = plt.axes(projection='3d')
plot_tri(ax, points1, points1hull, points2, points2hull)
plt.savefig('test.png')

### Write points to PyMol

points = np.concatenate((points1, points2), axis = 0)
points_overlap = np.concatenate((y, x), axis = 0)
names = ['NI' if points_overlap[i] else 'CU'for i in range(len(points))]
resnums = [i for i in range(len(points))]

chains = []
for i in range(len(points)):
    if i < len(points1):
        chains.append('A')
    elif i >= len(points1):
        chains.append('B')

def hull2pymol(points, names, resnums, chains, outdir, filename):
    mm = pr.AtomGroup('MetalMol')
    mm.setCoords(points)
    mm.setResnums(resnums)
    mm.setNames(names)
    mm.setResnames(names)
    mm.setChids(chains)

    pr.writePDB(outdir + filename + '.pdb', mm)

outdir = workdir
filename = 'Test_hull'
hull2pymol(points, names, resnums, chains, outdir, filename)

