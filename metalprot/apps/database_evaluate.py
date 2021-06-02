import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

def ev_atom_distribution(pdbs, centroid_pdb, atom_sel = metal_sel, align_sel = 'bb'):

    centroid_metal_coord = centroid_pdb.select(atom_sel)[0].getCoords()
    
    metal_coords = []

    dists = []

    for pdb in pdbs:
        if align_sel:
            pr.calcTransformation(pdb.select(align_sel), centroid_pdb.select(align_sel)).apply(pdb)
        
        metal_coords.append(pdb.select(atom_sel)[0].getCoords() - centroid_metal_coord) #transform to (0, 0, 0)

        dist = pr.calcDistance(pdb.select(atom_sel)[0], centroid_pdb.select(atom_sel)[0])
        dists.append(dist)

    return np.array(metal_coords), np.array(dists)

def plt_dist_his(dist, outplot_name = 'dist.png'):
    n, bins, patches = plt.hist(x=dist, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Distance with centroid')
    plt.ylabel('Frequency')
    plt.title('Dist Histogram')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outplot_name)
    plt.close()

def plt_3d_points(metal_coords, outplot_name = 'points.png'):

    fig = plt.figure(figsize=(8,8))

    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(metal_coords[:, 0], metal_coords[:, 1], metal_coords[:, 2])
    plt.savefig(outplot_name)
    plt.close()

