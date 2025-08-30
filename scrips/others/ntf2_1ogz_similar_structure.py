'''
The script tries to find similar design regions of 1ogz family protein.
'''
from numpy.core.fromnumeric import argmin
import prody as pr
from prody.proteins.pdbfile import writePDB
from sklearn.neighbors import NearestNeighbors
import numpy as np
import os

def calc_pairwise_neighbor(n_x, n_y, rmsd = 5):
    '''
    Use sklean NearestNeighbors, find every x point has how many neighbors from y.
    '''

    neigh_y = NearestNeighbors(radius= rmsd) 
    neigh_y.fit(n_y)

    x_in_y = neigh_y.radius_neighbors(n_x)
    return x_in_y[0], x_in_y[1]


workdir = '/mnt/e/DesignData/ligands/LigandBB/ntf2/'

#Manually defined logz design region
win_filter = [9, 10, 13, 17, 25, 29, 35, 37, 53, 54, 57, 77, 79, 81, 83, 85, 94, 96, 98, 100, 108, 111, 113, 115]

_1ogz = pr.parsePDB(workdir + '1ogz.pdb')

_1ogz.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getResnames()

### -------------------------------------------------------
### Based on pymol, the region of 3vsy agreed well with logz.

_3vsy = pr.parsePDB(workdir + '3vsy_align.pdb')

_3vsy.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getResnames()

coord_logz = _1ogz.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getCoords()

coord_3vsy = _3vsy.select('name CA').getCoords()

x_in_y_dist, x_in_y_ind = calc_pairwise_neighbor(coord_logz, coord_3vsy)

res_index = []

for i in range(len(x_in_y_ind)):
    if len(x_in_y_dist[i]) <= 0:
        continue
    min_dist_ind = np.argmin(x_in_y_dist[i])
    res_index.append(x_in_y_ind[i][min_dist_ind])

print(res_index)

#[9, 10, 13, 17, 25, 29, 35, 37, 53, 54, 57, 79, 81, 83, 85, 87, 96, 98, 100, 102, 110, 113, 115, 117]

### -------------------------------------------------------
### Based on pymol, the region of 6p44 doesn't agreed well with logz.

_6p44 = pr.parsePDB(workdir + 'align_1dmm.pdb')

coord_logz = _1ogz.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getCoords()

coord_6p44 = _6p44.select('name CA').getCoords()

x_in_y_dist, x_in_y_ind = calc_pairwise_neighbor(coord_logz, coord_6p44)

res_index = []

for i in range(len(x_in_y_ind)):
    if len(x_in_y_dist[i]) <= 0:
        continue
    min_dist_ind = np.argmin(x_in_y_dist[i])
    res_index.append(x_in_y_ind[i][min_dist_ind])

print(res_index)

# [8, 9, 12, 16, 24, 28, 34, 36, 52, 53, 56, 75, 77, 79, 81, 82, 86, 88, 90, 92, 100, 103, 105, 107]
    
### -------------------------------------------------------
### _5z3r

_5z3r = pr.parsePDB(workdir + '5z3r_align.pdb')

#_6p44.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getResnames()

coord_logz = _1ogz.select('name CA and resindex ' + ' '.join([str(x) for x in win_filter])).getCoords()

coord_5z3r = _5z3r.select('name CA').getCoords()

x_in_y_dist, x_in_y_ind = calc_pairwise_neighbor(coord_logz, coord_5z3r)

res_index = []

for i in range(len(x_in_y_ind)):

    min_dist_ind = np.argmin(x_in_y_dist[i])

    res_index.append(x_in_y_ind[i][min_dist_ind])

print(res_index)

# [6, 7, 10, 14, 22, 26, 32, 34, 57, 58, 61, 83, 85, 87, 89, 92, 101, 103, 105, 107, 115, 118, 120, 122]


### -------------------------------------------------------
### Superimpose all family_1ogz protein structures to 1ogz 13, 17, 25, 54

workdir = '/mnt/e/DesignData/ligands/LigandBB/ntf2/family_3vsy/'

structures = []
for s in os.listdir(workdir):
    if '.pdb' not in s:
        continue
    structures.append(s)

pdbs = []
for s in structures:
    pdb = pr.parsePDB(workdir + s)
    pdbs.append(pdb)

_target = [p for p in pdbs if '1ogz'in p.getTitle()][0]
t_sel = _target.select('protein and resindex 13:54 and bb')
for p in pdbs:
    p_sel = p.select('protein and resindex 14:55 and bb')

    transform = pr.calcTransformation(p_sel, t_sel).apply(p)
    pr.writePDB(workdir + 'align_' + p.getTitle(), p)
    