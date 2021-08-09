import itertools
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
from sklearn.neighbors import NearestNeighbors


metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

def calc_pairwise_neighbor(n_x, n_y, rmsd = 0.25):
    '''
    Use sklean NearestNeighbors
    '''

    neigh_y = NearestNeighbors(radius= rmsd) 
    neigh_y.fit(n_y)

    x_in_y = neigh_y.radius_neighbors(n_x)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])

    return x_in_y, x_has_y

def calc_points(pdbss):
    pointss = []
    for ps in pdbss:
        points = []
        for p in ps:
            points.append(p.select(metal_sel)[0].getCoords())
        pointss.append(np.array(points))
    return pointss

def calc_neighbor(pointss, metal_sel = 'ion or name NI MN ZN CO CU MG FE' ):
    '''
    '''
   
    neighbor_pair_dict = {}
    for i in range(len(pointss)):
        for j in range(i+1, len(pointss)):           
            x_in_y, x_has_y = calc_pairwise_neighbor(pointss[i], pointss[j])
            y_in_x, y_has_x = calc_pairwise_neighbor(pointss[j], pointss[i])
            if x_has_y and y_has_x:
                neighbor_pair_dict[(i, j)] = (x_in_y)
                neighbor_pair_dict[(j, i)] = (y_in_x)

    return neighbor_pair_dict

def calc_comb_neighbor(neighbor_pair_dict, wins):
    '''
    This method try to keep points in one win showed in all other wins.
    It only guarantee pair-wise exist.
    '''
    comb_dict = {}

    for x, y in itertools.permutations(wins, 2):
        if (x, y) not in neighbor_pair_dict.keys():
            return None
    
    keys = list(itertools.permutations(wins, 2))

    count = len(wins)
    for i in range(len(wins)):
        keygroup = keys[i* (count - 1): i* (count - 1) + (count - 1)]
        
        group_i_count = len(neighbor_pair_dict[keygroup[0]][1])
        exists = [True for j in range(group_i_count)]
        for g in range(group_i_count):

            for k in keygroup:
                if len(neighbor_pair_dict[k][1][g]) <= 0:
                    exists[g] = False              

        comb_dict[wins[i]] = exists

    return comb_dict

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
    return

def neighbor2pymol(pointss, comb_dict, outdir):

    for id in range(len(pointss)):
        filename = 'chain' + str(id)
        write2pymol(pointss[id], outdir, filename)

    for key in comb_dict.keys():
        exists = comb_dict[key]
        ps_vs = [pointss[key][v] for v in range(len(exists)) if exists[v]]
        filename = str(key) + '_psvs'
        write2pymol(ps_vs, outdir, filename)
    return 


### Start neighbor search

workdir = '/mnt/e/DesignData/ligands/CoiledCoil/C4_rosetta/Rosetta_Output/eva_out_CUNI_M1-1/mems/'

outdir = workdir + 'neighbor_visualization/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

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

#Get neighbor pair dict.

pointss = calc_points(pdbss)
neighbor_pair_dict = calc_neighbor(pointss)

wins = [0, 1, 2, 3]

#Calc overall neighbors.

comb_dict = calc_comb_neighbor(neighbor_pair_dict, wins)

neighbor2pymol(pointss, comb_dict, outdir)



