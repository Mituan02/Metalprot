import prody as pr
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import argmin
import itertools

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

### calc pairwise angle dist info.

def calc_pair_geometry(pdb):
    '''
    Calc paired query angle and distance of binding core. 
    '''
    #Get metal, binding atom for each binding atom
    cts = []

    metal = pdb.select(metal_sel)[0]
    _contact_aas = pdb.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal.getResindex()))
    #For each aa, only select one contact atom. 
    resindices = np.unique(_contact_aas.getResindices())
    for rid in resindices:
        #TO DO: Not the first one, but the closest one for the  such ASP contains two contact atoms.
        _ct = _contact_aas.select('resindex ' + str(rid))
        if len(_ct) > 1:
            dists = [None]*len(_ct)
            for i in range(len(_ct)):
                dist = pr.calcDistance(metal, _ct[i])
                dists[i] = dist
            ct = _ct[argmin(dists)]
        else:
            ct = _ct[0]
        cts.append(ct)

    ct_len = len(cts)

    pair_infos = []

    for i, j in itertools.combinations(range(ct_len), 2):  
        pair_info = [] 

        pair_info.append((cts[i].getResname(), cts[j].getResname()))
        pair_info.append(str(cts[i].getResindex()) + '-' + str(cts[j].getResindex()))
        dist = pr.calcDistance(cts[i], cts[j])
        pair_info.append(dist)
        angle = pr.calcAngle(cts[i], metal, cts[j])
        pair_info.append(angle)

        pair_infos.append(pair_info)
    return pair_infos

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

