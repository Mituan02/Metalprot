import prody as pr
import itertools
import os
import numpy as np
from scipy.spatial.distance import cdist

def get_all_pbd_prody(workdir):
    '''
    find all .pdb file in a folder and load them with prody.
    return a list of pdb_prody.
    '''
    pdbs = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            pdbs.append(pdb_prody)
        except:
            print('not sure')   
    return pdbs

def get_contact_map(target, win_filter = None):
    '''
    calculate contact map for 2aa_sep database.
    return the ordered distance array and resindex array.
    '''
    xyzs = []
    for c in target.select('protein and name CA').getCoords():
        xyzs.append(c)
    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    dist_array = []
    id_array = []

    for i in range(len(xyzs)):
        for j in range(i+1, len(xyzs)):
            if win_filter:
                if i not in win_filter or j not in win_filter:
                    continue
            dist_array.append(dists[i, j])  
            id_array.append((i, j))
    dist_array, id_array = zip(*sorted(zip(dist_array, id_array)))
    return dist_array, id_array, dists


def extract_aa_comb(pdb, id_array_dist, outdir):
    resinds = np.unique(pdb.select('resname HIS CYS GLU ASP').getResindices())
    print(len(resinds))
    all_cbs = []

    un_dist_dict = set()
    for cbs in itertools.combinations(resinds, 3):
        #print(cbs)
        dist_ok = True
        for x,y in itertools.combinations(cbs, 2):
            if (x,y) in un_dist_dict:
                dist_ok = False
                break
            if  id_array_dist[(x, y)]> 10:
                un_dist_dict.add((x, y))
                dist_ok = False
                break
        if dist_ok:
            all_cbs.append(cbs)

    for cbs in itertools.combinations(resinds, 4):
        print(cbs)
        dist_ok = True
        for x,y in itertools.combinations(cbs, 2):
            if (x,y) in un_dist_dict:
                dist_ok = False
                break
            if id_array_dist[(x, y)] > 10:
                dist_ok = False
                break
        if dist_ok:
            all_cbs.append(cbs)

    count = 0
    for cbs in all_cbs:
        pdb_sel = pdb.select('resindex ' + ' '.join([str(c) for c in cbs]))
        pr.writePDB(outdir + pdb.getTitle().split('.')[0] + str(count), pdb_sel)
        count += 1


def generate_decoy(pdb, outdir):
    '''
    #Test example
    workdir = '/mnt/e/DesignData/ligands/NoLigand/NoMetalProt_1-550/'

    outdir = workdir + 'output_1a2o/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    pdb = pr.parsePDB(workdir + '1a2o.pdb')

    generate_decoy(pdb, outdir)

    '''
    dist_array, id_array, dists = get_contact_map(pdb)

    id_array_dist = {}
    for i in range(len(id_array)):
        id = id_array[i]
        id_array_dist[id] = dist_array[i]

    extract_aa_comb(pdb, id_array_dist, outdir)

    return

def generate_decoy_all(workdir, outpath = 'metal_decoys', ):
    '''
    Generate decoy metal binding cores for all decoy pdbs. 
    '''

    pdbs = get_all_pbd_prody(workdir)

    outdir = workdir + outpath

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for pdb in pdbs:
        generate_decoy(pdb, outdir)

    return 

