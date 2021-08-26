import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
from . import cluster
from . import transformation
from . import core

@dataclass
class clu_info:
    Metal: str
    clu_type: str
    clu_rmsd: float
    total_num: int
    clu_rank: int
    clu_num: int
    score: float

    def to_tab_string(self):
        clu_info = self.Metal + '\t' + self.clu_type + '\t' + str(self.clu_rmsd) + '\t' + str(self.total_num) + '\t' + str(self.clu_rank) + '\t' + str(self.clu_num) + '\t'+ str(self.score) 
        return clu_info

# Manipulate rcsb file

def organize_rcsb_file(workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"):
    '''
    The .csv files downloaded from rcsb database will be combined first, 
    then generate tab deliminated txt file.
    '''
    all_lines = []
    for file in os.listdir(workdir):
        if file.endswith(".csv"):
            with open(workdir + file, 'r') as f:
                all_lines.extend(f.readlines())
    with open(workdir + 'all_rcsb.txt', 'w') as f:
        f.write('\t'.join(all_lines[0].split(',')))
        for r in all_lines:
            if 'Entry ID' not in r and r.split(',')[0]!= '':
                f.write('\t'.join(r.split(',')))

# download rcsb pdb files

def download_pdb(workdir, filename, resolution = 2.5):
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    all_pdbs = []
    with open(workdir + filename, 'r') as f: 
        for line in f.readlines():
            #print(line)
            r = line.split('\t')
            #print(r)
            if r[0] == '"Entry ID"': continue
            if r[0] == '' or r[4]== '' or (',' in r[4]) or float(r[4].split('"')[1]) > 2.5:
                continue
            all_pdbs.append(r[0].split('"')[1])

    exist_pdb = set()
    for file in os.listdir(workdir):
        if file.endswith(".pdb.gz"):
            exist_pdb.add(file.split('.')[0].upper())

    pr.pathPDBFolder(workdir)
    for p in all_pdbs:
        if p in exist_pdb: continue
        pr.fetchPDBviaFTP(p, compressed = False) 

    # #Then unzip them in linux with:  
    # cd /mnt/e/DesignData/ligands/NI_rcsb/
    # gunzip *.gz

# Basic function. 

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

def writepdb(cores, outdir):
    '''
    cores = list of (pdb name, pdb_prody)
    '''

    if not os.path.exists(outdir):
            os.mkdir(outdir)
    for c in cores:
        outfile = c[0].split('.')[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

def load_cores(workdir):
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            _core = core.Core(pdb_prody)
            cores.append(_core)
        except:
            print('not sure')   
    return cores

# Write Summary

def write_clu_info(filename, clu_infos):
    '''
    Write information of cluters.
    @ clu_infos: [clu_info]
    '''
    with open(filename, 'w') as f:
        f.write('Metal\tclu_type\tclu_rmsd\ttotal_num\tclust_rank\tclust_num\tscore\n')
        for r in clu_infos:
            f.write(r.to_tab_string() + '\n')  

def plot_clu_info(clu_infos, outplot):
    fig, (ax, ax1) =plt.subplots(2, 1, figsize=(12,8))
    
    x = list(range(1, len(clu_infos) + 1))
    y = [c.score for c in clu_infos]
    ax.plot(x, y)
    ax.hlines(y=0, xmin = 0, xmax = x[-1], color='r')
    ax.legend()
    #ax.set_xticks(x)
    ax.set_xlabel("Rank", fontsize = 12)
    ax.set_ylabel("vdM score", fontsize = 12)

    counts = [c.clu_num for c in clu_infos]
    ax1.bar(x, counts)

    plt.tight_layout()
    plt.savefig(outplot)
    plt.close()

# Prepare rcsb database. // extract seq within +-3 aa for each contact aa. 

def connectivity_filter(pdb_prody, ind, ext_ind):
    res1 = pdb_prody.select('protein and resindex ' + str(ind))
    res2 = pdb_prody.select('protein and resindex ' + str(ext_ind))
    if not res2:
        return False
    if res1[0].getResnum() - res2[0].getResnum() == ind - ext_ind and res1[0].getChid() == res2[0].getChid() and res1[0].getSegname() == res2[0].getSegname():
        return True
    return False

def extend_res_indices(inds_near_res, pdb_prody, extend = 4):
    extend_inds = []
    inds = set()
    for ind in inds_near_res:
        for i in range(-extend, extend + 1):
            the_ind = ind + i
            if the_ind not in inds and the_ind>= 0 and connectivity_filter(pdb_prody, ind, the_ind):
                extend_inds.append(the_ind)
                inds.add(the_ind)
    return extend_inds

def get_metal_core_seq(pdb_prody, metal_sel, extend = 4):
    metal = metal_sel.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb_prody.select('resname HIS GLU ASP CYS and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur') or len(all_near.select('nitrogen or oxygen or sulfur')) < 3:
            continue          
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        # all_near_res = pdb_prody.select('protein and resindex ' + ' '.join([str(ind) for ind in inds]))
        # if not all_near_res or len(np.unique(all_near_res.getResindices())) < 2:
        #     continue     
        # inds_near_res =  all_near_res.getResindices()
        # ext_inds = extend_res_indices(inds_near_res, pdb_prody, extend)
        ext_inds = extend_res_indices(inds, pdb_prody, extend)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + metal + '_'+ str(count), sel_pdb_prody))        
    return metal_cores

def get_2ndshell_indices(inds, pdb_prody, ni_index):
    _2nd_resindices = []
    for ind in inds:
        if pdb_prody.select('resindex ' + str(ind)).getResnames()[0] == 'HIS':
            dist1 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name ND1')[0])
            dist2 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name NE2')[0])
            if dist1 < dist2:
                index = pdb_prody.select('resindex ' + str(ind) + ' name NE2')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name NE2')[0].getResindex()
            else:
                index = pdb_prody.select('resindex ' + str(ind) + ' name ND1')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name ND1')[0].getResindex()                         
        elif pdb_prody.select('resindex ' + str(ind)).getResnames()[0] == 'ASP':
            dist1 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name OD1')[0])
            dist2 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name OD2')[0])
            if dist1 < dist2:
                index = pdb_prody.select('resindex ' + str(ind) + ' name OD2')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name OD2')[0].getResindex()
            else:
                index = pdb_prody.select('resindex ' + str(ind) + ' name OD1')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name OD1')[0].getResindex()               
        elif pdb_prody.select('resindex ' + str(ind)).getResnames()[0] == 'GLU':
            dist1 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name OE1')[0])
            dist2 = pr.calcDistance(pdb_prody.select('index ' + str(ni_index))[0], pdb_prody.select('resindex ' + str(ind) + ' name OE2')[0])
            if dist1 < dist2:
                index = pdb_prody.select('resindex ' + str(ind) + ' name OE2')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name OE2')[0].getResindex()
            else:
                index = pdb_prody.select('resindex ' + str(ind) + ' name OE1')[0].getIndex()
                resindex = pdb_prody.select('resindex ' + str(ind) + ' name OE1')[0].getResindex()
        else:
            continue     
        all_near = pdb_prody.select('protein and heavy and within 3.4 of index ' + str(index) + ' and not resindex ' + str(resindex))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur'):
            continue
        inds_2nshell = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        _2nd_resindices.extend(inds_2nshell) 

    return _2nd_resindices

def get_metal_core_seq_2ndshell(pdb_prody, metal_sel, extend = 4):
    metal = metal_sel.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb_prody.select('protein and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur') or len(all_near.select('nitrogen or oxygen or sulfur')) < 3:
            continue          
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        ext_inds = extend_res_indices(inds, pdb_prody, extend)
        _2ndshell_inds = get_2ndshell_indices(inds, pdb_prody, ni_index)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ ' '.join([str(ind) for ind in _2ndshell_inds]) + ' ' + str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + metal + '_'+ str(count), sel_pdb_prody))        
    return metal_cores

def extract_all_core_seq_from_path(workdir, metal_sel, extend = 4, extract_2ndshell = False):
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            if extract_2ndshell:
                core = get_metal_core_seq_2ndshell(pdb_prody, metal_sel, extend)
            else:
                core = get_metal_core_seq(pdb_prody, metal_sel, extend)
            if core:
                cores.extend(core)
            pdb_prody = None
        except:
            print('not sure')   
    return cores

def extract_all_core_seq(pdbs, metal_sel, extend = 4):
    cores = []
    for pdb in pdbs:
        core = get_metal_core_seq(pdb, metal_sel, extend)
        if core:
            cores.extend(core)
    return cores

def superimpose_core_and_writepdb(cores, first, metal_sel, outdir):
    '''
    Just supperimpose on the metal. And write pdb.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for c in cores:
        a_coords = first[1].select(metal_sel)[0].getCoords().reshape(1,3)
        b_coords = c[1].select(metal_sel)[0].getCoords().reshape(1,3)

        pr.calcTransformation(b_coords, a_coords).apply(c[1])

        outfile = c[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

# Prepare rcsb database. // Reduce duplication.

def reduce_dup(pdbs, metal_sel):
    '''
    Cluster pdbs based on the structure similarity.
    '''
    clusters = []
    clustered = set()
    for i in range(len(pdbs)):
        if i in clustered: continue

        cluster = [i]
        clustered.add(i)
        for j in range(i + 1, len(pdbs)):
            if j in clustered: continue
            if len(pdbs[i].select('name N CA C O')) != len(pdbs[j].select('name N CA C O')):
                continue
            #print(str(i) + '---' + str(j))
            sequence_equal = False
            if pdbs[i].select('name CA').getSequence() == pdbs[i].select('name CA').getSequence():
                sequence_equal = True
            #TO DO: The calcTransformation here will change the position of pdb. 
            #This will make the output pdb not align well. Current solved by re align.
            pr.calcTransformation(pdbs[j].select('name N CA C O'), pdbs[i].select('name N CA C O')).apply(pdbs[j])
            rmsd = pr.calcRMSD(pdbs[i].select('name N CA C O'), pdbs[j].select('name N CA C O'))

            if (sequence_equal and rmsd < 0.8) or rmsd < 0.25:
                cluster.append(j)
                clustered.add(j)

        clusters.append(cluster)
    return clusters

def write_dup_summary(workdir, pdbs, clusters):
    with open(workdir + 'dup_summary.txt', 'w') as f:
        for clu in clusters:
            line = [pdbs[ind].getTitle() for ind in clu]
            f.write('\t'.join(line) + '\n')

def extract_rep_and_writepdb(pdbs, clusters, metal_sel, outdir):
    '''
    select represent pdb of each cluster and copy into a subfolder.
    TO DO: It is more efficient to copy the file into the subfolder.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sel_pdbs = []

    for clu in clusters:
        sel_pdbs.append((pdbs[clu[0]].getTitle(), pdbs[clu[0]]))

    superimpose_core_and_writepdb(sel_pdbs, sel_pdbs[0], metal_sel, outdir)


### clustring

def superimpose_aa_core(pdbs, rmsd = 0.5, len_sel = 5, align_sel = 'name C CA N O NI', min_cluster_size = 0):
    '''
    There are so many ways to superimpose aa and metal.
    This method try to algin the C-CA-N_NI
    '''
    clu = cluster.Cluster()
    clu.rmsd_cutoff = rmsd
    clu.pdbs = []

    for pdb in pdbs:
        c = pdb.select(align_sel).getCoords()       
        if len(c)!= len_sel: 
            continue
        clu.pdb_coords.append(c)
        clu.pdbs.append(pdb)
    if len(clu.pdb_coords) <= 0:
        return 
    clu.pdb_coords = np.array(clu.pdb_coords, dtype = 'float32')

    clu.make_pairwise_rmsd_mat()  
    if not clu._square:
        clu.make_square()
    if not clu._adj_mat:
        clu.make_adj_mat()
    clu.cluster(min_cluster_size = min_cluster_size)

    return clu

def get_clu_info_write(outfile, pdbs, clu, rmsd = 0.5, metal_sel = 'NI', align_sel = 'name C CA N O NI'):
    clu_infos = []
    n_avg = sum([len(m) for m in clu.mems])/len(clu.mems)
    for i in range(len(clu.mems)):
        c = clu_info(Metal=metal_sel, clu_type = align_sel, 
            clu_rmsd=rmsd, total_num = len(pdbs), clu_rank = i, 
            clu_num=len(clu.mems[i]), score = np.log(len(clu.mems[i])/n_avg) )
        clu_infos.append(c)
    write_clu_info(outfile, clu_infos)
    return clu_infos

def print_cluster_pdbs(clu, outdir, rmsd = 0.5, tag = ''):

    for i in range(len(clu.mems)):
        cluster_out_dir = outdir + str(i) + '/'
        # if not os.path.exists(cluster_out_dir):
        #     os.mkdir(cluster_out_dir)
        _print_cluster_rank_pdbs(clu, i, cluster_out_dir, tag)

def _print_cluster_rank_pdbs(clu, rank, outdir='./', tag=''):
    try: os.makedirs(outdir)
    except: pass
    try:
        cent = clu.cents[rank]
        mems = clu.mems[rank]

        # Align backbone of cluster centroid to backbone of centroid of largest cluster.
        R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[cent],
                                        clu.pdb_coords[clu.cents[0]])
        cent_coords = np.dot((clu.pdb_coords[cent] - m_com), R) + t_com

        for i, mem in enumerate(mems):
            R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[mem], cent_coords)
            pdb = clu.pdbs[mem].copy()
            pdb_coords = pdb.getCoords()
            coords_transformed = np.dot((pdb_coords - m_com), R) + t_com
            pdb.setCoords(coords_transformed)
            is_cent = '_centroid' if mem == cent else ''
            pr.writePDB(outdir + tag + 'cluster_' + str(rank) + '_mem_' + str(i)
                         + is_cent + '_' + pdb.getTitle().split('.')[0] + '.pdb', pdb)

    except IndexError:
        print('Cluster', rank, 'does not exist.')

def check_metal_diff(clu, metals = ['CA', 'CO', 'CU', 'FE', 'MG', 'MN', 'NI', 'ZN']):
    metal_diffs = []
    for rank in range(len(clu.mems)):
        if not clu.mems[rank].any(): continue

        mems = clu.mems[rank]  
        metal_diff = [0]*len(metals)
        for i, mem in enumerate(mems):
            pdb = clu.pdbs[mem]

            if pdb.select('ion and name CA'):
                metal_diff[0]+=1
            elif pdb.select('name CO'):
                metal_diff[1]+=1
            elif pdb.select('name CU'):
                metal_diff[2]+=1
            elif pdb.select('name FE'):
                metal_diff[3]+=1
            elif pdb.select('name MG'):
                metal_diff[4]+=1                
            elif pdb.select('name MN'):
                metal_diff[5]+=1
            elif pdb.select('name NI'):
                metal_diff[6]+=1
            elif pdb.select('name ZN'):
                metal_diff[7]+=1

        metal_diffs.append(metal_diff)
    return metal_diffs    

def write_metal_diff(workdir, metal_diffs, metals = ['CA', 'CO', 'CU', 'FE', 'MG', 'MN', 'NI', 'ZN']):
    with open(workdir + 'metal_diff.txt', 'w') as f:
        f.write('\t'.join(metals) + '\n')
        for md in metal_diffs:
            f.write('\t'.join([str(x) for x in md]) + '\n')


# run cluster

def run_cluster(_pdbs, workdir, outdir, rmsd, metal_sel, len_sel, align_sel, min_cluster_size = 0, tag = ''):
    
    clu = superimpose_aa_core(_pdbs, rmsd = rmsd, len_sel = len_sel, align_sel = align_sel, min_cluster_size = min_cluster_size)
    
    if not clu or len(clu.mems) == 0: return
    
    print_cluster_pdbs(clu, workdir + outdir, rmsd, tag)

    metal_diffs = check_metal_diff(clu)

    write_metal_diff(workdir + outdir, metal_diffs)

    clu_infos = get_clu_info_write(workdir + outdir + '_summary.txt', _pdbs, clu, rmsd = rmsd, metal_sel = metal_sel, align_sel = align_sel)

    plot_clu_info(clu_infos, workdir + outdir + '_score.png')

    return
