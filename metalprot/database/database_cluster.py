import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..basic import cluster
from ..basic import transformation

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

def run_cluster(_pdbs, workdir, outdir, rmsd, metal_sel, len_sel, align_sel, min_cluster_size = 0, tag = '', is_self_center = False):
    if not is_self_center:
        clu = superimpose_aa_core(_pdbs, rmsd = rmsd, len_sel = len_sel, align_sel = align_sel, min_cluster_size = min_cluster_size)
    else:
        clu = cluster_selfcenter(_pdbs, rmsd = rmsd, len_sel = len_sel, align_sel = align_sel, min_cluster_size = min_cluster_size)
    
    if not clu or len(clu.mems) == 0: 
        print('clu is None')
        return
    
    print_cluster_pdbs(clu, workdir + outdir, rmsd, tag)

    metal_diffs = check_metal_diff(clu)

    write_metal_diff(workdir + outdir, metal_diffs)

    clu_infos = get_clu_info_write(workdir + outdir + '_summary.txt', _pdbs, clu, rmsd = rmsd, metal_sel = metal_sel, align_sel = align_sel)

    plot_clu_info(clu_infos, workdir + outdir + '_score.png')

    return


### new cluster method.

def cluster_selfcenter(pdbs, rmsd = 0.5, len_sel = 5, align_sel = 'name C CA N O NI', min_cluster_size = 0):
    '''
    selfcenter cluster method.
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

    #Here is a change of the cluster greedy funciton.
    #Basicly, we try to extract cent for every pdb.
    all_mems = []
 
    indices = np.arange(clu.adj_mat.shape[0])
    try:
        for cent in range(clu.adj_mat.shape[0]):
            row = clu.adj_mat.getrow(cent)
            tf = row.toarray().astype(bool)[0]
            mems = indices[tf]
            all_mems.append(mems)

    except KeyboardInterrupt:
        pass
    cent_mems = list(zip(*sorted(enumerate(all_mems), key=lambda x: len(x[1]), reverse=True)))
    clu.mems = list(cent_mems[1])
    clu.cents = list(cent_mems[0])
    return clu