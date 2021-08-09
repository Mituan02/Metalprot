import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import extract_vdm
from metalprot import hull
from metalprot import search

'''
python /mnt/e/GitHub_Design/Metalprot/metalprot/scrips/run_neighbor_search.py
'''

### Generate queryss
queryss = []

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/'

# Get query pdbs, Add the cluster metal points into the query.hull_points

centroid_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi_HIS'], score_cut = 1, clu_num_cut = 50)

#Prepare all the data, this should be optimized in the future.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel = 'name N CA C'
query_all_metal = centroid_querys[0].copy()
all_querys = []
all_coords = []
id_cluster_dict = {}
cluster_centroid_dict = {}

query_id = 0
for query in centroid_querys:
    inds = np.unique(query.query.getResindices())
    ind = inds[int(inds.shape[0]/2)-1]
    transform = pr.calcTransformation(query.query.select(align_sel + ' and resindex ' + str(ind)), query_all_metal.query.select(align_sel + ' and resindex ' + str(query_all_metal.contact_resind)))
    transform.apply(query.query)

    cluster_key = query.get_cluster_key()
    cluster_centroid_dict[cluster_key] = query
    cluster_coords = []

    clu = extract_vdm.get_vdm_cluster(query)
    clu.realign_by_CCAN(target = query, align_sel=align_sel)
    for q in clu.querys:
        all_querys.append(q)
        all_coords.append(q.query.select(metal_sel)[0].getCoords())
        cluster_coords.append(q.query.select(metal_sel)[0].getCoords())

        id_cluster_dict[query_id] = cluster_key
        query_id += 1
    cluster_centroid_dict[cluster_key].hull_ag = hull.transfer2pdb(cluster_coords)

query_all_metal.hull_ag = hull.transfer2pdb(all_coords)

print(len(all_querys))


#contact_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['M8_AtomContact4_clusters'], score_cut = 0, clu_num_cut = 2)
contact_querys = None
#_2nd_querys = extract_vdm.extract_all_centroid(query_dir + '20210608/', summary_name = '_summary.txt', file_name_includes = ['M7_AA2sMetal-HIS_clusters'], score_cut = 0, clu_num_cut = 0)
_2nd_querys = None

### run Search_struct

workdir = '/mnt/e/DesignData/ligands/LigandBB/MID1sc10/'

outdir = workdir + 'output_neighbor_test/'

target_path = workdir + '5od1_zn.pdb'

rmsd_cuts = 0.25

num_iter = 3

win_filter = None

'''
#For testing natural metal binding proteins.
_pdb = pr.parsePDB(target_path)

cores = ligand_database.get_metal_core_seq(_pdb, metal_sel = 'ion or name NI MN ZN CO CU MG FE' , extend = 4)

print(cores[0][1].select('name CA').getSequence())

#pr.writePDB(target_path + '_core.pdb', cores[0][1])

win_filter = [x for x in np.unique(cores[0][1].getResindices())][0:-1]

print(win_filter)
'''

#win_filter = [30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
win_filter = [34,  60,  64]

ss =  search.Search_vdM(target_path, outdir, all_querys, id_cluster_dict, cluster_centroid_dict, query_all_metal, num_iter, rmsd_cuts, win_filter)


ss.neighbor_generate_query_dict()


ss.neighbor_generate_pair_dict()


ss.neighbor_search_wins()

#Test Graph
'''
win_comb = [34, 60, 64]

graph = Graph(win_comb, len(ss.querys))

graph.calc_pair_connectivity(ss.neighbor_pair_dict)

graph.get_paths()

'''
ss.neighbor_extract_query()

ss.neighbor_calc_comb_score()

ss.neighbor_calc_geometry()

ss.neighbor_write()

ss.neighbor_write_summary()