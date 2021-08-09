import os
import sys
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot import extract_vdm
from metalprot import hull
from metalprot import search
import pickle


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/'

# Get query pdbs, Add the cluster metal points into the query.hull_points

centroid_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi'], score_cut = 0, clu_num_cut = 0)

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

with open(query_dir + 'pickle/AllMetalQuery.pkl', 'wb') as f:
    pickle.dump(query_all_metal, f)

with open(query_dir + 'pickle/AAMetalPhiPsi.pkl', 'wb') as f:
    pickle.dump(all_querys, f)

with open(query_dir + 'pickle/cluster_centroid_dict.pkl', 'wb') as f:
    pickle.dump(cluster_centroid_dict, f)

with open(query_dir + 'pickle/id_cluster_dict.pkl', 'wb') as f:
    pickle.dump(id_cluster_dict, f)