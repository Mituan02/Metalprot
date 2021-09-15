import os
import sys
from numpy.lib.shape_base import split
import prody as pr
import numpy as np
#You can either add the python package path.
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.search import extract_vdm
from metalprot.basic import hull
import pickle


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210911/'

# Get query pdbs, Add the cluster metal points into the query.hull_points

centroid_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi', 'cluster'], score_cut = 0, clu_num_cut = 0)


centroid_query_dict = {}
for i in range(len(centroid_querys)):
    splits = centroid_querys[i].query.getTitle().split('_')
    key = splits[1] + '_' + '_'.join(splits[7:12])
    centroid_query_dict[key] = i

#Prepare all the data, this should be optimized in the future.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 
align_sel = 'name N CA C'
query_all_metal = centroid_querys[0].copy()

all_coords = []
id_cluster_dict = {}
cluster_centroid_dict = {}

for query_id in range(len(centroid_querys)):
    query = centroid_querys[query_id].copy()
    inds = np.unique(query.query.getResindices())
    ind = inds[int(inds.shape[0]/2)-1]
    transform = pr.calcTransformation(query.query.select(align_sel + ' and resindex ' + str(ind)), query_all_metal.query.select(align_sel + ' and resindex ' + str(query_all_metal.contact_resind)))
    transform.apply(query.query)

    all_coords.append(query.query.select(metal_sel)[0].getCoords())

    cluster_key = query.get_cluster_key()
    cluster_centroid_dict[cluster_key] = query
    cluster_coords = []
    selfcenter_cluster_queryid = []

    id_cluster_dict[query_id] = cluster_key

    clu = extract_vdm.get_vdm_cluster(query)

    for q in clu.querys:
        transform.apply(q.query)
        cluster_coords.append(q.query.select(metal_sel)[0].getCoords())

        splits = q.query.getTitle().split('_')
        if 'centroid' in splits:
            key = splits[1] + '_' + '_'.join(splits[7:12])
        else:
            key = splits[1] + '_' + '_'.join(splits[6:11])
        
        selfcenter_cluster_queryid.append(centroid_query_dict[key])
    cluster_centroid_dict[cluster_key].hull_ag = hull.transfer2pdb(cluster_coords)
    cluster_centroid_dict[cluster_key].selfcenter_cluster_queryid = selfcenter_cluster_queryid

query_all_metal.hull_ag = hull.transfer2pdb(all_coords)

'''
#plot the coords of all metals.
points = query_all_metal.hull_ag.getCoords()
hull.write2pymol(points, query_dir, 'align_' + query_all_metal.query.getTitle())
'''

outdir = query_dir + 'pickle/'
os.makedirs(outdir, exist_ok= True)

with open(outdir + 'AllMetalQuery.pkl', 'wb') as f:
    pickle.dump(query_all_metal, f)

with open(outdir + 'AAMetalPhiPsi.pkl', 'wb') as f:
    pickle.dump(centroid_querys, f)

with open(outdir + 'cluster_centroid_dict.pkl', 'wb') as f:
    pickle.dump(cluster_centroid_dict, f)

with open(outdir + 'id_cluster_dict.pkl', 'wb') as f:
    pickle.dump(id_cluster_dict, f)


