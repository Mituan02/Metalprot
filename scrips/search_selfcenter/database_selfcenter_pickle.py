'''
Here we load all vdm database into pkl file. pkl file could be load fast.

The database should be selfcentered.

all_vdms = [vdm]
all_metal_vdm = vdm with all member metal position.
centroid_query_dict = [clu_key:vdm id] for example '(HIS, 1234): 4567' 

'''

import os
import prody as pr
from metalprot.search import extract_vdm
from metalprot.basic import prody_ext
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/search_selfcenter/database_selfcenter_pickle.py
'''

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211013_selfcenter/'
query_dir = '/mnt/e/DesignData/ligands/all/20220116_FE_MN_CO/20220116_selfcenter/'

centroid_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi', 'cluster'], file_name_not_includes=['@', 'CYS'], score_cut = 0, clu_num_cut = 0)


centroid_query_dict = {}
for i in range(len(centroid_querys)):
    splits = centroid_querys[i].query.getTitle().split('_')
    key = splits[1] + '_' + '_'.join(splits[7:12])
    centroid_query_dict[key] = i

#Prepare all the data, this should be optimized in the future.

metal_sel = 'name NI MN ZN CO CU MG FE' 
align_sel = 'name N CA C'
all_metal_vdm = centroid_querys[0].copy()

all_vdms = []
cluster_centroid_dict = {}

all_coords = []

max_clu_num = 0
for query_id in range(len(centroid_querys)):
    print(query_id)
    query = centroid_querys[query_id].copy()

    #transformation
    ind = query.contact_resind
    transform = pr.calcTransformation(query.query.select(align_sel + ' and resindex ' + str(ind)), all_metal_vdm.query.select(align_sel + ' and resindex ' + str(all_metal_vdm.contact_resind)))
    transform.apply(query.query)

    #obtain all member metal coords.
    all_coords.append(query.get_metal_coord())

    cluster_coords = []
    selfcenter_cluster_queryid = []

    mem_vdm_names = extract_vdm.get_mem_vdm_names(query)
    if 'cluster_0' in query.query.getTitle():
        max_clu_num = len(mem_vdm_names)
    for vn in mem_vdm_names:
        splits = vn.split('_')
        if 'centroid' in splits:
            key = splits[1] + '_' + '_'.join(splits[7:12])
        else:
            key = splits[1] + '_' + '_'.join(splits[6:11])

        id = centroid_query_dict[key]       
        q = centroid_querys[id].copy()
        ## align to heavy of query
        #transform.apply(q.query)
        ## align to 'N CA C' of query
        #pr.calcTransformation(q.query.select(align_sel + ' and resindex ' + str(ind)), all_metal_vdm.query.select(align_sel + ' and resindex ' + str(all_metal_vdm.contact_resind))).apply(q.query)
        pr.calcTransformation(q.query.select('heavy'), query.query.select('heavy')).apply(q.query)
        
        cluster_coords.append(q.get_metal_coord())
        selfcenter_cluster_queryid.append(id)

    query.id = query_id
    query.clu_rank = int(query.query.getTitle().split('_')[3]) 
    query.max_clu_num = max_clu_num
    query.metal_atomgroup = prody_ext.transfer2pdb(cluster_coords)
    query.clu_member_ids = selfcenter_cluster_queryid
    all_vdms.append(query)

    cluster_key = query.get_cluster_key()
    cluster_centroid_dict[cluster_key] = query.id


all_metal_vdm.metal_atomgroup = prody_ext.transfer2pdb(all_coords)

'''
#plot the coords of all metals.
points = all_metal_vdm.get_metal_mem_coords()
hull.write2pymol(points, query_dir, 'align_' + all_metal_vdm.query.getTitle())
'''

#query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210916_2017_2018_selfcenter_alignBB/'
outdir = query_dir + 'pickle_noCYS/'
os.makedirs(outdir, exist_ok= True)

with open(outdir + 'all_metal_vdm.pkl', 'wb') as f:
    pickle.dump(all_metal_vdm, f)

with open(outdir + 'AAMetalPhiPsi.pkl', 'wb') as f:
    pickle.dump(all_vdms, f)

with open(outdir + 'cluster_centroid_dict.pkl', 'wb') as f:
    pickle.dump(cluster_centroid_dict, f)



