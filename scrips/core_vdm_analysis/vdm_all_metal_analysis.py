'''
The purpose is to calculate the distribution of vdms from different metals along the number of clusters.
'''

import os
import prody as pr
from metalprot.search import extract_vdm
from metalprot.basic import hull
import pickle

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

query_dir = '/mnt/e/DesignData/ligands/all/20211207/20211209_selfcenter_nometal/'


centroid_querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['AAMetalPhiPsi', 'cluster', 'ASP'], file_name_not_includes=['@'], score_cut = 0, clu_num_cut = 0)

metals = [q.query.select(metal_sel)[0].getName()  for q in centroid_querys]

centroid_query_dict = {}
for i in range(len(centroid_querys)):
    splits = centroid_querys[i].query.getTitle().split('_')
    key = splits[1] + '_' + '_'.join(splits[7:12])
    centroid_query_dict[key] = i

#Prepare all the data, this should be optimized in the future.

vdm_metal_dict = {}

for query_id in range(len(centroid_querys)):
    print(query_id)
    query = centroid_querys[query_id].copy()

    metal_dict = {}
    #obtain all member metal coords.
    mem_vdm_names = extract_vdm.get_mem_vdm_names(query)
    for vn in mem_vdm_names:
        splits = vn.split('_')
        if 'centroid' in splits:
            key = splits[1] + '_' + '_'.join(splits[7:12])
        else:
            key = splits[1] + '_' + '_'.join(splits[6:11])

        id = centroid_query_dict[key]  
        metal = metals[id]

        if metal in metal_dict.keys():
            metal_dict[metal] += 1
        else:
            metal_dict[metal] = 1
    vdm_metal_dict[query_id] = metal_dict 



with open(query_dir + 'all_vdm_metal_ASP.tsv', 'w') as f:
    metal_names = ['MN', 'FE', 'CO', 'NI', 'CU', 'ZN']
    f.write('id\tcentroid_metal\ttotal_number\tMN\tFE\tCO\tNI\tCU\tZN\n')
    for query_id in range(len(centroid_querys)):
        metal_dict = vdm_metal_dict[query_id]  
        x = '\t'.join([str(metal_dict[k]) if k in metal_dict.keys() else '0' for k in metal_names])   
        f.write(str(query_id) + '\t' + metals[query_id] + '\t' + str(sum(metal_dict.values())) + '\t' + x + '\n')

# Calc metal count in the vdM lib.
metal_names = ['MN', 'FE', 'CO', 'NI', 'CU', 'ZN']
metal_count_dict = {}
for m in metals:
    if m in metal_count_dict.keys():
        metal_count_dict[m] += 1
    else:
        metal_count_dict[m] = 1
print(metal_count_dict)