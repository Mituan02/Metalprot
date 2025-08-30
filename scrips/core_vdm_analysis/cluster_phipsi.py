'''
Check each cluster, the distribution of delta phi psi compared with the centroid vdM. 
Different clustering method may have different delta phi psi distribution.
'''

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


query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/20210916_2017_2018_selfcenter/pickle_noCys/'


with open(query_dir + 'AAMetalPhiPsi.pkl', 'rb') as f:
    all_querys = pickle.load(f)

with open(query_dir + 'cluster_centroid_dict.pkl', 'rb') as f:
    cluster_centroid_dict = pickle.load(f)

phipsi_dict = {}

for query in all_querys:
    if not ('cluster_0_' in query.query.getTitle() or 'cluster_1_' in query.query.getTitle() or 'cluster_2_' in query.query.getTitle()):
        continue
    deta_phipsi = []
    cluster_key = query.get_cluster_key()
    for id in cluster_centroid_dict[cluster_key].selfcenter_cluster_queryid:
        q = all_querys[id]
        dphi = query.phi - q.phi
        if dphi < -180:
            dphi = dphi + 360
        if dphi > 180:
            dphi = dphi - 360
        dpsi = query.psi - q.psi
        if dpsi < -180:
            dpsi = dpsi + 360
        if dpsi > 180:
            dpsi = dpsi - 360
        deta_phipsi.append((dphi, dpsi))

    phipsi_dict[query.query.getTitle()] = deta_phipsi

len(list(phipsi_dict))

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/reason/'


with open(workdir + 'cluster_phi_psi.tsv', 'w') as f:
    f.write('title\tdphi\tdpsi\n')
    for k in phipsi_dict.keys():
        for xs in phipsi_dict[k]:     
            f.write(k + '\t' + str(round(xs[0], 2)) + '\t' + str(round(xs[1])) + '\n')
