import os
import numpy as np
import itertools
import prody as pr
from scipy.spatial.distance import cdist
from .ligand_database import clu_info
from .search_struct import Query


def read_cluster_info(file_path):
    clu_infos = []
    if not os.path.exists(file_path):
        return clu_infos
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        for line in lines[1:]:
            ts = line[:-1].split('\t')
            if len(ts)!=7: continue
            clu_infos.append(clu_info(Metal=ts[0], clu_type = ts[1], clu_rmsd= ts[2], total_num=float(ts[3]), clu_rank= int(ts[4]), clu_num= int(ts[5]), score=float(ts[6])))
    return clu_infos

def filter_cluter_info(clu_infos, score_cut = 0, clu_num_cut = 10):
    filtered_infos = []
    for info in clu_infos:
        if info.score >= score_cut or info.clu_num >= clu_num_cut:
            filtered_infos.append(info)
    return filtered_infos

def extract_query(workdir, file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10):
    clu_infos = read_cluster_info(workdir + file_path)
    filtered_infos = filter_cluter_info(clu_infos, score_cut, clu_num_cut)

    querys = []
    if len(filtered_infos) == 0:
        return querys

    for info in filtered_infos:
        centroid = [file for file in os.listdir(workdir  + str(info.clu_rank)) if 'centroid' in file][0]
        pbd = pr.parsePDB(workdir + str(info.clu_rank) + '/' + centroid)
        query = Query(pbd, info.score, info.clu_num, info.total_num)
        query.path = workdir + str(info.clu_rank) + '/'
        querys.append(query)

    return querys

def extract_centroid_pdb_in_clu(workdir, file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10):  
    clu_infos = read_cluster_info(workdir + file_path)
    filtered_infos = filter_cluter_info(clu_infos, score_cut, clu_num_cut)
    pdb_paths = []
    if len(filtered_infos) == 0:
        return pdb_paths

    for info in filtered_infos:
        centroid = [file for file in os.listdir(workdir  + str(info.clu_rank)) if 'centroid' in file][0]
        pdb_paths.append(workdir + str(info.clu_rank) + '/' + centroid)
    return pdb_paths


def extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster'], score_cut = 0, clu_num_cut = 2):
    
    subfolders_with_paths = [f.path for f in os.scandir(query_dir) if f.is_dir()]

    querys = []

    for subfolder in subfolders_with_paths:
        exist = True
        for n in file_name_includes:
            if n not in subfolder:
               exist = False
        if exist: 
            qs = extract_query(subfolder + '/', file_path = '_summary.txt', score_cut = score_cut, clu_num_cut = clu_num_cut)
            querys.extend(qs)


    return querys




