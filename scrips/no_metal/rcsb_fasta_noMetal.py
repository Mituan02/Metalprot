# To only download pdbs with sequence similarity < 70%. We need to use cluster the fasta sequences.
# One can run MMseq2 with the downloaded fasta data. Or use the clustered data provided by rcsb.
# The rcsb provide MMseqs2 cluster data at: https://www.rcsb.org/docs/programmatic-access/file-download-services#sequence-data

import os
import prody as pr

### Load the cluster 
workdir = '/mnt/e/DesignData/ligands/NoLigand/'

id_dict = {} # {PDB_ID: PDB_ID extended for fasta}
clu_dict = {} # {PDB_ID ext: clusterID}

with open(workdir + 'bc-30.out', 'r') as f:
    count = 0
    for line in f.readlines():
        count += 1
        entrys = line.split(' ')
        for en in entrys:
            e = en.split('\n')[0]
            if len(e)<4:
                continue
            clu_dict[e] = count
            k = e[0:4]
            if k in id_dict.keys():
                id_dict[k].append(e)
            else:
                id_dict[k] = [e]

#list(id_dict.values())[1:10]  # For manual check

### load targets
targets = []

with open(workdir + 'pdb_list.txt', 'r') as f:
    for line in f.readlines():
        e = line.split('\n')[0]
        if len(e)<4:
            continue
        targets.append(e)

#targets[0:10]

### remove duplicate based on cluster. Each cluster only keep one members.

target_dict = {}  # {PDB_ID: [cluster_id]}
for t in targets:
    if t not in id_dict.keys():
        print(t)
        continue
    ks = id_dict[t]
    for k in ks:
        cluid = clu_dict[k]
        if t in target_dict.keys():
            if cluid not in target_dict[t]:
                target_dict[t].append(cluid)
        else:
            target_dict[t] = [cluid]

#list(target_dict.values())[1:10]

### remove any target that has 2 cluster id. And create the new dict with {clu_id: [targets]}

clu_t_dict = {}  #{clu_id: [targets]}
for t in target_dict.keys():
    if len(target_dict[t]) > 1:
        continue
    id = target_dict[t][0]
    if id in clu_t_dict.keys():
        clu_t_dict[id].append(t)
    else:
        clu_t_dict[id] = [t]

len(list(clu_t_dict))

### Write the clu_t_dict 
with open(workdir + 'reduced.txt', 'w') as f:
    for id in clu_t_dict.keys():
        f.write(str(id) + '\t' + '\t'.join(clu_t_dict[id]) + '\n')

### Download the first represent pdb for each cluster.
def download_pdb(workdir, all_pdbs):
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    exist_pdb = set()
    for file in os.listdir(workdir):
        if file.endswith(".pdb.gz"):
            exist_pdb.add(file.split('.')[0].upper())

    pr.pathPDBFolder(workdir)
    for p in all_pdbs:
        if p in exist_pdb: continue
        pr.fetchPDBviaFTP(p, compressed = False) 


all_pdbs = [pdbs[0] for pdbs in clu_t_dict.values()]

download_pdb(workdir + 'no_metal_pres/', all_pdbs)
