import os
import prody as pr
from metalprot import ligand_database



### load all the csv
def load_csv(workdir):
    all_lines = []
    for file in os.listdir(workdir):
        if file.endswith(".csv"):
            with open(workdir + file, 'r') as f:
                all_lines.extend(f.readlines())

    with open(workdir + 'all_rcsb.tsv', 'w') as f:
        f.write('\t'.join(all_lines[0].split(',')))
        for r in all_lines:
            if 'Entry ID' not in r:
                f.write('\t'.join(r.split(',')))
    return all_lines


### only keep the ligand contain porphyrin
def extract_pdb_ligand(workdir, all_lines):
    pdbs = {}

    key = ''
    porphyrin = ''
    porphyrinId = ''
    for r in all_lines:
        if 'Entry ID' in r:
            continue
        items =  r.split(',')
        #print(r)
        if items[0]:
            key = items[0][1:-1]

        if len(items) >= 5 and 'PORPHYRIN' in items[4]:
            porphyrin = items[4][1:-1]
            porphyrinId = items[3][1:-1]

            if key in pdbs.keys():
                pdbs[key].append((porphyrinId, porphyrin))
            else:
                pdbs[key] = [(porphyrinId, porphyrin)]

    with open(workdir + 'all_porphyrin.tsv', 'w') as f:
        f.write('Entry ID\tligand\tliand name\n')
        for k in pdbs.keys():
            for x in pdbs[k]:
                f.write(k + '\t' + x[0] + '\t' + x[1] + '\n')

    return pdbs




workdir = '/mnt/e/DesignData/ligands/porphyrin/porphyrin_search/'

all_lines = load_csv(workdir)

pdbs = extract_pdb_ligand(all_lines)

list(pdbs)[0:10]
pdbs['5JPR']

### 
workdir = '/mnt/e/DesignData/ligands/porphyrin/pdbs/'

all_lines = load_csv(workdir)

ligands = ['HEM', 'HNI', 'COH', 'HEB', 'FDE', 'ZNH', 'HEC', 'HEA', 'HAS', 'MNH', 'MNR']
    

### only keep the ligand contain porphyrin
def extract_pdb_ligand2(workdir, all_lines, ligands):

    pdbs = {}

    key = ''
    porphyrin = ''
    porphyrinId = ''
    for r in all_lines:
        if 'Entry ID' in r:
            continue
        items =  r.split(',')
        #print(r)
        if items[0]:
            key = items[0][1:-1]

        if len(items) >= 8 and items[7][1:-1] in ligands:
            porphyrin = items[8][1:-1]
            porphyrinId = items[7][1:-1]

            if key in pdbs.keys():
                pdbs[key].append((porphyrinId, porphyrin))
            else:
                pdbs[key] = [(porphyrinId, porphyrin)]

    with open(workdir + 'all_porphyrin.tsv', 'w') as f:
        f.write('Entry ID\tligand\tliand name\n')
        for k in pdbs.keys():
            for x in pdbs[k]:
                f.write(k + '\t' + x[0] + '\t' + x[1] + '\n')

    return pdbs

pdbs = extract_pdb_ligand2(workdir, all_lines, ligands)

exist_ligands = ['HEM', 'HNI', 'FDE', 'ZNH', 'HEC', 'HEA', 'HAS', 'MNH', 'MNR']

### only download pdbs with sequence similarity less than 70%
targets = set()
with open(workdir + 'all_porphyrin.tsv', 'r') as f:
    for line in f.readlines():
        e = line.split('\t')[0]
        if len(e)!=4:
            continue
        targets.add(e)


id_dict = {} # {PDB_ID: PDB_ID extended for fasta}
clu_dict = {} # {PDB_ID ext: clusterID}

with open(workdir + 'bc-70.out', 'r') as f:
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

list(id_dict.keys())[1:10] 
list(id_dict.values())[1:10]  # For manual check

list(clu_dict.keys())[1:10] 
list(clu_dict.values())[1:10]  # For manual check


### remove duplicate based on cluster. Each cluster only keep one members.

no_in_id_dict = []

target_dict = {}  # {PDB_ID: [cluster_id]}
for t in targets:
    if t not in id_dict.keys():
        print(t)
        no_in_id_dict.append(t)
        continue
    ks = id_dict[t]
    for k in ks:
        cluid = clu_dict[k]
        if t in target_dict.keys():
            if cluid not in target_dict[t]:
                target_dict[t].append(cluid)
        else:
            target_dict[t] = [cluid]


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


all_pdbs = [pdbs[0] for pdbs in clu_t_dict.values()]

all_pdbs.extend(no_in_id_dict)


with open(workdir + 'reduced.txt', 'w') as f:
    for pdb in all_pdbs:
        f.write(pdb +  '\n')

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


download_pdb(workdir + 'reduced_pdbs/', all_pdbs)


download_pdb(workdir + 'all_pdbs/', targets)