import os
import prody as pr


# Manipulate rcsb file

def organize_rcsb_file(workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"):
    '''
    The .csv files downloaded from rcsb database will be combined first, 
    then generate tab deliminated txt file.
    '''
    all_lines = []
    for file in os.listdir(workdir):
        if file.endswith(".csv"):
            with open(workdir + file, 'r') as f:
                all_lines.extend(f.readlines())
    with open(workdir + 'all_rcsb.txt', 'w') as f:
        f.write('\t'.join(all_lines[0].split(',')))
        for r in all_lines:
            if 'Entry ID' not in r and r.split(',')[0]!= '':
                f.write('\t'.join(r.split(',')))

# download rcsb pdb files

def download_pdb(workdir, filename, resolution = 2.5):
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    all_pdbs = []
    with open(workdir + filename, 'r') as f: 
        for line in f.readlines():
            #print(line)
            r = line.split('\t')
            #print(r)
            if r[0] == '"Entry ID"': continue
            if r[0] == '' or r[4]== '' or (',' in r[4]) or float(r[4].split('"')[1]) > 2.5:
                continue
            all_pdbs.append(r[0].split('"')[1])

    exist_pdb = set()
    for file in os.listdir(workdir):
        if file.endswith(".pdb.gz"):
            exist_pdb.add(file.split('.')[0].upper())

    pr.pathPDBFolder(workdir)
    for p in all_pdbs:
        if p in exist_pdb: continue
        pr.fetchPDBviaFTP(p, compressed = False) 

    # #Then unzip them in linux with:  
    # cd /mnt/e/DesignData/ligands/NI_rcsb/
    # gunzip *.gz

