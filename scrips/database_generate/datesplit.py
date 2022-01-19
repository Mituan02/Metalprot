'''
Split the metal vdM database into 2 folder by date. (before ~ 2017) (2018 ~ after)
'''

import os
import sys
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
#from metalprot.database import database_extract


#----------------------------------------------------------
#extract core with date attached in the front as the name. Here, I just add the date in the existed cores.
workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/"

filename='all_rcsb.txt'

pdb_date = {}

if not os.path.exists(workdir):
    os.mkdir(workdir)
with open(workdir + filename, 'r') as f:
    for line in f.readlines():
        r = line.split('\t')
        if r[0] == '"Entry ID"': continue
        if r[0] == '' or r[6]== '' or (',' in r[6]) or float(r[6].split('"')[1]) > 2.5:
            continue
        pdb_date[r[0].split('"')[1]] = r[3].split('"')[1].split('-')[0]


workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell/"
outdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date/"
os.makedirs(outdir, exist_ok=True)

for file in os.listdir(workdir):
    if file.endswith(".pdb"):
        name = file.split('_')[0].upper()
        date = pdb_date[name]

        shutil.copy(workdir + file, outdir + date + '_' + file)


#---------------------------------------------------------
#Run database_preparation script to remove duplicates ordered by date. 
#Split the cores by date: 1. _Seq_core_date2017_reps. 2. _Seq_core_date2018_reps
#Cores in _Seq_core_date2017_reps are used for vdM extraction and clustering. 


#----------------------------------------------------------
#Copy all the original full pdbs of _Seq_core_date2018_reps into testing folder. 

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/"

test_pdbs = set()

for file in os.listdir(workdir + '_Seq_core_date2017_reps/'):
    if file.endswith(".pdb"):
        name = file.split('_')[1].lower()
        test_pdbs.add(name)

origi_file_dir = "/mnt/e/DesignData/ligands/ZN_rcsb/all_rcsb2/"
for file in os.listdir(origi_file_dir):
    print(file)
    if file.split('.')[0] in test_pdbs:
        print('---')
        shutil.copy(origi_file_dir + file, workdir + '_train_full_pdbs/' + file)


#----------------------------------------------------------
#Copy all the original full pdbs of _Seq_core_date2018_reps into testing folder. 
