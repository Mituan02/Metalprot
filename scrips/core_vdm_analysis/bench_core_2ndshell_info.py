import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import pickle

database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220116_selfcenter_all2ndshell/pickle_noCYS/'

with open(database_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)

with open(database_dir + 'all_2ndshell_vdm_id_dict.pkl', 'rb') as f:
    all_2ndshell_vdm_id_dict = pickle.load(f)

#Transform the key so that the bench core could be indexed.
list(all_2ndshell_vdm_id_dict.keys())[0]
vdm_id_dict = {}
for key in all_2ndshell_vdm_id_dict.keys():
    newkeysStr = key.split('_centroid_')[1].split('_')
    newkey = '_'.join([newkeysStr[i] for i in [0, 1, 2, 3, 5, 6]])
    vdm_id_dict[newkey] = all_2ndshell_vdm_id_dict[key]


bbdatabase_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220116_selfcenter_all2ndshell/pickle_noCYS/'

with open(bbdatabase_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    bb_2ndshell_vdms = pickle.load(f)

with open(bbdatabase_dir + 'all_2ndshell_vdm_id_dict.pkl', 'rb') as f:
    bb_2ndshell_vdm_id_dict = pickle.load(f)

#Transform the key so that the bench core could be indexed.
list(bb_2ndshell_vdm_id_dict.keys())[0]
bbvdm_id_dict = {}
for key in bb_2ndshell_vdm_id_dict.keys():
    newkeysStr = key.split('_centroid_')[1].split('_')
    newkey = '_'.join([newkeysStr[i] for i in [0, 1, 2, 3, 5, 6]])
    bbvdm_id_dict[newkey] = bb_2ndshell_vdm_id_dict[key]



workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_3contact_noCYS/"

cores = database_extract.load_cores(workdir)

### extract reps and cluster
with open(workdir + '_summary_2ndshell3.tsv', 'w') as f:
    f.write('Title\tAAs\tTotal2nd\tTotalbb2nd\tEmpty\tTotal2ndScore\tTotalbb2ndScore\n')

    for c in cores:
        c.atomGroupDict.clear()
        c.generate_AA_2ndShell_Metal(key = 'all_2ndshell', filter_AA =False, AA = '', only_bb_2ndshell = False)
        c.generate_AA_2ndShell_Metal(key = 'all_bb2ndshell', filter_AA =False, AA = '', only_bb_2ndshell = True)

        namekeys = []
        notempty = set()

        if 'all_2ndshell' in c.atomGroupDict.keys():
            for ag in c.atomGroupDict['all_2ndshell']:
                #ag[1] here is supposed to be tag in the title.
                nkey = c.full_pdb.getTitle() + ag[1][:-1]
                if not nkey in vdm_id_dict.keys():
                    continue
                notempty.add(ag[1].split('_')[1])
                namekeys.append(nkey)

        total2ndshell = len(namekeys)
        if 'all_bb2ndshell' in c.atomGroupDict.keys():
            totalbb2ndshell = len(c.atomGroupDict['all_bb2ndshell'])
        else:
            totalbb2ndshell = 0
        empty2ndshell = 3- len(notempty)
        total2ndshellScore = 0
        for nkey in namekeys:
            if not nkey in vdm_id_dict.keys():
                continue
            vdm_id = vdm_id_dict[nkey]
            total2ndshellScore += all_2ndshell_vdms[vdm_id].score

        bbnamekeys = []
        if 'all_bb2ndshell' in c.atomGroupDict.keys():
            for ag in c.atomGroupDict['all_bb2ndshell']:
                #ag[1] here is supposed to be tag in the title.
                bbnkey = c.full_pdb.getTitle() + ag[1][:-1]
                if not nkey in bbvdm_id_dict.keys():
                    continue
                bbnamekeys.append(bbnkey) 
        totalbb2ndshellScore = 0
        for bbnkey in bbnamekeys:
            if not bbnkey in bbvdm_id_dict.keys():
                continue
            vdm_id = bbvdm_id_dict[bbnkey]
            totalbb2ndshellScore += bb_2ndshell_vdms[vdm_id].score                   

        f.write(c.full_pdb.getTitle() + '\t')
        f.write('_'.join([a.getResname() for a in c.contact_aas]) + '\t')
        f.write(str(total2ndshell) + '\t')
        f.write(str(totalbb2ndshell) + '\t')
        f.write(str(empty2ndshell) + '\t')
        f.write(str(total2ndshellScore) + '\t')
        f.write(str(totalbb2ndshellScore) + '\n')
    

### Benchmark 3 contact core 2ndshell from probe.

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
import pickle
import pandas as pd

database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/pickle_noCYS/'

with open(database_dir + 'all_2ndshell_vdms.pkl', 'rb') as f:
    all_2ndshell_vdms = pickle.load(f)

with open(database_dir + 'all_2ndshell_vdm_id_dict.pkl', 'rb') as f:
    all_2ndshell_vdm_id_dict = pickle.load(f)

#Transform the key so that the bench core could be indexed.
list(all_2ndshell_vdm_id_dict.keys())[0]
vdm_id_dict = {}
for key in all_2ndshell_vdm_id_dict.keys():
    newkeysStr = key.split('_centroid_')[1]
    vdm_id_dict[newkeysStr] = all_2ndshell_vdm_id_dict[key]
list(vdm_id_dict.keys())[0]

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_3contact_noCYS/"
files = []
for file in os.listdir(workdir):
    if '.pdb' not in file:
        continue
    files.append(file)


infos = []

database_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps_probe2ndshell/probe_infos/'
for file in files:
    name = file.split('.')[0]
    df_hbs = pd.read_csv(database_dir + name + '_probe2ndshell.csv', header=0,index_col=0)

    total_bb = 0
    score = 0
    seen = set()
    for i in range(len(df_hbs)):
        chain1 = df_hbs.iloc[i]['chain1']
        resnum1 = int(df_hbs.iloc[i]['resnum1'])
        chain2 = df_hbs.iloc[i]['chain2']
        resnum2 = int(df_hbs.iloc[i]['resnum2'])
        is_bb_or_sc = '_sc_'
        if df_hbs.iloc[i]['name2'] == 'N' or df_hbs.iloc[i]['name2'] == 'O' or df_hbs.iloc[i]['name2'] == 'H':
            is_bb_or_sc = '_bb_'
            total_bb += 1
        tag = is_bb_or_sc +  '_'.join([str(x) for x in df_hbs.iloc[i]])
        key = name + tag

        seen.add((chain1, resnum1))
        try:
            id = vdm_id_dict[key]
            vdm = all_2ndshell_vdms[id]
            score += vdm.score
        except:
            print(key + ' is not in vdm_id_dict.')

        empty = 3 - len(seen)
        total = len(df_hbs)

    infos.append((name, total, total_bb, empty, score))


with open(workdir + '_probeHB_summary.tsv', 'w') as f:
    f.write('Title\tTotal2nd\tTotalbb2nd\tEmpty\tTotal2ndScore\n')

    for info in infos:
        f.write('\t'.join([str(x) for x in info]) + '\n')