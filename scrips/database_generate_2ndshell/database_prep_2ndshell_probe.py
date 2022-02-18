'''
Prepare 2nd shell vdm from probe calculated Hbond. 
We need the probe output and the metal binding core to run the function.
After run the script, one can use the 'run_cluster_2ndshell_probe.py' to get the clustered database.
*hb* (H-bond), then *wh* (weak H-bond), then *cc* (close contact).  Contacts
    labeled *so* (small overlap), *bo* (bad overlap), *wc* (wide contact)
'''

from itertools import chain
import os
import pandas as pd
import prody as pr
from metalprot.basic import probe
from metalprot.database import database_extract

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/database_generate/database_prep_2ndshell_probe.py
'''

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"

probe_dir = workdir + '_Seq_core_2ndshell_date_reps_probe/'

outdir = workdir + '_Seq_core_2ndshell_date_reps_probe2ndshell/'
os.makedirs(outdir , exist_ok=True)

df_probes = {}
for probefile in os.listdir(probe_dir):
    if '.probe' not in probefile or '.csv' in probefile:
        continue
    df = probe.probe2csv(probe_dir, probefile)
    df_probes[probefile.split('.')[0]] = df


cores = database_extract.load_cores(workdir + '_Seq_core_2ndshell_date_reps/') 

# core = cores[2]
# df = df_probes[core.full_pdb.getTitle()]
# df_hbs = extract_2ndshell_probe(core, df)

for core in cores:
    name = core.full_pdb.getTitle()
    if name not in df_probes.keys():
        continue
    df = df_probes[name]

    df_hbs = probe.extract_2ndshell_probe(core, df)
    df_hbs.to_csv(outdir + name + '_probe2ndshell.csv')

    for i in range(0, len(df_hbs)):
        try:
            chain1 = df_hbs.iloc[i]['chain1']
            resnum1 = int(df_hbs.iloc[i]['resnum1'])
            chain2 = df_hbs.iloc[i]['chain2']
            resnum2 = int(df_hbs.iloc[i]['resnum2'])
            if chain1 == chain2 and resnum1 == resnum2:
                print(name + ' 2nd shell should not happen in the same position.')
                continue
            resind1 = core.full_pdb.select('chid ' + chain1 + ' and resnum ' + str(resnum1)).getResindices()[0]
            resind2 = core.full_pdb.select('chid ' + chain2 + ' and resnum ' + str(resnum2)).getResindices()[0]

            _2ndshellVdm = core.full_pdb.select('resindex ' + str(resind1)  + ' ' + str(resind2) + ' ' + str(core.metal_resind))

            is_bb_or_sc = '_sc_'
            if df_hbs.iloc[i]['name2'] == 'N' or df_hbs.iloc[i]['name2'] == 'O' or df_hbs.iloc[i]['name2'] == 'H':
                is_bb_or_sc = '_bb_'

            tag = is_bb_or_sc +  '_'.join([str(x) for x in df_hbs.iloc[i]])

            resname = df_hbs.iloc[i]['resname1']
            outdir_aa = outdir + resname + '/'
            os.makedirs(outdir_aa, exist_ok=True)
            pr.writePDB(outdir_aa + name + tag + '.pdb', _2ndshellVdm)
        except:
            print(name)
            print(df_hbs)