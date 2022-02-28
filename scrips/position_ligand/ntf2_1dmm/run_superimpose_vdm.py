import os
import sys
import pdb
import prody as pr
import pickle
import pandas as pd
import numpy as np
from metalprot.basic import constant
from metalprot.basic import utils
from metalprot.combs import search_lig_indep

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_2ndRound/output_F55D_sel/'
file = 'o1_1dmm_876_F55D_tts_657.pdb'
target_path = workdir + file

predefined_win_filters = [15, 19, 27]

outdir = workdir + '1dmm_tts_results/output_superimpose_vdm/'
os.makedirs(outdir, exist_ok=True)
target = search_lig_indep.prepare_rosetta_target(outdir, target_path, predefined_win_filters)
abples, phipsi = utils.seq_get_ABPLE(target)


resnum = 55
cg = 'phenol'
aa = 'ASP'
pos = target.select('resnum ' + str(resnum))
abple = abples[pos.getResindices()[0]]
tag = cg + '_' + aa + '_'

path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'
df_vdm = search_lig_indep.load_old_vdm(path_to_database, cg, aa)
df_vdm_filter = search_lig_indep.filter_db(df_vdm, abple=abple)
grs = df_vdm_filter.groupby(['CG', 'rota', 'probe_name'])
count = 0
for i, (n, gr) in enumerate(grs):
    v = gr.copy()
    ag = search_lig_indep.df2ag(v)
    pr.writePDB(outdir + tag + str(count) + '.pdb', ag)
    count += 1
print(count)


    

    
