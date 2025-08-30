'''
Calc pairwise clash dist of core pdb. 
'''
import os
import prody as pr
from metalprot.database import database_extract, database_evaluate
import numpy as np
import matplotlib.pyplot as plt

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211008/_Seq_core_filter/'

'''
pdb = pr.parsePDB(workdir + '1982_4cpa_ZN_1.pdb')

min_dists = database_evaluate.calc_bb_clash_min_dist(pdb)

#min_dists = [3.386951136346675, 3.363647127746903, 3.076913713447294, 3.2841622980601897]
'''
pdbs = database_extract.get_all_pbd_prody(workdir)

all_title = []
all_min_dist = []

for pdb in pdbs:
    min_dists = database_evaluate.calc_bb_clash_min_dist(pdb)
    titles = [pdb.getTitle() for i in range(len(min_dists))]
    all_title.extend(titles)
    all_min_dist.extend(min_dists)


outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211008/reason/'

with open(outdir + 'bb_clash_min_dist.tsv', 'w') as f:
    for i in range(len(all_title)):
        f.write(all_title[i] + '\t' + str(all_min_dist[i])+ '\n')


def plt_dist_his(dist, outplot_name = 'dist.png'):
    n, bins, patches = plt.hist(x=dist, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.2)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('dist')
    plt.ylabel('Frequency')
    plt.title('vdM sc&bb clash Dist Histogram')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outplot_name)
    plt.close()

plt_dist_his(all_min_dist, outdir + 'bb_clash_min_dist.png')

### Extract vdM sc clash min dist.

all_title = []
all_min_dist = []

for pdb in pdbs:
    min_dists = database_evaluate.calc_vdm_clash_min_dist(pdb)
    titles = [pdb.getTitle() for i in range(len(min_dists))]
    all_title.extend(titles)
    all_min_dist.extend(min_dists)

with open(outdir + 'vdM_sc_clash_min_dist.tsv', 'w') as f:
    for i in range(len(all_title)):
        f.write(all_title[i] + '\t' + str(all_min_dist[i])+ '\n')

plt_dist_his(all_min_dist, outdir + 'vdM_sc_clash_min_dist.png')