import os
import sys
import prody as pr
from metalprot import search_struct, extract_vdm
from metalprot import ligand_database as ldb

### Set working directories.

workdir = '/mnt/e/DesignData/ligands/LigandBB/Design_Matt/New/Metal/'

outdir = workdir + 'output_test/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

#target_path = workdir + '3into4_helix_assembly_renum.pdb'

query_dir = '/mnt/e/DesignData/ligands/FE_rcsb/'


### Get query pdbs 

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M1-1'], score_cut = 0, clu_num_cut = 1)

### Get target VdMs

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

_pdbs = ldb.get_all_pbd_prody(workdir)

cores = ldb.extract_all_core_seq_from_path(workdir, metal_sel, extend = 4, extract_2ndshell = True)

ldb.writepdb(cores, outdir + 'target_core/')

_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core/')

all_aas = ['HIS', 'ASP', 'GLU', 'CYS']
for aa in all_aas:
    aa_cores = ldb.extract_all_core_aa(_pdbs, metal_sel, aa = 'resname '+ aa, consider_phipsi = False, extention_out = 0, extract2aa = False, extention = 0, filter_aas = False, aas = None, extract2aa_sep = False, extract2ndshell = False)
    ldb.writepdb(aa_cores, outdir + 'target_core_aa/')

### Match target 

_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core_aa/')

pairs = []

for _pdb in _pdbs:
    for query in querys:
        if len(_pdb.select('heavy')) != len(query.query.select('heavy')):
                continue
        pr.calcTransformation(_pdb.select('heavy'), query.query.select('heavy')).apply(_pdb)
        rmsd = pr.calcRMSD(_pdb.select('heavy'), query.query.select('heavy'))

        if rmsd <= 0.5:
            pairs.append((_pdb, query))
            continue

for pair in pairs:
    pr.writePDB(outdir + pair[0].getTitle(), pair[0])
    pr.writePDB(outdir + pair[1].query.getTitle(), pair[1].query)

with open(outdir + '_summary.txt', 'w') as f:
    #f.write('Metal\tclu_type\tclu_rmsd\ttotal_num\tclust_rank\tclust_num\tscore\n')
    for pair in pairs:
        f.write(pair[0].getTitle() + pair[1].to_tab_string() + '\n')  