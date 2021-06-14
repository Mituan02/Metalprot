import os
import sys
from metalprot.apps.core import Core
import prody as pr
from metalprot import search_struct, extract_vdm
from metalprot import ligand_database as ldb

### Set working directories.

workdir = '/mnt/e/DesignData/ligands/LigandBB/zn_bench_mark/output_test_7cid_vdm_analyze/'

outdir = workdir + 'output_test_7cid_vdm_analyze/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

#target_path = workdir + '3into4_helix_assembly_renum.pdb'

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'


### Get query pdbs 

querys = extract_vdm.extract_all_centroid(query_dir, summary_name = '_summary.txt', file_name_includes = ['cluster', 'M1-1_AAMetalSc'], score_cut = 0, clu_num_cut = 1)

for q in querys:
    print(q.query.getResnames()[0] + '-' + str(len(q.query.select('heavy'))))
    #print(q.query.getNames())
print(len(querys))
### Get target VdMs

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

_pdb = ldb.get_all_pbd_prody(workdir)

cores = ldb.extract_all_core_seq_from_path(workdir, metal_sel, extend = 4, extract_2ndshell = True)

ldb.writepdb(cores, outdir + 'target_core/')

_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core/')

all_aas = ['HIS', 'ASP', 'GLU', 'CYS']
for aa in all_aas:
    for p in _pdbs:
        core = Core(p)
        core.generate_AA_Metal(aa, key = 'AAMetal_' + aa)
        core.write_vdM(outdir + 'target_core_aa/', key = 'AAMetal_' + aa)

### Match target 

_pdbs = ldb.get_all_pbd_prody(outdir + 'target_core_aa/')

pairs = []

for _pdb in _pdbs:
    for query in querys:
        if len(_pdb.select('heavy')) == len(query.query.select('heavy')):
            _pdb_copy = _pdb.copy()
            print('----------')
            print(_pdb.getResnames()[0] + '-' + str(len(_pdb.select('heavy'))) )
            print(query.query.getTitle() + '-' + str(len(query.query.select('heavy'))))
            pr.calcTransformation(_pdb_copy.select('heavy'), query.query.select('heavy')).apply(_pdb_copy)
            rmsd = pr.calcRMSD(_pdb_copy.select('heavy'), query.query.select('heavy'))

            if rmsd <= 1:
                pairs.append((_pdb_copy, query, rmsd))

count = 0
for pair in pairs:
    pr.writePDB(outdir + str(count) + '_' + pair[0].getTitle(), pair[0])
    pr.writePDB(outdir + str(count) + '_' + pair[1].query.getTitle(), pair[1].query)
    count += 1

with open(outdir + '_summary.txt', 'w') as f:
    #f.write('Metal\tclu_type\tclu_rmsd\ttotal_num\tclust_rank\tclust_num\tscore\n')
    for pair in pairs:
        f.write(str(pair[2]) + '\t' + pair[0].getTitle() + '\t' + pair[1].to_tab_string() + '\n')  