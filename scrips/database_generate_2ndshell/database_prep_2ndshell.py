import os
import sys
import prody as pr
import shutil
#sys.path.append(r'/mnt/e/GitHub_Design/MetalDesign')
from metalprot.database import database_extract

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/"
metal_sel = 'name ZN'

# To get the core reps with 2nd shell, extract the core with 2nd shell first. Then match the name from the 1st shell core reps. 

cores_2nd = database_extract.extract_all_core_seq_from_path('/mnt/e/DesignData/ligands/ZN_rcsb/all_rcsb/', metal_sel, extend = 4, extract_2ndshell = True, _2ndshell_extend = 4)
os.makedirs(workdir + '20220116_2ndshell/_Seq_core_2ndshell/', exist_ok=True)
database_extract.superimpose_core_and_writepdb(cores_2nd, cores_2nd[0], metal_sel, workdir + '20220116_2ndshell/_Seq_core_2ndshell/')

# Check datesplit.py to add the date.


### Remove duplicate by 2nd shell.

metal_sel = 'name NI MN ZN CO CU MG FE' 

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/"

database_extract.rm_core(workdir + '_Seq_core_2ndshell_date/', workdir + '_Seq_core_2ndshell_date_filter/', min_contact_aa_num=2)

pdbs = database_extract.get_all_pbd_prody(workdir + '_Seq_core_2ndshell_date_filter/')

clusters = database_extract.reduce_dup(pdbs, metal_sel)

database_extract.write_dup_summary(workdir, pdbs, clusters)

database_extract.extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '_Seq_core_2ndshell_date_reps/')




#>>> Extract 2ndshell from all metal.
import os
from metalprot.database import database_extract

def extract_core_2nd(indir, outdir, metal_sel):
    #indir = '/mnt/e/DesignData/ligands/CO_rcsb/all_rcsb/'
    #metal_sel = 'name CO'
    cores_2nd = database_extract.extract_all_core_seq_from_path(indir, metal_sel, extend = 4, extract_2ndshell = True, _2ndshell_extend = 4)
    database_extract.superimpose_core_and_writepdb(cores_2nd, cores_2nd[0], metal_sel, outdir)
    return

def extrat_mvdm2nd_reps(indir, outdir_filter, outdir_reps, metal_sel):

    database_extract.rm_core(indir, outdir_filter, min_contact_aa_num=2)

    pdbs = database_extract.get_all_pbd_prody(outdir_filter)

    clusters = database_extract.reduce_dup(pdbs, metal_sel)

    database_extract.write_dup_summary(outdir_reps, pdbs, clusters)

    database_extract.extract_rep_and_writepdb(pdbs, clusters, metal_sel, outdir_reps)

    return

def run_ext_all_mvdm2ndshell():
    
    workdir = '/mnt/e/DesignData/ligands/'

    mts = ['CO', 'CU', 'FE', 'MN', 'NI']

    for mt in mts:
        metal_sel = 'name ' + mt
        for folder in os.listdir(workdir + mt + '_rcsb/'):
            if not 'all_rcsb' in folder or not os.path.isdir(workdir + mt + '_rcsb/' + folder):
                continue
            indir = workdir + mt + '_rcsb/' + folder + '/'
            outdir = workdir + 'all/20220710/_Seq_core_2ndshell_' + mt + '/'
            os.makedirs(outdir, exist_ok=True)
            extract_core_2nd(indir, outdir, metal_sel)

        outdir_filter = workdir + 'all/20220710/_Seq_core_2ndshell_filter_' + mt + '/'
        os.makedirs(outdir_filter, exist_ok=True)
        outdir_reps = workdir + 'all/20220710/_Seq_core_2ndshell_reps_' + mt + '/'
        os.makedirs(outdir_reps, exist_ok=True)

        extrat_mvdm2nd_reps(outdir, outdir_filter, outdir_reps, metal_sel)
    return

run_ext_all_mvdm2ndshell()