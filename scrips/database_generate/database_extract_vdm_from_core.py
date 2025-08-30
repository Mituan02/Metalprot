'''
For the generation of single vdMers. 
This is the second step. To run this, run the 'database_preparation.py' first. And then run 'run_cluster_single.py'
'''

import os
import sys
import prody as pr
#sys.path.append(r'/mnt/e/GitHub_Design/Metalprot')
from metalprot.database import database_extract
from metalprot.database import database_vdMAtomOrder

'''
python /Users/lonelu/GitHub_Design/Metalprot/scrips/database_generate/database_extract_vdm_from_core.py
'''

def extract_single_vdm(cores, outdir, AA, key, basic = True, extention = None, n = None, key_out = None, phipsi = True):
    for c in cores:   

        if AA == 'CYS':
            ## Filter 4 CYS ones. The reason is that there are too many CYS and I think the all-CYS binding is very special.
            print('------')
            print(list(c.contact_aas.getResnames()))
            if list(c.contact_aas.getResnames()).count('CYS') >= len(list(c.contact_aas.getResnames())):
                print(list(c.contact_aas.getResnames()).count('CYS'))
                continue

        if basic:
            c.generate_AA_Metal(AA, key)
        elif extention:
            c.generate_AA_ext_Metal(AA, extention, key)          
        elif n:
            c.generate_AA_kMetal(AA, n, key, key_out) 
        elif phipsi:
            c.generate_AA_phipsi_Metal(AA, key)
        c.write_vdM(outdir, key = key)
    return


def extract_vdm(coredir, outdir, aas):
    cores = database_extract.load_cores(coredir)

    for AA in aas:
        ### extract core
        extract_single_vdm(cores, outdir + 'AAMetalPhiPsi_' + AA + '_reps/', AA = AA, key = 'AAMetalPhiPsi_' + AA, basic = False, phipsi=True)

        ### shift contacting Oxygen
        if AA == 'ASP' or AA == 'GLU':
            _pdbs = []
            for file in os.listdir(outdir + 'AAMetalPhiPsi_' + AA + '_reps/'):
                #print(file)
                if '.pdb' not in file:
                    continue
                pdb = pr.parsePDB(outdir + 'AAMetalPhiPsi_' + AA + '_reps/' + file)
                database_vdMAtomOrder.asp_glu_oxy_shift(pdb)
                _pdbs.append(pdb)

            _outdir =  outdir + 'AAMetalPhiPsi_' + AA + '_reps_shift/'
            os.makedirs(_outdir, exist_ok=True)
            for pdb in _pdbs:
                pr.writePDB(_outdir + pdb.getTitle(), pdb)
    return

workdir = "/Users/lonelu/DesignData/ligands_metal/ZN_rcsb_datesplit/20211013/_Seq_core_date_reps/"
#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/_Seq_core_2ndshell_date_reps/'
#workdir = '/mnt/e/DesignData/DL/_Seq_core_date_reps/'

outdir = '/Users/lonelu/DesignData/ligands_metal/ZN_rcsb_datesplit/20211013/20231020_cys_vdm_filtered/'
#outdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20220116_2ndshell/20220128_1stshell/'
#outdir = '/mnt/e/DesignData/DL/20220511_vdm_reps/'
os.makedirs(outdir, exist_ok = True)

#aas = ['HIS', 'GLU', 'ASP', 'CYS']
aas = ['CYS']

extract_vdm(workdir, outdir, aas)