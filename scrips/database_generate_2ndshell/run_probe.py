import os 
import sys
sys.path.append(r'/wynton/home/degradolab/lonelu/GitHub_Design/Combs2')
import combs2

path_to_probe = '/wynton/home/degradolab/lonelu/software/molprobity/molprobity/bin/linux/probe'


db_H_dir = '/wynton/home/degradolab/lonelu/DesignData/_Seq_core_2ndshell_date_reps_withH/'
outdir = '/wynton/home/degradolab/lonelu/DesignData/probe_out/'
os.makedirs(outdir, exist_ok = True)


for pdb in os.listdir(db_H_dir):
    if not pdb.endswith('.pdb'):
        continue
        
    x = combs2.design.probe.parse_probe(db_H_dir + pdb, path_to_probe = path_to_probe, outdir = outdir)

    x.to_csv(outdir + pdb + '_probe.csv')

    '''
    cmd = [path_to_probe,'-U -SEGID -CON -NOFACE', '-Explicit', '-MC', '-WEAKH -DE32', '-' + str(4), '-SE', '"', 'ALL', 'NOT METAL', '"', db_H_dir + pdb]
    os.system(' '.join(cmd) + ' > ' + outdir + pdb.split('.')[0])
    '''