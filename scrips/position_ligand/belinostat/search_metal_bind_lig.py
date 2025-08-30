'''
Search the mvdms on the bb, then find the pre-posed lig that are close to the metal.
'''
import os
import prody as pr

workdir = '/mnt/e/DesignData/Metalloenzyme/'

lig_indir = workdir + 'outdir/'

metal_indir = workdir + 'targets/output_01_f63440_nick_ala__20220616-225348/represents/'

metal_cands = []
for file in os.listdir(metal_indir):
    if not '.pdb' in file or 'idealgeo' in file:
        continue
    metal_cand = pr.parsePDB(metal_indir + file)
    metal_cands.append(metal_cand)

lig_cands = []
for file in os.listdir(lig_indir):
    if not '.pdb' in file or 'vdm' in file:
        continue
    lig_cand = pr.parsePDB(file)
