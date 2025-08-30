'''
It turned out the pair vdm strategy is slow and not working now!
Move back to the old method. Search the metal binding first, then attach the ligand on metal with all possible conformations. Then search the vdMs.
The script here is to pre generate ligands attached on metals of all conformations.
It first superimpose the ligand on X-axis, then rotate by Y, then rotate by Z.
The final result is dumped into pickle file for quick loading.
'''
import os
import sys
import prody as pr
import numpy as np
import math
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
import pickle

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/prepare_belinostat_on_metal.py
'''

def _generate_rots_byX(outdir, std, std_cp_sel = 'all', tag_pre = '50g_md_0_x'):
    '''
    Rotate by x_axis.
    '''
    #>>> lig to x-axis with ZN on [0, 0, 0]
    dist = pr.calcDistance(std.select('name ZN'), std.select('name O4'))[0]
    xcoord = np.array([[0, 0, 0], [dist , 0, 0]])
    tf = pr.calcTransformation(std.select('name O4 ZN').getCoords()[::-1], xcoord)
    tf_rv = pr.calcTransformation(xcoord, std.select('name O4 ZN').getCoords()[::-1])
    ligs = []
    for i in range(120):
        std_cp = std.select(std_cp_sel).copy()
        tf.apply(std_cp)
        v = std_cp.getCoords()
        theta = 2*math.pi/120 * i
        axis = xcoord[1]/norm(xcoord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        #new_v = tf_rv.apply(new_v) 
        std_cp.setCoords(new_v)
        std_cp.setTitle(tag_pre + str(i))
        #pr.writePDB(outdir + tag_pre + str(i) +'.pdb', std_cp)
        ligs.append(std_cp)
    return ligs


def _generate_rots_byXYZ(outdir, std, _coord, tag_pre = '50g_md_0_x_y'):
    '''
    After rotate by X, the ligand can be shift a little bit
    '''
    ligs = []
    for i in range(-5, 6):
        std_cp = std.copy()
        v = std_cp.getCoords()
        theta = 2*math.pi/120 * i
        axis = _coord[1]/norm(_coord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        #new_v = tf_rv.apply(new_v) 
        std_cp.setCoords(new_v)
        std_cp.setTitle(tag_pre + str(i))
        #pr.writePDB(outdir + tag_pre + str(i) + '.pdb', std_cp)
        ligs.append(std_cp)
    return ligs


def main(name = '50g_md_0.pdb'):
    workdir = '/mnt/e/DesignData/Metalloenzyme/belinostat/ligs/meo_50g_amber14eht_md_out/'   
    #>>> Please change name for each ligand.
    lig = pr.parsePDB(workdir + name)

    outdir = workdir + lig.getTitle() + '_std_ligs/'
    #os.makedirs(outdir, exist_ok = True)

    tag_pre = lig.getTitle() + '_x'
    ligs_x = _generate_rots_byX(outdir, lig, tag_pre = tag_pre)

    ligs_xy = []
    for std in ligs_x:
        #std = pr.parsePDB(outdir+ '50g_md_0_x0.pdb')
        tag_pre = std.getTitle() + '_y'
        _coord = np.array([[0, 0, 0], [0, 1, 0]])
        _ligs = _generate_rots_byXYZ(outdir, std, _coord, tag_pre = tag_pre)
        ligs_xy.extend(_ligs)

    ligs_xyz = []
    ligs_xyz.append(ligs_x[0]) # Make sure the first one is the starting one. 
    for std in ligs_xy:
        #std = pr.parsePDB(outdir+ '50g_md_0_x0.pdb')
        tag_pre = std.getTitle() + '_z'
        _coord = np.array([[0, 0, 0], [0, 0, 1]])
        _ligs = _generate_rots_byXYZ(outdir, std, _coord, tag_pre = tag_pre)
        ligs_xyz.extend(_ligs)

    all_coords = []
    labels = []

    for i in range(len(ligs_xyz)):
        coords = ligs_xyz[i].select('heavy').getCoords()
        all_coords.extend(coords)
        labels.append(ligs_xyz[i].getTitle())

    with open(workdir + lig.getTitle() + '_coord.pkl', 'wb') as f:
        pickle.dump((ligs_x[0], all_coords, labels), f)

    return

for i in range(17):
    main(name = '50g_md_' + str(i) + '.pdb')

'''
import os
import pickle
import prody as pr

workdir = '/mnt/e/DesignData/Metalloenzyme/belinostat/ligs/meo_50g_amber14eht_md_out/'  
with open(workdir + '50g_md_0.pkl', 'rb') as f:
    all_ligs = pickle.load(f)
lig = all_ligs[0]
outdir = workdir + lig.getTitle() + '_std_ligs1/'
os.makedirs(outdir, exist_ok = True)
[pr.writePDB(outdir + _lig.getTitle(), _lig) for _lig in all_ligs[0:5]]

with open(workdir + '50g_md_1.pkl', 'rb') as f:
    all_ligs = pickle.load(f)
lig = all_ligs[0]
outdir = workdir + lig.getTitle() + '_std_ligs1/'
os.makedirs(outdir, exist_ok = True)
[pr.writePDB(outdir + _lig.getTitle(), _lig) for _lig in all_ligs[0:5]]
len(all_ligs)
'''