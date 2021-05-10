import os
import sys
import prody as pr
import numpy as np
import math as m
import string
from . import constant

def Rx(theta):
    return np.matrix([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])
  
def Ry(theta):
    return np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]])
  
def Rz(theta):
    return np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

def getRotation(phi, theta, psi):
    # phi = #m.pi/2*i
    # theta = m.pi/2*i
    # psi = #m.pi/2
    # print("phi =", phi)
    # print("theta  =", theta)
    # print("psi =", psi)       
    R = Rz(psi) * Ry(theta) * Rx(phi)
    #print(np.round(R, decimals=2))
    return R

def angle(v1, v2):
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle

def _generate_cn_symmetry_helix(outdir, target, name = '', metal_sel ='name NI', n = 4, x_rotations = [0], x2_rotations = [[0, 0, 0, 0]], z_rotation = None):
    '''
    To form a cn symetry helix, One can transform the metal of the target to [0, 0, 0] and contacting atom to [contact_dist, 0, 0] (align on x xias) to simplify rotation.
    Example, in 4 helix bundles, the Chain A B C D, it will form a squre binding geometry.
    x_rotations = [0]. the target can rorate based on x itself.
    x2_rotations = [[0, 2, 0, -2]]. Example in 4 helix bundles, A C should have the same x2_rotaiton, B D should have opposite x2_rotation. So it still keep the squre binding geometry.

    '''
    if name == '':
        name = target.GetTitle()

    metal = target.select(metal_sel)[0]

    nearest_atom = target.select('protein and within 2.83 of index ' + str(metal.getIndex()))[0]

    dist = pr.calcDistance(metal, nearest_atom)

    xyzs = np.array([metal.getCoords(), nearest_atom.getCoords()])

    xyzs_tox = np.array([[0,0,0], [dist, 0, 0]])
    
    target_t1 = target.copy()

    tf1 = pr.calcTransformation(xyzs, xyzs_tox).apply(target_t1) 

    #TO DO: decide the helix direction in y axis metal-contact-helix. 
    #nearest_atom_aa_c = target_t1.select('name C and resindex ' + str(nearest_atom.getResindex()))[0]
    #nearest_atom_aa_o = target_t1.select('name O and resindex ' + str(nearest_atom.getResindex()))[0]

    nearest_atom_aa_c = target_t1.select('name C')
    nearest_atom_aa_o = target_t1.select('name O')

    helix_not_aligned = True

    while helix_not_aligned:
        #v1 = nearest_atom_aa_c.getCoords() - nearest_atom_aa_o.getCoords()
        v1 = np.sum(nearest_atom_aa_c.getCoords() - nearest_atom_aa_o.getCoords(), axis=0)/nearest_atom_aa_c.getCoords().shape[0]
        ag = angle([0, v1[1], v1[2]], [0, 0, 1]) # Projection of 3D Vector onto 2D Plane xyz -> yz
        #print(180* ag/np.pi)
        if (180* ag/np.pi) <= 6:
            helix_not_aligned = False
        else:
            R = getRotation(np.pi*(3/180), 0, 0)
            pr.Transformation(R, np.zeros(3)).apply(target_t1)

    for ix in range(len(x_rotations)):
        target_t2 = target_t1.copy()
        x = x_rotations[ix]

        if x != 0:
            R = getRotation(m.pi*(x/180), 0, 0)
            pr.Transformation(R, np.zeros(3)).apply(target_t2)

        for ix2 in range(len(x2_rotations)):
            
            ags = []
            for i in range(0, n):  
                target_t3 = target_t2.copy()
                ags.append(target_t3)

                R = getRotation(0, 0, m.pi/(n/2)*i)  
                pr.Transformation(R, np.zeros(3)).apply(ags[i])  

                x2 = x2_rotations[ix2][i]                  
                R = getRotation(m.pi*(x2/180), 0, 0)
                pr.Transformation(R, np.zeros(3)).apply(ags[i]) 

            #rotate around z again for special purpose, check 'metalprot/scripts/run_gss_c2.py'.    
            if z_rotation:
                for i in range(0, n):  
                    R = getRotation(0, 0, m.pi*(z_rotation/180))  
                    pr.Transformation(R, np.zeros(3)).apply(ags[i])  

            #print(ags)
            MergeAtomGroup(outdir, ags, string.ascii_uppercase[0:n], name + '_x1_' + str(ix) + '_x2_' + str(ix2) + '_t' + '.pdb')


def MergeAtomGroup(outdir, atomGroups, chains, name, keep_metal = False):
    '''
    The function assume each atomGroup has only one chain. 
    '''
    all_ags = atomGroups[0].select('protein').toAtomGroup()
    all_ags.setChids(chains[0])
    
    for i in range(1, len(atomGroups)):
        ag = atomGroups[i].select('protein').toAtomGroup()
        ag.setChids(chains[i])
        all_ags = all_ags + ag

    if keep_metal:
        all_ags.add(atomGroups[0].select('ion or name NI MN ZN CO CU MG FE').toAtomGroup())

    pr.writePDB(outdir + name + '.pdb', all_ags)


def MergeAtomGroupAuto(outdir, atomGroups, name, keep_metal = False):
    chain_num = 0
    all_ags = atomGroups[0].select('protein and chindex 0').toAtomGroup()
    all_ags.setChids(string.ascii_uppercase[chain_num])
    for i in range(len(atomGroups)):
        ag = atomGroups[i]
        for chindi in np.unique(ag.getChindices()):
            if i == 0 and chindi == 0:
                continue
            chain_num += 1
            chain = ag.select('protein and chindex ' + str(chindi)).toAtomGroup()
            chain.setChids(string.ascii_uppercase[chain_num])
            all_ags = all_ags + chain

    if keep_metal:
        all_ags.add(atomGroups[0].select('ion or name NI MN ZN CO CU MG FE').toAtomGroup())
    
    pr.writePDB(outdir + name + '.pdb', all_ags)
    


def generate_c4s(outdir, target_path, metal_sel, pre_flag = 'target_c4_'):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    target = pr.parsePDB(target_path)

    x_rotations = []
    x2_rotations = []
    for i in range(-10, 11):
        x_rotations.append(i*3)

    for j in range(-5, 6):
        x2_rotations.append([0, j*3, 0, j*3])

    _generate_cn_symmetry_helix(outdir, target, name = pre_flag, metal_sel = metal_sel, n = 4, x_rotations = x_rotations, x2_rotations = x2_rotations )


def cal_phipsi(pr_pdb):
    seq = []
    phi_180 = []
    psi_180 = []
    for p in pr_pdb.iterResidues():
        seq.append(p.getResname())
        try:
            phi_180.append(pr.calcPhi(p))
        except:
            phi_180.append(None)
        try:
            psi_180.append(pr.calcPsi(p))
        except:
            psi_180.append(None)
    return phi_180, psi_180, seq


def cal_dssp(phi, psi, seq):
    character = []
    for i in range(len(seq)):
        s = seq[i]
        if phi[i] and psi[i]:  
            phi_ind = int((phi[i]+ 180)/10 -1)    
            psi_ind = int(35 - (psi[i]+ 180)/10)  
            if s not in constant.APBEL_DICT.keys():
                s = 'ALA'     
            character.append(constant.APBEL_DICT[s][psi_ind, phi_ind])
        else:
            character.append(None)
    return character
    
                

            