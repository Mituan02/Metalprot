import os
import sys
import prody as pr
import numpy as np
import math as m
import string


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


def generate_cn_symmetry_helix(workdir, target, name = '', metal_sel ='name NI', n = 4, x_rotations = [0], x2_rotations = [[0, 0, 0, 0]]):
    '''
    To form a cn symetry helix, One can transform the metal of the target to [0, 0, 0] and contacting atom to [contact_dist, 0, 0] (align on x xias) to simplify rotation.
    Example, in 4 helix bundles, the Chain A B C D, it will form a squre binding geometry.
    x_rotations = [0]. the target can rorate based on x itself.
    x2_rotations = [[0, 2, 0, -2]]. Example in 4 helix bundles, A C should have the same x2_rotaiton, B D should have opposite x2_rotation. So it still keep the squre binding geometry.

    '''
    if name == '':
        name = target.GetTitle()

    metal = target.select(metal_sel)[0]

    first_near = target.select('protein and within 2.83 of index ' + str(metal.getIndex()))[0]

    dist = pr.calcDistance(metal, first_near)

    xyzs = np.array([metal.getCoords(), first_near.getCoords()])

    xyzs_tox = np.array([[0,0,0], [dist, 0, 0]])
    
    target_t1 = target.copy()

    tf1 = pr.calcTransformation(xyzs, xyzs_tox).apply(target_t1)


    for ix in range(len(x_rotations)):
        target_t2 = target_t1.copy()
        x = x_rotations[ix]

        if x != 0:
            R = getRotation(m.pi*2*(x/360), 0, 0)
            pr.Transformation(R, zeros(3)).apply(target_t2)

        for ix2 in range(len(x2_rotations)):
            
            ags = []
            for i in range(0, n):  
                target_t3 = target_t2.copy()
                ags.append(target_t3)

                R = getRotation(0, m.pi/(n/2)*i, 0)
                pr.Transformation(R, zeros(3)).apply(ags[i])   

                x2 = x2_rotations[ix2][i]                  
                R = getRotation(m.pi*2*(x2/360), 0, 0)
                pr.Transformation(R, zeros(3)).apply(ags[i])

            #print(ags)
            MergeAtomGroup(workdir, ags, string.ascii_uppercase[0:n], name + '_x1_' + str(ix) + '_x2_' + str(ix2) + '_t' + '.pdb')




def MergeAtomGroup(workdir, atomGroups, chains, name, keep_metal = False):
    all_ags = atomGroups[0].select('protein').toAtomGroup()

    
    for i in range(1, len(atomGroups)):
        ag = atomGroups[i].select('protein').toAtomGroup()
        ag.setChids(chains[i])
        all_ags = all_ags + ag

    if keep_metal:
        all_ags.add(atomGroups[0].select('ion or name NI MN ZN CO CU MG FE').toAtomGroup())

    pr.writePDB(workdir + name + '.pdb', all_ags)




    
                

            