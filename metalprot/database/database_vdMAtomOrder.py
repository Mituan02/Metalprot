'''
The two oxygens of ASP or GLU vdMs is randomly assigned. 
The function here is meant to order them in consistant way.
'''

import os
import prody as pr
from numpy.core.fromnumeric import argmin

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


def get_contact_info(pdb, sel = 'protein and oxygen'):
    metal = pdb.select(metal_sel)[0]
    _contact_aas = pdb.select(sel + ' and within 2.83 of resindex ' + str(metal.getResindex()))
    if len(_contact_aas) > 1:
        dists = [None]*len(_contact_aas)
        for i in range(len(_contact_aas)):
            dist = pr.calcDistance(metal, _contact_aas[i])
            dists[i] = dist
        contact_aa = _contact_aas[argmin(dists)]
    else:
        contact_aa = _contact_aas[0]

    contact_ind = contact_aa.getIndex()
    contact_name = pdb.select('index ' + str(contact_ind))[0].getName()
    contact_resind = contact_aa.getResindex()
    contact_resname = pdb.select('index ' + str(contact_ind))[0].getResname()

    return contact_ind, contact_name, contact_resind, contact_resname


def asp_glu_oxy_shift(pdb):
    '''
    The Oxygen atom in metal-asp or metal-glu binding vdMs are not ordered in a consistant way.
    The function is to order them consistantly with -O named O1, =O named O2.
    '''

    contact_ind, contact_name, contact_resind, contact_resname = get_contact_info(pdb)

    if contact_resname == 'ASP':
        # For ASP, the sc contacting oxygen is OD1 or OD2. Shift OD2 to OD1 if OD2 is the contacting atom.
        if contact_name == 'OD1':
            print(pdb.getTitle() + ' shifted.')
            contact_coord = pdb.select('index ' + str(contact_ind))[0].getCoords()

            other_o_ind = pdb.select('resindex ' + str(contact_resind) + ' and name OD2')[0].getIndex()
            other_coord = pdb.select('index ' + str(other_o_ind))[0].getCoords()
            
            pdb.select('index ' + str(contact_ind))[0].setCoords(other_coord)
            pdb.select('index ' + str(other_o_ind))[0].setCoords(contact_coord)

    if contact_resname == 'GLU':
        # For GLU, the sc contacting oxygen is OE1 or OE2. Shift OE2 to OE1 if OE2 is the contacting atom.
        if contact_name == 'OE1':
            print(pdb.getTitle() + ' shifted.')
            contact_coord = pdb.select('index ' + str(contact_ind))[0].getCoords()

            other_o_ind = pdb.select('resindex ' + str(contact_resind) + ' and name OE2')[0].getIndex()
            other_coord = pdb.select('index ' + str(other_o_ind))[0].getCoords()
            
            pdb.select('index ' + str(contact_ind))[0].setCoords(other_coord)
            pdb.select('index ' + str(other_o_ind))[0].setCoords(contact_coord)

    return