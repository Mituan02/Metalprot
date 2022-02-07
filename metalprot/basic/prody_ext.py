'''
To help manipulate pdb with prody functions.
'''

import prody as pr
import string
import numpy as np
from . import constant

def transfer2pdb(points, names = None, resnums = None, resnames = None, title = 'MetalMol'):
    '''
    Speically used for Metal binding geometry pdb. 
    '''
    if not names:
        names = ['NI' for i in range(len(points))]
    if not resnums:
        resnums = [0 for i in range(len(points))]
    chains = ['A' for i in range(len(points))]
    mm = pr.AtomGroup(title)
    mm.setCoords(points)
    mm.setResnums(resnums)
    mm.setNames(names)
    mm.setResnames(names)
    mm.setChids(chains)
    if resnames:
        mm.setResnames(resnames)
    return mm


def write2pymol(points, outdir, filename, names = None):
    '''
    Speically used for Metal binding geometry pdb output.
    '''
    mm = transfer2pdb(points, names)
    pr.writePDB(outdir + filename + '.pdb', mm)


def ordered_sel(mobile, target, mob_sel, target_sel):
    '''
    The prody sel always follow the order of atom index. 
    Here is to follow the selection order for alignment.
    mobile: prody object.
    target, prody object.
    mob_sel: list of atom names.
    target_sel: list of atom names.
    '''
    mobile_sel_coords = []
    for s in mob_sel.split(' '):
        mobile_sel_coords.append(mobile.select('name ' + s).getCoords()[0])
    
    target_sel_coords = []
    for s in target_sel.split(' '):
        target_sel_coords.append(target.select('name ' + s).getCoords()[0])

    return np.array(mobile_sel_coords), np.array(target_sel_coords)


def ordered_sel_transformation(mobile, target, mob_sel, target_sel):
    mobile_sel_coords, target_sel_coords = ordered_sel(mobile, target, mob_sel, target_sel)
    transformation = pr.calcTransformation(mobile_sel_coords, target_sel_coords)
    transformation.apply(mobile)
    rmsd = pr.calcRMSD(mobile_sel_coords, target_sel_coords)
    return rmsd


def ordered_sel_rmsd(mobile, target, mob_sel, target_sel):
    mobile_sel_coords, target_sel_coords = ordered_sel(mobile, target, mob_sel, target_sel)
    rmsd = pr.calcRMSD(mobile_sel_coords, target_sel_coords)
    return rmsd


def combine_vdm_into_ag(vdms, tag, geometry, overlapScore = 0, cluScore = 0, ideal_geometry = None):
    '''
    Merge all vdms into one prody AtomGroup.
    Generally for CombInfo.centroid_dict.
    '''
    ag = pr.AtomGroup(tag)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []
    betas = []
    occu = []
    chain_num = 0
    for v in vdms:
        c = v.query.select('not name NI MN ZN CO CU MG FE' )
        c.setChids(string.ascii_uppercase[chain_num])
        coords.extend(c.getCoords())
        chids.extend(c.getChids())
        names.extend(c.getNames())
        resnames.extend(c.getResnames())
        resnums.extend(c.getResnums())
        betas.extend([overlapScore for x in range(len(c))])
        occu.extend([cluScore for x in range(len(c))])
        chain_num += 1

    if not geometry:
        metal_center = pr.calcCenter([v.select('name NI MN ZN CO CU MG FE')[0].getCoords() for v in vdms])
        geometry = transfer2pdb(metal_center, ['NI'])
    geometry.setChids(string.ascii_uppercase[chain_num])
    _geo = geometry.select('name NI')
    coords.extend(_geo.getCoords())
    chids.extend(_geo.getChids())
    names.extend(_geo.getNames())
    resnames.extend(_geo.getResnames())
    resnums.extend(_geo.getResnums())
    betas.append(overlapScore)
    occu.append(cluScore)
    chain_num += 1

    if ideal_geometry:     
        ideal_geometry.setChids(string.ascii_uppercase[chain_num])
        coords.extend(ideal_geometry.getCoords())
        chids.extend(ideal_geometry.getChids())
        names.extend(ideal_geometry.getNames())
        resnames.extend(ideal_geometry.getResnames())
        resnums.extend(ideal_geometry.getResnums())
        betas.extend([overlapScore for i in range(len(ideal_geometry))])
        occu.extend([cluScore for i in range(len(ideal_geometry))])

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    ag.setBetas(betas)
    ag.setOccupancies(occu)
    return ag


def combine_vdm_target_into_ag(target, resind_vdm_dict, write_geo, geometry, title, aa = 'ALA'):
    '''
    combine vdms into target and mutate all other aa 2 ala or gly for ligand position.
    '''

    ag = pr.AtomGroup(title)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []

    if aa == 'GLY':
        bb_sel = 'name N CA C O H'
    elif aa == 'ALA':
        bb_sel = 'name N CA C O H CB'
    else:
        bb_sel = ''

    for i in target.select('protein and name C').getResindices():
        c = target.select('resindex ' + str(i)  + ' ' + bb_sel)

        if i in resind_vdm_dict.keys():
            v = resind_vdm_dict[i]
            query = v.query.select('resindex ' + str(v.contact_resind))
            query.setChids([c.getChids()[0] for x in range(len(query))])
            query.setResnums([c.getResnums()[0] for x in range(len(query))])
            c = query

        else:
            if aa == 'GLY':
                c.setResnums(['GLY' for x in range(len(ala))])
            elif aa == 'ALA':
                '''
                Add CB for gly in original structure.
                '''
                ala = constant.ideal_ala.copy()
                pr.calcTransformation(ala.select('name N CA C'), c.select('name N CA C')).apply(ala)
                ala.setChids([c.getChids()[0] for x in range(len(ala))])
                ala.setResnums([c.getResnums()[0] for x in range(len(ala))])
                c = ala.select(bb_sel)
            
        coords.extend(c.getCoords())
        chids.extend(c.getChids())
        names.extend(c.getNames())
        resnames.extend(c.getResnames())
        resnums.extend(c.getResnums())

    if write_geo:
        if not geometry:
            metal_center = pr.calcCenter([v.select('name NI MN ZN CO CU MG FE')[0].getCoords() for v in resind_vdm_dict.values()])
            geometry = transfer2pdb(metal_center, ['NI'])
        geometry.setChids(string.ascii_uppercase[len(np.unique(target.select('protein').getChids()))])
        _geo = geometry.select('name NI')
        coords.extend(_geo.getCoords())
        chids.extend(_geo.getChids())
        names.extend(_geo.getNames())
        resnames.extend(_geo.getResnames())
        resnums.extend(_geo.getResnums())

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    return ag


def target_to_all_gly_ala(target, title, win_no_mutation = [], aa = 'ALA'):
    '''
    For the prody object target, mutate all the aa into ala or gly. 
    '''
    ag = pr.AtomGroup(title)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []
    betas = []
    occu = []

    if aa == 'GLY':
        bb_sel = 'name N CA C O H'
    if aa == 'ALA':
        bb_sel = 'name N CA C O H CB'
    
    for i in target.select('protein and name C').getResindices():
        if i in win_no_mutation:
            c = target.select('resindex ' + str(i))
        else:
            c = target.select('resindex ' + str(i)  + ' ' + bb_sel)
            if aa == 'GLY':
                c.setResnames(['GLY' for i in range(len(c))])
            elif aa == 'ALA':
                '''
                Add CB for gly in original structure.
                '''
                ala = constant.ideal_ala.copy()
                pr.calcTransformation(ala.select('name N CA C'), c.select('name N CA C')).apply(ala)
                ala.setChids([c.getChids()[0] for x in range(len(ala))])
                ala.setResnums([c.getResnums()[0] for x in range(len(ala))])
                c = ala.select(bb_sel)
        coords.extend(c.getCoords())
        chids.extend(c.getChids())
        names.extend(c.getNames())
        resnames.extend(c.getResnames())
        resnums.extend(c.getResnums())
        betas.extend([0 for x in range(len(c))])
        occu.extend([0 for x in range(len(c))])

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    ag.setBetas(betas)
    ag.setOccupancies(occu)
    return ag