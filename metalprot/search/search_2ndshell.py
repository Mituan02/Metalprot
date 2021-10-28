import os
import numpy as np
import prody as pr

from ..basic.vdmer import get_contact_atom
#from ..basic.vdmer import metal_sel
metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

def construct_pseudo_2ndshellVdm(target, vdM, w):
    '''
    Find the resind of the vdm on target. Then extract resinds of atoms within a distance. 
    Followed by extracting the vdm resind and the atoms resind pairs together with the metal. 
    '''
        
    nearby_aas = target.select('protein and not carbon and not hydrogen and within 10 of resindex ' + str(w))
    nearby_aa_resinds = np.unique(nearby_aas.getResindices())    


    ags = []
    count = 0
    for resind in nearby_aa_resinds:
        if vdM.win and resind in vdM.win:
            continue
        neary_aas_coords = []
        neary_aas_coords.extend(target.select('name N C CA O and resindex ' + str(resind)).getCoords())
        neary_aas_coords.extend(vdM.query.select('name N C CA O and resindex ' + str(vdM.contact_resind)).getCoords())
        neary_aas_coords.extend(vdM.query.select(metal_sel).getCoords())
        coords = np.array(neary_aas_coords)

        names = []
        names.extend(target.select('name N C CA O and resindex ' + str(resind)).getNames())
        names.extend(vdM.query.select('name N C CA O and resindex ' + str(vdM.contact_resind)).getNames())
        names.extend(vdM.query.select(metal_sel).getNames())
        
        atom_contact_pdb = pr.AtomGroup('nearby_bb' + str(count))
        atom_contact_pdb.setCoords(coords)
        atom_contact_pdb.setNames(names)
        ags.append(atom_contact_pdb)
        count +=1

    return ags

    


def supperimpose_2ndshell(ag, vdm, _2ndvdm, rmsd_cut):
    '''
    supperimpose query to ag. 
    '''
    #print('supperimpose_2ndshell ' + query_2nd.query.getTitle())

    if not ag:
        #print('ag None')
        return None
    if not _2ndvdm.ag_2ndshell:
        #print('_2ndvdm ' + _2ndvdm.query.getTitle() + ' None')
        return None
    
    #if _2ndvdm.aa_type != vdm.aa_type:
    if get_contact_atom(_2ndvdm.query).getResname() != vdm.aa_type:
        return None

    if len(_2ndvdm.ag_2ndshell) != len(ag):
        #print('_2ndvdm ' + _2ndvdm.query.getTitle() + ' len: ' + str(len(_2ndvdm.ag_2ndshell)))
        #print('ag ' + ' len: ' + str(len(ag)))
        return None

    _2ndvdm = _2ndvdm.copy()
    transform = pr.calcTransformation(_2ndvdm.ag_2ndshell, ag)
    transform.apply(_2ndvdm.ag_2ndshell)
    transform.apply(_2ndvdm.query)
    rmsd = pr.calcRMSD(ag, _2ndvdm.ag_2ndshell)

    if rmsd <= rmsd_cut:
        candidate = _2ndvdm.copy()
        return (candidate, rmsd)
    return None


'''
Inheritated from Search_vdM. 
Search the 2nd shell h-hond. 
'''

def run_search_2ndshell(comb_dict, target, secondshell_vdms, rmsd_2ndshell):
    '''
    
    '''
    print('run search 2nd-shell vdM.')
    if not comb_dict:
        print('No 1st shell metal-binding vdM found. No need to search 2ndshell.')

    for key in comb_dict.keys():
        search_2ndshell(comb_dict, key, target, secondshell_vdms, rmsd_2ndshell)

    return 


def search_2ndshell(comb_dict, key, target, secondshell_vdms, rmsd_2ndshell):
    '''
    
    '''
    for w in key[0]:
        vdm = comb_dict[key].centroid_dict[w]
        ags = construct_pseudo_2ndshellVdm(target, vdm, w)

        #TO DO: Try to use nearest neighbor.
        candidates = []
        for ag in ags:   
            for _2ndvdm in secondshell_vdms:
                candidate = supperimpose_2ndshell(ag, vdm, _2ndvdm, rmsd_2ndshell)
                if candidate:
                    candidates.append(candidate)

        comb_dict[key].secondshell_dict[w] = candidates

    return


def write_2ndshell(workdir, comb_dict):
    '''
    #Could be combined in write_comb_info.
    '''
    for key in comb_dict.keys():
        outdir = workdir + 'win_' + '-'.join(str(k) for k in key[0]) + '/'
        
        tag = 'win_' + '-'.join([str(k) for k in key[0]]) + '_clu_' + '-'.join(k[0] + '-' + str(k[1]) for k in key[1])
        
        outdir += tag + '_2ndS/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for w in key[0]:
            for c in comb_dict[key].secondshell_dict[w]:
                out_file = outdir + tag + '_w_' + str(w) + '_2ndS_' + c[0].query.getTitle()
                pr.writePDB(out_file, c[0].query)

        with open(outdir + tag + '_2ndshell_summary.tsv', 'w') as f:
            f.write('name\twin\trmsd\n')
            for w in key[0]:
                for c in comb_dict[key].secondshell_dict[w]:
                    name = tag + '_w_' + str(w) + '_2ndS_' + c[0].query.getTitle()
                    f.write(name + '\t' + str(w) + '\t' + str(round(c[1], 2)) + '\n')
    
    return 