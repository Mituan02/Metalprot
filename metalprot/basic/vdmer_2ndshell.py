import prody as pr
import numpy as np


metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

def _generate_2ndshellVdm(query, _2ndshell, _1stshell, metal):
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []

    neary_aas_coords = []
    neary_aas_coords.extend([x for x in _2ndshell.getCoords()])
    neary_aas_coords.extend([x for x in _1stshell.getCoords()])
    neary_aas_coords.append(metal.getCoords())
    coords = np.array(neary_aas_coords)

    chids.extend(['A' for i in range(len(coords))])

    names.extend(_2ndshell.getNames())
    names.extend(_1stshell.getNames())
    names.append(metal.getName())

    resnames.extend(_2ndshell.getResnames())
    resnames.extend(_1stshell.getResnames()) #In a weird case '1996_1dvf_ZN_1', the 1st and 2nd resnames are the same.
    resnames.append(metal.getResname())

    resnums.extend(_2ndshell.getResnums())
    resnums.extend(_1stshell.getResnums())
    resnums.append(metal.getResnum())

    ag = pr.AtomGroup(query.getTitle())
    ag.setCoords(coords)
    ag.setNames(names)
    ag.setResnums(resnums)
    ag.setResnames(resnames)
    ag.setChids(chids)
    return ag

def organize_2ndshellVdm(query):
    '''
    A query for 2nd shell can have different atom orders as an prody.atomGroup.
    The function is supposed to get the atomGroup of the query with the right atom order.
    '''

    metal = query.select(metal_sel)[0]
    metal_resind = metal.getResindex()

    contact_aa = query.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal_resind))
    _1stshell = query.select('protein and resindex ' + ' '.join([str(x) for x in contact_aa.getResindices()]))
    _1stshell_inds = [_1stshell.getResindices()[0]]
    all_resinds = query.select('protein').getResindices()
    _2ndshell_resinds = [x for x in all_resinds if x not in _1stshell_inds]
    #In certain situation, two contact aas are 2ndshell to each other will be ignored.
    if len(_2ndshell_resinds) == 0: 
        print('Not a 2ndshell: ' + query.getTitle())
        return None
    _2ndshell = query.select('protein and resindex ' + ' '.join([str(x) for x in _2ndshell_resinds]))

    ag = _generate_2ndshellVdm(query, _2ndshell, _1stshell, metal)

    return ag


def organize_2ndshellVdm_probe(query):
    '''
    A query for 2nd shell can have different atom orders as an prody.atomGroup.
    The function is supposed to get the atomGroup of the query with the right atom order.

    Different from the previous function, this function tried to apply to 2ndshell vdms from probe calculated Hbond.
    '''

    metal = query.select(metal_sel)[0]

    itags = query.getTitle().split('_')
    if itags[5]!= 'hb':
        print('Please check the naming of the 2nd shell vdm.')
        return 
    
    
    _1stshell = query.select('chid ' + itags[6] + ' resnum ' + str(itags[7]))

    if not _1stshell:
        print(query.getTitle() +  '_1stshell is None')

    _2ndshell = query.select('chid ' + itags[11] + ' resnum ' + str(itags[12]))
    
    
    if not _2ndshell:
        print(query.getTitle() + '_2ndshell is None')
    ag = _generate_2ndshellVdm(query, _2ndshell, _1stshell, metal)
   
    return ag



class SecondShellVDM():
    '''
    The 2ndshell metal-binding vdM class that reprsent H-Bond.
    '''
    def __init__(self, query, id=-1, clu_rank=-1, score=0, clu_num=0, max_clu_num=0, clu_total_num=0, ag_2ndshell = None):
        self.query = query # The prody vdM pdb.
        
        self.id = id # Each vdM has a unique id, for index the vdM library. 
        self.clu_rank = clu_rank # The rank of cluster in the same type of vdM.

        # vdM info
        self.score = score
        self.clu_num = clu_num
        self.max_clu_num = max_clu_num
        self.clu_total_num = clu_total_num

        self.ag_2ndshell = ag_2ndshell
        #self.ag_2ndshell = organize_2ndshellVdm(query)

    #     #Calc contact_resind
    #     self.contact_resind = None
    #     self.aa_type = None

    #     self.set_vdm()


    # def set_vdm(self):
    #     cen_inds = np.unique(self.query.getResindices())
    #     self.contact_resind = cen_inds[int(cen_inds.shape[0]/2)-1]
    #     self.aa_type = self.query.select('resindex ' + str(self.contact_resind)).getResnames()[0] # what the contacting amino acid.

    #     return

    def copy(self):            
        ag_2ndshell = None
        if self.ag_2ndshell:
            ag_2ndshell = self.ag_2ndshell.copy()
        _2ndvdm = SecondShellVDM(self.query.copy(), self.id, self.clu_rank, self.score, self.clu_num, self.max_clu_num, self.clu_total_num, ag_2ndshell)
        return _2ndvdm
