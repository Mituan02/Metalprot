from math import comb, log
import os
from typing import Dict
from numba.cuda import target
import numpy as np
import itertools
from numpy.lib.function_base import extract
import prody as pr
from prody.measure.transform import superpose
from prody.proteins.pdbfile import writePDB
import scipy.spatial
from scipy.spatial.distance import cdist, dice
import datetime
import math

from sklearn import neighbors

from ..basic import hull
from ..basic import utils
from ..database import core
from ..database import database_extract
from ..basic.vdmer import metal_sel

from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

from .search import Search_vdM, supperimpose_target_bb
from .graph import Graph
from .comb_info import CombInfo

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
        neary_aas_coords.extend(vdM.query.select('bb or ion or name NI MN ZN CO CU MG FE').getCoords())
        coords = np.array(neary_aas_coords)

        names = []
        names.extend(target.select('name N C CA O and resindex ' + str(resind)).getNames())
        names.extend(vdM.query.select('bb or ion or name NI MN ZN CO CU MG FE').getNames())
        
        atom_contact_pdb = pr.AtomGroup('nearby_bb' + str(count))
        atom_contact_pdb.setCoords(coords)
        atom_contact_pdb.setNames(names)
        ags.append(atom_contact_pdb)
        count +=1

    return ags


def supperimpose_2ndshell(ag, _2ndvdm, rmsd_cut):
    '''
    supperimpose query to ag. 
    '''
    #print('supperimpose_2ndshell ' + query_2nd.query.getTitle())
    _2ndvdm = _2ndvdm.copy()
    transform = pr.calcTransformation(_2ndvdm.ag_2ndshell, ag)
    transform.apply(_2ndvdm.ag_2ndshell)
    transform.apply(_2ndvdm.query)
    rmsd = pr.calcRMSD(ag, _2ndvdm.ag)

    if rmsd <= rmsd_cut:
        candidate = _2ndvdm.copy()
        return candidate
    return None

class Search_2ndshell(Search_vdM):
    '''
    Inheritated from Search_vdM. 
    Search the 2nd shell h-hond. 
    '''

    def run_search_2ndshell(self, comb_dict):
        '''
        
        '''
        print('run search 2nd-shell vdM.')
        if not comb_dict:
            print('No 1st shell metal-binding vdM found. No need to search 2ndshell.')

        for key in comb_dict.keys():
            self.search_2ndshell(key)

        return 


    def search_2ndshell(self, key):
        '''
        
        '''
        for w in key[0]:
            vdm = self.best_aa_comb_dict[key].centroid_dict[w]
            ags = construct_pseudo_2ndshellVdm(self.target, vdm, w)

            #TO DO: Try to use nearest neighbor.
            for ag in ags:
                candidates = []
                for _2ndvdm in self.secondshell_vdms:
                    candidate = supperimpose_2ndshell(ag, _2ndvdm, self.rmsd_2ndshell)
                    if candidate:
                        candidates.append(candidate)

            self.best_aa_comb_dict[key].second_dict[w] = candidates

        return


    def write_2ndshell(self, rank, outpath = '/mem_combs/'):
        '''
        #Could be combined in write_comb_info.
        '''
        outdir = self.workdir + outpath
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        comb = self.combs[rank]
        count = 1
        for query in comb.querys:
            count2 = 1
            for _2ndshell in query._2nd_shells:
                pdb_path = outdir + str(rank + 1) + '_' + str(count) + '_2ndshell_' + str(count2) + '_' + str(round(_2ndshell.score, 2)) + '_' + _2ndshell.query.getTitle()
                pr.writePDB(pdb_path, _2ndshell.query)
                count2+=1
            count+=1

        return 