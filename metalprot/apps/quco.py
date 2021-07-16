import prody as pr
import itertools

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 


class Query:
    def __init__(self, query, score = 0, clu_num = 0, clu_total_num = 0, is_bivalent = False, win = None, path = None, ag = None):
        self.query = query
        self.score = score
        self.clu_num = clu_num
        self.clu_total_num = clu_total_num
        self.is_bivalent = is_bivalent

        #Extra properties for special usage.
        self.win = win
        self.path = path
        self.ag = ag
        self._2nd_shells = []
        if self.win is not None and '_win_' not in self.query.getTitle():
            self.query.setTitle(self.query.getTitle().split('.pdb')[0] + '_win_' + '_'.join([str(w) for w in self.win]) + '.pdb')
        
    def set_win(self, win):
        self.win = win
        if self.win is not None and '_win_' not in self.query.getTitle():
            self.query.setTitle(self.query.getTitle().split('.pdb')[0] + '_win_' + '_'.join([str(w) for w in self.win]) + '.pdb')
        

    def to_tab_string(self):
        query_info = self.query.getTitle() + '\t' + str(round(self.score, 2)) + '\t' + str(self.clu_num)  + '\t'+ str(self.clu_total_num)
        return query_info

    def copy(self):
        return Query(self.query.copy(), self.score, self.clu_num, self.clu_total_num, self.is_bivalent, self.win, self.path, self.ag)
    
    def win_str(self):
        return '-'.join([str(w) for w in self.win])

    def write(self, outpath):
        pr.writePDB(outpath, self.query)


class Comb:
    def __init__(self, querys, min_contact_query = None, min_contact_rmsd = None):
        self.querys = querys       
        self.total_score = sum([q.score for q in querys])
        self.total_clu_number= sum([q.clu_num for q in querys])        
        self.scores = [q.score for q in querys]
        self.clu_nums= [q.clu_num for q in querys]
        self.min_contact_query = min_contact_query
        self.min_contact_rmsd = min_contact_rmsd

        self.pair_dists = None
        self.pair_angles = None

    def to_tab_string(self):
        query_names = '||'.join([q.query.getTitle() for q in self.querys])
        query_scores = '||'.join([str(round(n, 2)) for n in self.scores])
        query_clu_nums = '||'.join([str(s) for s in self.clu_nums])
        wins = '||'.join([q.win_str() for q in self.querys])
        vdm_info = str(round(self.total_score, 2)) + '\t' + str(self.total_clu_number) + '\t' + query_names + '\t' + query_scores + '\t' + query_clu_nums + '\t' + wins
        if self.min_contact_query:
            vdm_info += '\t' + self.min_contact_query.query.getTitle() + '\t' + str(round(self.min_contact_query.score, 2)) + '\t' + str(round(self.min_contact_rmsd, 2))
            
        if self.pair_dists:
            vdm_info += '\t' + '||'.join([str(round(d, 2)) for d in self.pair_dists]) + '\t' + '||'.join([str(round(a, 2)) for a in self.pair_angles])
        return vdm_info


    def calc_pair_geometry(self):
        '''
        Calc paired query angle and distance. 
        '''
        #Get metal, binding atom for each binding atom
        mts = []
        cts = []
        for q in self.querys:
            metal = q.query.select(metal_sel)[0]
            _contact_aas = q.query.select('protein and not carbon and not hydrogen and within 2.83 of resindex ' + str(metal.getResindex()))
            #For each aa, only select one contact atom. 
            resindices = _contact_aas.getResindices()
            for rid in resindices:
                #TO DO: Not the first one, but the closest one for the  such ASP contains two contact atoms.
                ct = _contact_aas.select('resindex ' + str(rid))[0]
                mts.append(metal)
                cts.append(ct)

        ct_len = len(cts)

        dist_pair = []
        angle_pair = []
        for i, j in itertools.combinations(range(ct_len), 2):   
            dist = pr.calcDistance(cts[i], cts[j])
            dist_pair.append(dist)
            angle = pr.calcAngle(cts[i], mts[i], cts[j])
            angle2 = pr.calcAngle(cts[i], mts[j], cts[j])
            angle_pair.append((angle + angle2)/2)
        self.pair_dists = dist_pair
        self.pair_angles = angle_pair   
