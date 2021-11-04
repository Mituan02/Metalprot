import itertools


class Graph:
    def __init__(self, wins, all_pos_len):
        #print('graph.')
        self.wins = wins
        self.win_inds = list(range(len(wins)))
        self.all_pos_len = all_pos_len

        self.pair_dict = {}
        for c, c2 in itertools.combinations(self.win_inds, 2):
            self.pair_dict[(c, c2)] = {}

        self.all_paths = []


    def calc_pair_connectivity(self, neighbor_pair_dict):

        for c, c2 in itertools.combinations(self.win_inds, 2):
            wx = self.wins[c]
            wy = self.wins[c2]    
            for r in range(self.all_pos_len):
                connect = neighbor_pair_dict[(wx, wy)][r]
                #print('The len of connect in {} is {}'.format(r, len(connect)))
                if len(connect) <= 0:
                    continue        
                self.pair_dict[(c, c2)][r] = set(connect)
        return

    
    def get_paths(self):

        c = 0
        for r in range(self.all_pos_len):
            self.get_path_helper(c, r, [r])
        return


    def get_path_helper(self, c, r, temp):
        '''
        Dynamic programming. 
        '''
        if len(temp) == len(self.wins):
            self.all_paths.append([t for t in temp])
            return
        #print('c ' + str(c) + ' r ' + str(r))

        if r not in self.pair_dict[(c, c+1)].keys():
            return

        rs = self.pair_dict[(c, c+1)][r]

        for _r in rs:
            #print('col ' + str(c) + ' temp ' + '_'.join([str(t) for t in temp]))
            satisfy = True
            for ci in range(len(temp)-1):
                _r_ci = temp[ci]
                if _r_ci not in self.pair_dict[(ci, c+1)].keys():
                    satisfy = False
                else:
                    if _r not in self.pair_dict[(ci, c+1)][_r_ci]:
                        satisfy = False
                    #break
            if satisfy:
                _temp = temp.copy()
                _temp.append(_r)
                self.get_path_helper(c+1, _r, [t for t in _temp])

        return