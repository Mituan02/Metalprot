import itertools


class Graph:
    def __init__(self, wins, all_pos_lens):
        #print('graph.')
        self.wins = wins
        self.win_inds = list(range(len(wins)))
        self.all_pos_lens = all_pos_lens

        self.pair_dict = {}
        for c, c2 in itertools.permutations(self.win_inds, 2):
            self.pair_dict[(c, c2)] = set()

        self.all_paths = []


    def calc_pair_connectivity(self, neighbor_pair_dict):

        for c, c2 in itertools.permutations(self.win_inds, 2):
            wx = self.wins[c]
            wy = self.wins[c2]    
            for r in range(self.all_pos_lens[c]):
                connect = neighbor_pair_dict[(wx, wy)][r]
                #print('The len of connect in {} is {}'.format(r, len(connect)))
                if len(connect) <= 0:
                    continue     
                for y in connect:
                    self.pair_dict[(c, c2)].add((r, y))
        return

    
    def get_paths(self):

        c = 0
        for r in range(self.all_pos_lens[c]):
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
        rs = []
        for i in range(self.all_pos_lens[c+1]):
            if (r, i) in self.pair_dict[(c, c+1)]:
                rs.append(i)

        if len(rs) == 0:
            return      
        for _r in rs:
            #print('col ' + str(c) + ' temp ' + '_'.join([str(t) for t in temp]))
            satisfy = True
            for ci in range(len(temp)-1):
                tv = temp[ci]
                if not (tv, _r) in self.pair_dict[(ci, c+1)]:
                    satisfy = False
                    #break
            if satisfy:
                _temp = temp.copy()
                _temp.append(_r)
                self.get_path_helper(c+1, _r, [t for t in _temp])

        return