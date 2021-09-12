from .search import Graph, CombInfo, Search_vdM

class Search_selfcenter(Search_vdM):
    '''
    Inheritated from Search_vdM. 
    Except here is searching the selfcenter vdM database. 
    '''

    def run_neighbor_search(self):
        '''
        All functions need to run the neighbor search.
        '''
        print('run_neighbor_search')

        #TO DO: where should I apply filters: win filter, query_metal filter, phipsi, etc.
        self.neighbor_generate_query_dict()

        self.neighbor_generate_pair_dict()

        if self.parallel:
            self.neighbor_search_wins_pool()
        # else:
        #     self.neighbor_search_wins()

        # self.neighbor_write_represents()

        # self.neighbor_write_summary()

        return


    def neighbor_search_wins_pool(self):


        return 