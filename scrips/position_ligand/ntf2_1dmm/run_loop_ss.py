import os
import sys
from rvrsprot.loop import loop_ss
import datetime 


'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/ntf2_1dmm/run_loop_ss.py

'''
class Para():
    workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/search_unknow_result/output__helix6a_10-14-139_HDH.pdb_20220624-222024/'
    outdir = None
    target_file = 'helix6a_10-14-139_tts_0.pdb'


    # loop_topo_sels = {
    #     # Toplogy
    #     'ABC_frt' : {
    #         'A-B': [('A', 21, 3, 'B', 4, 3), ('A', 18, 3, 'B', 7, 3)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
    #         'B-C': [('B', 21, 3, 'C', 6, 3), ('B', 18, 3, 'C', 9, 3)],
    #         'C-D': [('C', 21, 3, 'D', 5, 3), ('C', 18, 3, 'D', 8, 3)],
    #         'D-E': [('D', 22, 3, 'E', 10, 3)],
    #         'E-F': [('E', 21, 3, 'F', 4, 3), ('E', 18, 3, 'F', 7, 3)],
    #         'F-A': [('F', 21, 3, 'A', 4, 3), ('F', 18, 3, 'A', 7, 3)],
    #     },
    #     'ABC_rev' : {
    #         'A-F': [('A', 21, 3, 'F', 5, 3), ('A', 18, 3, 'F', 8, 3)],
    #         'F-E': [('F', 21, 3, 'E', 4, 3), ('F', 18, 3, 'E', 7, 3)],
    #         'E-D': [('E', 21, 3, 'D', 10, 3)],
    #         'D-C': [('D', 21, 3, 'C', 6, 3), ('D', 18, 3, 'C', 9, 3)],
    #         'C-B': [('C', 21, 3, 'B', 6, 3), ('C', 18, 3, 'B', 9, 3)],
    #         'B-A': [('B', 21, 3, 'A', 4, 3), ('B', 18, 3, 'A', 7, 3)],
    #     },
    # }

    loop_topo_sels = {
        # Toplogy
        'ABC_frt' : {
            #'A12': [('A', 20, 26, 'A', 48, 54)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
            #'A23': [('A', 68, 73, 'A', 76, 81)],
            #'A34': [('A', 95, 101, 'A', 122, 128)],
            'A41': [('A', 138, 144, 'A', 6, 12)]
        },
    }

    ### Master search para
    loop_range = [3, 18]
    loop_target_list = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt' # gpu:'/mnt/e/GitHub_Design/master_db/list'
    rmsdCut = 0.5
    master_query_loop_top = 500
    cluster_rmsd = 0.75
    rm_duplicate = True

    ### Output para
    cluster_count_cut = 10
    select_min_rmsd_pdb= False


def main():
    para = Para()
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = para.workdir + 'loop_ss_' + time_tag + '/'
    if para.outdir:
        outdir = para.workdir + para.outdir
    os.makedirs(outdir, exist_ok=True)
    loop_ss.run_loop_ss(outdir, para.workdir + para.target_file, para.loop_topo_sels, para)

    return

if __name__=='__main__':
    main()