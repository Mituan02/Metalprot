'''
The basic functions to search 2ndshell.
'''

import prody as pr


def extract_2ndshell_cgs(outdir, target, win_filters):
    cgs = []
    for w in win_filters:
        cg = target.select('protein and heavy and sc and chid ' + w[0] + ' and resnum ' + str(w[1])).toAtomGroup()
        cg.setTitle('w_' + w[0] + str(w[1]))
        cgs.append(cg)
        pr.writePDB(outdir + cg.getTitle() + '.pdb', cg)
    return cgs


def calc_chidres_around_pose(target, chidres, predefined_win_filters, dist = 18):
    '''
    To find a 2ndshell for His Asp Glu, we only need to search the positions around the selected position.
    chidres: The selected position. Ex: ('A', 5)
    dist: The CA distance between the selected position and the candidate position.
    '''
    resinds = target.select('name CA and within ' + str(dist) + ' of (name CG and chid ' + chidres[0] + ' and resnum ' + str(chidres[1]) + ')').getResindices()

    chidreses = []
    for resind in resinds:
        chid = target.select('name CA and resindex ' + str(resind)).getChids()[0]
        res = target.select('name CA and resindex ' + str(resind)).getResnums()[0]
        if (chid, res) in predefined_win_filters:
            continue
        chidreses.append((chid, res))
    return chidreses
