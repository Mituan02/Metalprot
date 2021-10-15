'''
After eval_search the hundreds of benchmark pdbs, using this script to analyze the summaries.
'''

import os
import pandas


#workdir = '/mnt/e/DesignData/ligands/LigandBB/zn_eval_bench_mark/2013_new_dist45_phipsi15/'
#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_Seq_core_date_3contact_B45/eval_selfcenter/'

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_3contact/output_eval/' 

#test = 'output_2013_4bm9_ZN_1'

#file = workdir + test + '/_summary.tsv'

#df = pandas.read_csv(file, sep='\t', lineterminator='\n', index_col=False)

NotExist = []

Failed = []

infos = []

max_infos = []

for _dir in os.listdir(workdir):
    # if '.pdb' in _dir or '.tsv' in _dir:
    #     continue
    file = workdir + _dir + '/_summary.tsv'
    if not os.path.isfile(file):
        print('Not exist: ' + _dir)
        NotExist.append(_dir)
        continue

    df = pandas.read_csv(file, sep='\t', lineterminator='\n', index_col=False)

    if len(df) <= 0:
        print('Fail: ' + _dir)
        Failed.append(_dir)
        continue

    #Extract info
    info = []
    info.append(_dir)
    x = df.index[df['eval_is_origin']== True]
    
    info.append(len(df))
    
    if len(x) > 0:
        info.append(len(df['Wins'][x[0]].split('_')))
        info.append(True)
        info.append(df['TotalVdMScore'][x[0]])
        info.append(df['FracScore'][x[0]])
        info.append(df['MultiScore'][x[0]])
        info.append(df['overlap#'][x[0]])
        info.append(df['overlaps#'][x[0]])
        info.append(df['total_clu#'][x[0]])
        info.append(df['clu_nums'][x[0]])
    else:
        info.append(0)
        info.append(False)
        info.append(0)
        info.append(0)
        info.append(0)
        info.append(0)
        info.append('0||0||0')
        info.append(0)
        info.append('0||0||0')

    # ind = df['TotalVdMScore'].idxmax(axis = 1)
    # info.append(df['TotalVdMScore'][ind])
    # ind = df['FracScore'].idxmax(axis = 1)
    # info.append(df['FracScore'][ind])
    # ind = df['MultiScore'].idxmin(axis = 1)
    # info.append(df['MultiScore'][ind])
    # info.append(df['overlap#'][ind])
    # info.append(df['clu_nums'][x[0]])

    infos.append(info)

    #Extract max info
    maxinfo = []
    y = df['overlap#'].idxmax()
    maxinfo.append(_dir)
    maxinfo.append(len(df))
    if len(x) > 0:
        maxinfo.append(len(df['Wins'][x[0]].split('_')))
        maxinfo.append(x == y)
        maxinfo.append(df['TotalVdMScore'][x[0]])
        maxinfo.append(df['FracScore'][x[0]])
        maxinfo.append(df['MultiScore'][x[0]])
        maxinfo.append(df['overlap#'][x[0]])
        maxinfo.append(df['overlaps#'][x[0]])
        maxinfo.append(df['total_clu#'][x[0]])
        maxinfo.append(df['clu_nums'][x[0]])
    else:
        maxinfo.append(0)
        maxinfo.append(False)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append('0||0||0')
        maxinfo.append(0)
        maxinfo.append('0||0||0')


with open(workdir + '_extract.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tOverlap\tOverlaps\tclu_num\tclu_nums\n')
    for info in infos:
        f.write('\t'.join([str(o) for o in info]) + '\n')

with open(workdir + '_extract_max.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tOverlap\tOverlaps\tclu_num\tclu_nums\n')
    for info in infos:
        f.write('\t'.join([str(o) for o in info]) + '\n')

with open(workdir + '_notexist.tsv', 'w') as f:
    for a in NotExist:
        f.write(a + '\n')

with open(workdir + '_failed.tsv', 'w') as f:
    for a in Failed:
        f.write(a + '\n')


with open(workdir + '_extract_split.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tOverlap\tOverlap1\tOverlap2\tOverlap3\tclu_num\tclu_num1\tclu_num2\tclu_num3\n')
    for info in infos:
        title = info[0] 
        overlap = str(info[7])
        overlaps = sorted([int(x) for x in info[8].split('||')])
        clu_num = str(info[9])
        clu_nums = sorted([int(x) for x in info[10].split('||')])
        f.write(title + '\t' + overlap + '\t' + '\t'.join([str(o) for o in overlaps]) + '\t' + clu_num + '\t' + '\t'.join([str(o) for o in clu_nums]) + '\n')


