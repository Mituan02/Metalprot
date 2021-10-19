'''
After eval_search the hundreds of benchmark pdbs, using this script to analyze the summaries.
'''

import os
import pandas


#workdir = '/mnt/e/DesignData/ligands/LigandBB/zn_eval_bench_mark/2013_new_dist45_phipsi15/'
#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20210624/_Seq_core_date_3contact_B45/eval_selfcenter/'

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_date_3contact/20211015_out_eval2/' 

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
        info.append(df['vdm_scores'][x[0]])
        info.append(df['FracScore'][x[0]])
        info.append(df['MultiScore'][x[0]])
        info.append(df['overlap#'][x[0]])
        info.append(df['overlaps#'][x[0]])
        info.append(df['total_clu#'][x[0]])
        info.append(df['clu_nums'][x[0]])
        info.append(df['max_clu_nums'][x[0]])
    else:
        info.append(0)
        info.append(False)
        info.append(0)
        info.append(0)
        info.append(0)
        info.append(0)
        info.append(0)
        info.append('0||0||0')
        info.append(0)
        info.append('0||0||0')
        info.append('0||0||0')

    infos.append(info)

    #Extract max info
    maxinfo = []
    y = df['overlap#'].idxmax()
    maxinfo.append(_dir)
    maxinfo.append(len(df))
    if len(x) > 0:
        maxinfo.append(len(df['Wins'][y].split('_')))
        maxinfo.append(x[0] == y)
        maxinfo.append(df['TotalVdMScore'][y])
        maxinfo.append(df['vdm_scores'][y])
        maxinfo.append(df['FracScore'][y])
        maxinfo.append(df['MultiScore'][y])
        maxinfo.append(df['overlap#'][y])
        maxinfo.append(df['overlaps#'][y])
        maxinfo.append(df['total_clu#'][y])
        maxinfo.append(df['clu_nums'][y])
        maxinfo.append(df['max_clu_nums'][y])
    else:
        maxinfo.append(0)
        maxinfo.append(False)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append(0)
        maxinfo.append('0||0||0')
        maxinfo.append(0)
        maxinfo.append('0||0||0')
        maxinfo.append('0||0||0')
    max_infos.append(maxinfo)

with open(workdir + '_extract.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tscores\tOriginFracScore\tOriginMultiScore\tOverlap\tOverlaps\tclu_num\tclu_nums\tmax_clu\n')
    for info in infos:
        # if info[8] <= 0:
        #     continue
        f.write('\t'.join([str(o) for o in info]) + '\n')

with open(workdir + '_extract_max.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tscores\tOriginFracScore\tOriginMultiScore\tOverlap\tOverlaps\tclu_num\tclu_nums\tmax_clu\n')
    for info in max_infos:
        f.write('\t'.join([str(o) for o in info]) + '\n')

with open(workdir + '_notexist.tsv', 'w') as f:
    for a in NotExist:
        f.write(a + '\n')

with open(workdir + '_failed.tsv', 'w') as f:
    for a in Failed:
        f.write(a + '\n')


with open(workdir + '_extract_split.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tOverlap\tOverlap1\tOverlap2\tOverlap3\tclu_num\tclu_num1\tclu_num2\tclu_num3\tTotalScore\tscore1\tscore2\tscore3\n')
    for info in infos:
        if info[8] <= 0:
            continue
        title = info[0] 
        overlap = str(info[8])
        overlaps = sorted([int(x) for x in info[9].split('||')])
        f.write(title + '\t' + overlap + '\t' + '\t'.join([str(o) for o in overlaps]) + '\t')
        
        clu_num = str(info[10])
        clu_nums = sorted([int(x) for x in info[11].split('||')])
        f.write(clu_num + '\t' + '\t'.join([str(o) for o in clu_nums]) + '\t')

        score = str(info[4])
        scores = sorted([float(x) for x in info[5].split('||')])
        f.write(score + '\t' + '\t'.join([str(o) for o in scores]) + '\n')
        
with open(workdir + '_extract_split_ratio.tsv', 'w') as f:
    #f.write('Title\tTotalSolutions\twin_size\tIsOriginVdm\tOriginTotalVdMScore\tOriginFracScore\tOriginMultiScore\tMaxTotalVdMScore\tMaxFracScore\tMinMultiScore\n')
    f.write('Title\tclu_num1\tclu_num2\tclu_num3\tmax_clu1\tmax_clu2\tmax_clu3\tratio1\tratio2\tratio3\todratio1\todratio2\todratio3\n')
    for info in infos:
        if info[8] <= 0:
            continue
        title = info[0] 
        f.write(title + '\t')

        clu_nums = [int(x) for x in info[11].split('||')]
        f.write( '\t'.join([str(o) for o in clu_nums]) + '\t')

        max_clus = [int(x) for x in info[12].split('||')]
        f.write('\t'.join([str(o) for o in max_clus]) + '\t')

        ratios = [clu_nums[i]/max_clus[i] for i in range(len(clu_nums))]
        f.write('\t'.join([str(o) for o in ratios]) + '\t')

        ratio_ordered = sorted(ratios)
        f.write('\t'.join([str(o) for o in ratio_ordered]) + '\n')


        