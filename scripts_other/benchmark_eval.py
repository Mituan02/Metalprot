import os
import pandas

workdir = '/mnt/e/DesignData/ligands/LigandBB/zn_eval_bench_mark/2013_heavy_dist45_phipsi15/'



#test = 'output_2013_4bm9_ZN_1'

#file = workdir + test + '/_summary.tsv'

#df = pandas.read_csv(file, sep='\t', lineterminator='\n', index_col=False)

NotExist = []

Failed = []

infos = []

for _dir in os.listdir(workdir):
    if '.pdb' in _dir or '.tsv' in _dir:
        continue
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
    info = []
    info.append(_dir)
    x = df.index[df['eval_is_origin']== True]
    info.append(len(df))
    
    if len(x) > 0:
        info.append(True)
        info.append(df['TotalVdMScore'][x[0]])
        info.append(df['FracScore'][x[0]])
        info.append(df['MultiScore'][x[0]])
    else:
        info.append(False)
        info.append(0)
        info.append(0)
        info.append(0)


    ind = df['TotalVdMScore'].idxmax(axis = 1)
    info.append(df['TotalVdMScore'][ind])
    ind = df['FracScore'].idxmax(axis = 1)
    info.append(df['FracScore'][ind])
    ind = df['MultiScore'].idxmin(axis = 1)
    info.append(df['MultiScore'][ind])
    infos.append(info)

with open(workdir + '_extract.tsv', 'w') as f:
    for info in infos:
        f.write('\t'.join([str(o) for o in info]) + '\n')

with open(workdir + '_notexist.tsv', 'w') as f:
    for a in NotExist:
        f.write(a + '\n')

with open(workdir + '_failed.tsv', 'w') as f:
    for a in Failed:
        f.write(a + '\n')




    
