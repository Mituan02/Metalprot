import os
import pandas as pd

#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/_Seq_core_3contact_all_dup/_20211202_eval_dup/'
#workdir = '/mnt/e/DesignData/ligands/ZN_rcsb_new/20211119/_Seq_core_3contact_new/_20211202_eval_new/'
workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/family_3vsy_220118_result/'

outfile = workdir + '_summary_new.tsv'

df_all = pd.DataFrame()

no_answer = []
for dir in os.listdir(workdir):
    #print(dir)
    if not os.path.isdir(workdir + dir):
        continue

    exist = False
    for v in os.listdir(workdir + dir):

        file_path = workdir + dir + '/' + v
        #print(v)
        
        if ('_summary' not in v) or ('.tsv' not in v):
            continue
        
        df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)

        if df.shape[0] < 1:
            continue
        exist = True
        df['Title'] = '_'.join(dir.split('_')[0:-1])
        df_all = df_all.append(df.iloc[0])
    if not exist:
        no_answer.append('_'.join(dir.split('_')[0:-1]))

df_all.to_csv(outfile, sep='\t')

with open(workdir + '_summary_new_noanswer.tsv', 'w') as f:
    for nn in no_answer:
        f.write(nn + '\n')

