import os
import pandas as pd

workdir = ''

outfile = ''

df_all = pd.DataFrame()


for dir in os.listdir(workdir):
    
    if not os.path.isdir(dir):
        continue
    for v in os.listdir(workdir + dir):

        file_path = workdir + dir + '/' + v
        
        if ('_summary') not in v or ('.tsv' not in v):
            continue
    
        
        df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)

        if df.shape[0] < 1:
            continue

        df['Title'] = v.split('_')[0:-1]
        df_all = df_all.append(df)


df_all.to_csv(outfile, sep='\t')