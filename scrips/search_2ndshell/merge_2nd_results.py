import os 
import pandas as pd                                                                                                                                               

workdir = '/mnt/e/DesignData/ligands/ligandBB/huong/output_selfcenter_20220126-230900/'                                                                           

outfile = workdir + '_2nd_new.tsv'                                                                                                                                

df_all = pd.DataFrame()

for dir in os.listdir(workdir):
    #print(dir)
    if not os.path.isdir(workdir + dir):
        continue
    

    for dir2 in os.listdir(workdir + dir):

        dir2_path = workdir + dir + '/' + dir2 + '/'
        #print(v)

        if not os.path.isdir(dir2_path):
            continue

        for file in os.listdir(dir2_path):

            if ('_2ndshell_summary' not in file) or ('.tsv' not in file):
                continue
            file_path = dir2_path + file

            df = pd.read_csv(file_path, sep='\t', lineterminator='\n', index_col=False)

            if df.shape[0] < 1:
                continue

            df_all = df_all.append(df.iloc[0])

df_all.to_csv(outfile, sep='\t')
