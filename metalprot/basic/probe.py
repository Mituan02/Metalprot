import os
from collections import defaultdict
import pandas as pd
import prody as pr


#dtypes of dataframe columns in :func:`parse_probe`
probe_dtype_dict = {'interaction': 'category',
                    'chain1': 'category',
                    'resnum1': int,
                    'resname1': 'category',
                    'name1': 'category',
                    'atomtype1': 'category',
                    'chain2': 'category',
                    'resnum2': int,
                    'resname2': 'category',
                    'name2': 'category',
                    'atomtype2': 'category'}


# dataframe columns in :func:`parse_probe`
probe_col_names = list(probe_dtype_dict.keys())

def _parse_probe_line(line):
    """ """
        
    spl = line.split(':')[1:]
    INTERACTION = spl[1]
    CHAIN1 = spl[2][:2].strip()
    RESNUM1 = int(spl[2][2:6])
    RESNAME1 = spl[2][6:10].strip()
    NAME1 = spl[2][10:15].strip()
    ATOMTYPE1 = spl[12]
    CHAIN2 = spl[3][:2].strip()
    RESNUM2 = int(spl[3][2:6])
    RESNAME2 = spl[3][6:10].strip()
    NAME2 = spl[3][10:15].strip()
    ATOMTYPE2 = spl[13]
    if RESNAME1 == 'HOH':
        NAME1 = 'O'
    if RESNAME2 == 'HOH':
        NAME2 = 'O'
    return (INTERACTION, CHAIN1,
            RESNUM1, RESNAME1, NAME1, ATOMTYPE1,
            CHAIN2, RESNUM2, RESNAME2,
            NAME2, ATOMTYPE2)

def probe2csv(outdir, probe_file, ignore_bo=True, include_wc = False, write2file = False):

    probe_data = []
    #with os.popen(cmd) as probefile:
    with open(outdir + probe_file, 'r') as f:
        lines = f.readlines()
        #for line in probefile:
        for line in lines:
            line_data = _parse_probe_line(line)
            probe_data.append(line_data)
    probe_df = pd.DataFrame(probe_data, columns=probe_col_names)
    if write2file:
        probe_df.to_csv(outdir + probe_file + '_probe.csv')
    return probe_df


def extract_2ndshell_probe(core, df):
    '''
    extract 2ndshell from probe dataframe
    '''
    
    df_hbs = pd.DataFrame()
    for resind in core.contact_aa_resinds:
        chid = core.full_pdb.select('resindex ' + str(resind)).getChids()[0]
        resnum = core.full_pdb.select('resindex ' + str(resind)).getResnums()[0]
        # (df['atomtype1'] != 'C') & (df['name1'] != 'O') & (df['name1'] != 'N')
        df_hb = df[(df['interaction'] == 'hb') & (df['chain1'] == chid) & (df['resnum1'] == resnum) & (df['name1'] != 'H') & (df['name1'] != 'O') & (df['name1'] != 'N') & (df['atomtype1'] != 'C') & (df['atomtype2'] != 'C') & (df['atomtype1'] != 'Car') & (df['atomtype2'] != 'Car')]
        df_hb_filter = df_hb[(df_hb['chain1'] != df_hb['chain2']) | (df_hb['resnum1'] != df_hb['resnum2'])]
        df_hbs = pd.concat([df_hbs, df_hb_filter])
        #df_hbs = pd.concat([df_hbs, df_hb.drop_duplicates(subset=['chain1', 'resnum1', 'chain2', 'resnum2'])])

    # # The interaction type 'so' does not like hbond.
    # df_so = df[(df['interaction'] == 'so') & (df['chain1'] == chid) & (df['resnum1'] == resnum) & (df['name1'] != 'O') & (df['name1'] != 'N') & (df['atomtype1']!= 'C') & (df['atomtype2'] != 'C') & (df['atomtype1'] != 'Car') & (df['atomtype2'] != 'Car')]
    # df_hbs = pd.concat([df_hbs, df_so])
    
    return df_hbs