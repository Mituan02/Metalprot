'''
Here we want to filter the core we already extracted.
Filter is based on B factor, occupancy and number of contacting aa. 
'''
from metalprot.database import database_extract

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/database_generate/database_filter.py
'''

workdir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20210624/_Seq_core_date/'

outdir = '/mnt/e/DesignData/ligands/Zn_rcsb_datesplit/20211013/_Seq_core_filter/'

database_extract.rm_core(workdir, outdir)