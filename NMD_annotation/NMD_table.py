#!/usr/bin/python3

import sys
import pandas as pd
import re
from tqdm import tqdm
from operator import add

pos1 = int(sys.argv[1])
pos2 = int(sys.argv[2])

nmdA_annot_df = pd.read_table("/workspace/users/msanchezg/notebooks/NMDetective/nmdA.annot.enst.tsv.gz", sep='\t')
genes_ns_fs_df = pd.read_table("/workspace/users/msanchezg/notebooks/genes_ns_fs.tsv", sep='\t')
genes_ns_fs_df = genes_ns_fs_df.rename(columns={'Feature':'Transcript_ens'})
nmdA_annot_df = pd.merge(nmdA_annot_df,genes_ns_fs_df,how='inner')
nmdA_annot_df['exon_range'] = nmdA_annot_df['End'] - nmdA_annot_df['Start']
c_df = pd.read_table("/workspace/projects/cptac_analysis/data/ensembl_canonical_transcripts.tsv", header = None,sep='\t')
c_df = c_df.rename(columns={0:'Gene_ens', 1:'Transcript_ens',2:'Symbol'})
nmdA_annot_c_df = pd.merge(nmdA_annot_df,c_df,how='inner')

#fragment the whole table
df = nmdA_annot_c_df[pos1:pos2]

#convert nmd_scores column to lists
frames = df.nmd_score.str.split(",",expand=False)
frames_list = frames.tolist()

#convert start positions to list
start = df['Start'].tolist()

#convert exon_range to list
exon_range = df['exon_range'].tolist()

#convert chrom column to lists
chrom = df['CHROM'].tolist()

#convert enst column to lists
enst = df['Transcript_ens'].tolist()

nmd_list = [[],[],[],[]]
for i in tqdm(range(len(df))):

    #create a list for each start position
    start_pos = [start[i]]*(exon_range[i]+1)

    #create a list for each chrom
    chrom_pos = [chrom[i]]*(exon_range[i]+1)

    #create a list for each chrom
    enst_pos = [enst[i]]*(exon_range[i]+1)

    #create a list with all lists
    loc1 = list(range(0,(exon_range[i]+1)))
    location = list(map(add,loc1, start_pos))
    nmd = [chrom_pos,location,enst_pos,frames_list[i]]
    
    nmd_list = list(map(lambda x, y:x + y, nmd_list, nmd))

#convert list to dataframe
nmd_df = pd.DataFrame(nmd_list)
nmd_df = nmd_df.transpose()
nmd_df.columns = ['chrom','location','Feature','nmd_score']

tsv_path = '/workspace/users/msanchezg/notebooks/NMDetective/nmd.scores.ccle.'+str(pos2)+'.tsv'
nmd_df.to_csv(tsv_path, header=True, index=None, sep='\t')
