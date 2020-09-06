import pandas as pd
import re
import numpy as np
from tqdm.notebook import trange, tqdm

def cptac_table():
    ####################
    #Load CPTAC dataset
    ####################
    print('Load CPTAC dataset')
    CPTAC = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'LUAD', 'OV', 'UCEC']
    cptac_df = pd.DataFrame(columns=['gene', 'cpct_aliquot', 'protein_expression', 'aliquot', 'sample',
           'Type', 'median', 'stdev', 'norm_protein_expression', 'log2fpkm',
           'log10fpkm', '#Uploaded_variation', 'Location', 'Allele', 'Feature',
           'protein_mutation', 'Protein_position', 'Amino_acids', 'Consequence',
           'Phenotype', 'Ubiquitinases_Mutated', 'Altered_E3_Ligases',
           'Raw_Residual', 'Stability_Change', 'ABS_Stability_Change'])
    for cancer in tqdm(CPTAC, total=len(CPTAC),desc='Load CPTAC data'):
        df = pd.read_csv("/workspace/projects/cptac_analysis/data/"+cancer+"/dataset_irls.gz", sep='\t')
        Dataset = [cancer] * len(df)
        df['Dataset'] = Dataset
        frames = [cptac_df, df]
        cptac_df = pd.concat(frames)
    cptac_df.drop_duplicates(keep='first',inplace=True) 

    ######################################################################################################
    #Calculate the mean of rna values for the repeated mutations (log2fpkm, log10fpkm, Stability_Change)
    ######################################################################################################
    print('Eliminate duplicated rna measures')
    dupl_samples = ['C3L-00908', 'C3N-00545', 'C3N-01825']

    cptac_df['ID'] = cptac_df['gene'].astype(str)+cptac_df['Feature'].astype(str)+cptac_df['sample'].astype(str)+cptac_df['#Uploaded_variation'].astype(str)+cptac_df['Location'].astype(str)

    cptac_dupl_df = cptac_df[cptac_df['sample'].isin(dupl_samples)]
    cptac_dupl2_df = cptac_dupl_df.groupby(['ID'],as_index=False).mean()
    cptac_dupl3_df = cptac_dupl_df.drop(['log2fpkm','log10fpkm','Raw_Residual','Stability_Change','ABS_Stability_Change'],axis=1)
    cptac_dupl3_df.drop_duplicates(subset='ID',keep='first',inplace=True)
    cptac_dupl4_df = pd.merge(cptac_dupl3_df,cptac_dupl2_df,how='left')
    cptac2_df = cptac_df[~(cptac_df['sample'].isin(dupl_samples))]
    cptac_df = pd.concat([cptac2_df,cptac_dupl4_df],ignore_index=True)

    ################
    #Load AF tables
    ################
    
    print('Load AF tables')
    
    CPTAC = ['BRCA','CCRCC','COAD', 'GBM','HNSC','LUAD','OV','UCEC']
    AF_df = pd.DataFrame(columns=['CHROM','POS','ID','REF', 'ALT','QUAL','FILTER','AF_normal',
           'AF_cancer','DONOR'])
    for cancer in tqdm(CPTAC, total=len(CPTAC),desc='Load AF data'):
        df = pd.read_csv("/workspace/projects/cptac_analysis/data/"+cancer+"/AF_"+cancer+".tsv.gz", sep='\t')
        Dataset = [cancer] * len(df)
        df['Dataset'] = Dataset
        frames = [AF_df, df]
        AF_df = pd.concat(frames)

    AF_max_df = AF_df[['DONOR','Dataset','AF_cancer']].groupby('DONOR', as_index=False).max().sort_values(by='AF_cancer',ascending=True)
    AF_max_df = AF_max_df.rename(columns={'AF_cancer':'AF_max'})

    AF_df = pd.merge(AF_df, AF_max_df, how='left')
    AF_df['AF_cancer2'] = AF_df['AF_cancer'] / AF_df['AF_max']

    #Prepare AF_df table for the merge
    AF_df.drop_duplicates(keep='first',inplace=True) 
    AF_df['Del']=AF_df['REF'].str.len()>1
    AF_df['Ins']=AF_df['ALT'].str.len()>1
    AF_df = AF_df.replace({'Del':{True:1,False:0}})
    AF_df = AF_df.replace({'Ins':{True:-1,False:0}})
    AF_list = AF_df[['REF','ALT','Del','Ins']].values.tolist()
    ref2_list=[]
    alt2_list=[]
    for mutation in tqdm(AF_list,desc='Fix Ref and Alt from AF_df'):
        ref = mutation[0]
        alt = mutation[1]
        Del = mutation[2]
        Ins = mutation[3]
        if Del==1:
            ref_list = list(ref)
            ref2 = ''.join(ref_list[1:len(ref_list)])
            alt2 = '-'
        if Ins == -1:
            ref2 = '-'
            alt_list = list(alt)
            alt2 = ''.join(alt_list[1:len(alt_list)])
        if (Del==0)&(Ins==0):
            ref2 = ref
            alt2 = alt
        ref2_list.append(ref2)
        alt2_list.append(alt2)
    AF_df['REF2'] = ref2_list
    AF_df['ALT2'] = alt2_list
    AF_df['POS'] = AF_df['POS'].astype(int)
    AF_df['POS2'] = AF_df['POS']+AF_df['Ins']+AF_df['Del']
    AF_df['#Uploaded_variation'] = AF_df[['DONOR', 'REF2','ALT2']].agg('__'.join, axis=1)
    AF_df["POS2"] = AF_df["POS2"].astype('str')
    AF_df['Location'] = AF_df[['CHROM','POS2']].agg(':'.join,axis=1)
    AF_df = AF_df.rename(columns={'DONOR':'sample','Location':'location1'})
    AF2_df = AF_df[['sample','#Uploaded_variation','location1', 'AF_cancer', 'AF_max']]
    AF3_df = AF2_df.groupby(['sample','#Uploaded_variation','location1'],as_index=False).mean()

    #Prepare cptac_df table for the merge
    cptac_df['location1'] = cptac_df['Location'].str.split('-',1).str[0]

    #Merge cptac_df and AF2_df
    cptac_df = pd.merge(cptac_df, AF3_df, how='left')

    ####################
    #Load nmd_df table
    ####################
    print('Load NMD-score data and merge')
    nmd_df = pd.read_csv("/workspace/users/msanchezg/notebooks/NMDetective/nmd.scores.cptac.tsv.gz", sep='\t')
    nmd_df["location"] = nmd_df["location"].astype('str')
    nmd_df = nmd_df.rename(columns={'location':'location2'})

    #Prepare cptac_df table for marge with nmd_df table
    cptac_mut_df = cptac_df[['ID','Location']][cptac_df['Phenotype']!='WT']
    cptac_mut_df['chrom'] = cptac_mut_df['Location'].str.split(':').str[0]
    cptac_mut_df['location2'] = cptac_mut_df['Location'].str.split(':').str[1]
    cptac_mut_df['location2'] = cptac_mut_df['location2'].str.split('-').str[0]
    cptac_df = pd.merge(cptac_df,cptac_mut_df,how='left')

    #Merge nmd_df with cptac_df
    cptac_df = pd.merge(cptac_df,nmd_df,how='left')

    ########################################################
    #Calculate fold-change log10fpkm mutants vs wt (median)
    ########################################################
    print('Calculate rna fold-change')
    #1) Calculate the median of log10fpkm for the wt per each gene

    cptac_wt_rna_df = cptac_df[['gene','log10fpkm','Dataset']][cptac_df['Phenotype']=='WT']
    cptac_wt_rna_median_df = cptac_wt_rna_df[['gene','log10fpkm','Dataset']].groupby(['gene','Dataset'],as_index=False).median()
    cptac_wt_rna_median_df = cptac_wt_rna_median_df.rename(columns= {'log10fpkm':'log10fpkm_wt_median'})

    #2) Add this value to cptac_df df:

    cptac_df = pd.merge(cptac_df,cptac_wt_rna_median_df, how='left')

    #3) Create a column with the fold change (log10fpkm-log10fpkm_wt_median)

    cptac_df['Fold_change_rna'] = cptac_df['log10fpkm'] - cptac_df['log10fpkm_wt_median']

    ###########################################################
    #Load cdegron table, merge and create cterm_degron column
    ###########################################################
    print('Load cdegron table and merge')
    cdegron_df = pd.read_csv('/workspace/users/msanchezg/notebooks/cdegron_wtnsfs_annot_cptac.tsv', sep = '\t')
    cptac_df = pd.merge(cptac_df,cdegron_df,how='left')

    cptac_df['ID'] = cptac_df['gene'].astype(str)+cptac_df['Feature'].astype(str)+cptac_df['sample'].astype(str)+cptac_df['#Uploaded_variation'].astype(str)+cptac_df['Location'].astype(str)
    df = cptac_df[['ID','GG', 'RG', 'PG', 'XR', 'RXXG', 'EE', 'RXX', 'VX', 'AX', 'A']][~(cptac_df['GG'].isnull())]
    cdegron_list = ['GG', 'RG', 'PG', 'XR', 'RXXG', 'EE', 'RXX', 'VX', 'AX', 'A']
    for cdegron in cdegron_list:
        df[cdegron] = df[cdegron].astype(int)

    df.set_index('ID',inplace=True)

    df = df[df==1].stack().reset_index().drop(0,1)

    df = df.rename(columns={'level_1':'Cterm_degron'})
    cptac_df = pd.merge(cptac_df,df,how='left')

    ########################
    #Add table with last aa
    ########################
    print('Load last_aa table and merge')
    aa_df = pd.read_csv(r'/workspace/users/msanchezg/notebooks/last_aa_annot_cptac.tsv', sep = '\t')
    cptac_df = pd.merge(cptac_df,aa_df,how='left')

    ########################################################################################
    #Add table with relative protein positions (synonymous, missense, nonsense, frameshift)
    ########################################################################################
    print('Load protein_positions_relative table and merge')
    prots_df = pd.read_csv("/workspace/users/msanchezg/notebooks/pos_rel_cptac.tsv", sep='\t')
    prots_df['Protein_position_relative'] = pd.cut(prots_df.Protein_position2,bins=[0,0.25,0.50,0.75,1],labels=['Q1','Q2','Q3','Q4'])
    prots_df = prots_df.drop_duplicates(keep='first')
    prots_df = prots_df[['ID','Protein_position2','Protein_position_relative']]
    wt_df = cptac_df[['ID']][cptac_df['Phenotype']=='WT']
    wt_df['Protein_position_relative'] = 'WT'
    wt_df = wt_df.drop_duplicates(keep='first')
    prots_wt_mut_df = pd.concat([wt_df,prots_df],ignore_index=True)
    prots_wt_mut_df = prots_wt_mut_df.drop_duplicates(keep='first')
    cptac_df = pd.merge(cptac_df,prots_wt_mut_df,how='left')

    ################
    #Creating cterm degron
    ################
    print('Add creating cterm degron column')
    cptac_df['Creating_cterm_degron'] = ~cptac_df['Cterm_degron'].isnull()
    cptac_df = cptac_df.drop_duplicates(subset=['ID'])

    return cptac_df

def ccle_table():
    ####################
    #Load CCLE dataset
    ####################

    # Upload CCLE table (Cancer Cell Line Encyclopedia)
    print('Load CCLE dataset')
    ccle_df = pd.read_csv("/workspace/projects/cptac_analysis/data/CCLE/dataset_irls.gz", sep='\t')

                #############################################################################################################################
    #Calculate the mean of protein values for the repeated mutations (protein_expression, norm_protein_expression, Raw_Residual, Stability_Change, ABS_Stability_Change)
            #############################################################################################################################
    print('Eliminate duplicated protein measures')
    df1 = ccle_df[ccle_df['Phenotype']=='WT'].groupby(['gene','sample','Phenotype'],as_index=False).mean()
    df2 = ccle_df[ccle_df['Phenotype']!='WT'].groupby(['gene','Feature','sample','#Uploaded_variation','Location','Phenotype'],as_index=False).mean()
    ccle2_df = pd.concat([df1,df2])
    ccle3_df = ccle_df.drop(['protein_expression', 'norm_protein_expression', 'Raw_Residual', 'Stability_Change', 'ABS_Stability_Change'],axis=1)
    ccle4_df = pd.merge(ccle3_df,ccle2_df,how='left')
    ccle_df = ccle4_df.drop_duplicates(keep='first')

    ####################
    #Load nmd_df table
    ####################
    print('Load NMD-scores data and merge')
    nmd_cptac_df = pd.read_csv("/workspace/users/msanchezg/notebooks/NMDetective/nmd.scores.cptac.tsv.gz", sep='\t')
    nmd_ccle_df = pd.read_csv("/workspace/users/msanchezg/notebooks/NMDetective/nmd.scores.ccle.tsv.gz", sep='\t')
    frames = [nmd_cptac_df,nmd_ccle_df]
    nmd_df = pd.concat(frames)
    nmd_df["location"] = nmd_df["location"].astype('str')
    nmd_df = nmd_df.rename(columns={'location':'location2'})

    #Prepare ccle_df table for marge with nmd_df table
    ccle_mut_df = ccle_df[~ccle_df['Feature'].isnull()]
    ccle_mut_df['chrom'] = ccle_mut_df['Location'].str.split(':').str[0]
    ccle_mut_df['location2'] = ccle_mut_df['Location'].str.split(':').str[1]
    ccle_mut_df['location2'] = ccle_mut_df['location2'].str.split('-').str[0]
    ccle_df = pd.merge(ccle_df,ccle_mut_df,how='left')

    #Merge nmd_df with ccle_df
    ccle_df = pd.merge(ccle_df,nmd_df,how='left')

    ###############
    #Filer out CNA
    ###############
    print('Filter out CN high/low')
    ccle_df = ccle_df[(ccle_df['cna']<2)&(ccle_df['cna']>-2)]

            ######################################################################################################################################
    #Calculate fold-change logrsem mutants vs wt (median of subtypes in samples with 10 samples or more, and overall median for sample groups of less than 10 samples)
    ######################################################################################################################################
    print('Calculate fold-change rna')
    #1) Calculate the median of log10fpkm for the wt per each gene

    ccle_df['sample_type'] = ccle_df['sample'].str.split("_",1).str[1]

    ccle_wt_rna_df = ccle_df[['gene','logrsem','sample_type']][ccle_df['Phenotype']=='WT']
    ccle_wt_rna_median_df = ccle_wt_rna_df[['gene','logrsem','sample_type']].groupby(['gene'],as_index=False).median()                                               
    ccle_wt_rna_median_df = ccle_wt_rna_median_df.rename(columns= {'logrsem':'logrsem_wt_median'})
    ccle_wt_rna_median2_df = ccle_wt_rna_df[['gene','logrsem','sample_type']].groupby(['gene','sample_type'],as_index=False).median()
    ccle_wt_rna_median2_df = ccle_wt_rna_median2_df.rename(columns= {'logrsem':'logrsem_wt_median_sample_type'})

    #3) Add this value to cptac_df df:

    ccle_df = pd.merge(ccle_df,ccle_wt_rna_median_df, how='left')
    ccle_df = pd.merge(ccle_df,ccle_wt_rna_median2_df, how='left')

    #4) Create a column with the fold change (log10fpkm-log10fpkm_wt_median)

    sample_types_more10 = ['BREAST', 'CENTRAL_NERVOUS_SYSTEM', 'ENDOMETRIUM',
       'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KIDNEY', 'LARGE_INTESTINE',
       'LIVER', 'LUNG', 'OESOPHAGUS', 'OVARY', 'PANCREAS', 'SKIN', 'STOMACH',
       'UPPER_AERODIGESTIVE_TRACT', 'URINARY_TRACT']

    Fold_change_rna = []
    ccle_list = ccle_df[['sample_type','logrsem','logrsem_wt_median_sample_type','logrsem_wt_median']].values.tolist()
    for item in ccle_list:
        if item[0] in sample_types_more10:
            fc = item[1] - item[2]
            Fold_change_rna.append(fc)
        else:
            fc = item[1] - item[3]
            Fold_change_rna.append(fc)

    ccle_df['Fold_change_rna'] = Fold_change_rna

    ################################################################
    #Import cterm degrons table, merge and add cterm_degron column
    ################################################################
    print('Load cdegron table and merge')
    cdegron_df = pd.read_csv('/workspace/users/msanchezg/notebooks/cdegron_wtnsfs_annot_ccle.tsv', sep = '\t')
    ccle_df = pd.merge(ccle_df,cdegron_df,how='left')

    ccle_df['ID'] = ccle_df['gene'].astype(str)+ccle_df['Feature'].astype(str)+ccle_df['sample'].astype(str)+ccle_df['#Uploaded_variation'].astype(str)+ccle_df['Location'].astype(str)
    df = ccle_df[['ID','GG', 'RG', 'PG', 'XR', 'RXXG', 'EE', 'RXX', 'VX', 'AX', 'A']][~(ccle_df['GG'].isnull())]
    cdegron_list = ['GG', 'RG', 'PG', 'XR', 'RXXG', 'EE', 'RXX', 'VX', 'AX', 'A']
    for cdegron in cdegron_list:
        df[cdegron] = df[cdegron].astype(int)

    df.set_index('ID',inplace=True)

    df = df[df==1].stack().reset_index().drop(0,1)

    df = df.rename(columns={'level_1':'Cterm_degron'})
    ccle_df = pd.merge(ccle_df,df,how='left')

    ################################
    #Import last aa table and merge
    ################################
    print('Load last_aa table and merge')
    aa_df = pd.read_csv('/workspace/users/msanchezg/notebooks/last_aa_annot_ccle.tsv', sep = '\t')
    ccle_df = pd.merge(ccle_df,aa_df,how='left')

    #########################################################################################
    #Add table with relative protein positions (synonymous, missense, nonsense, frameshift)
    #########################################################################################
    print('Load protein_positions_relative table and merge')
    prots_df = pd.read_csv("/workspace/users/msanchezg/notebooks/pos_rel_ccle.tsv", sep='\t')
    prots_df['Protein_position_relative'] = pd.cut(prots_df.Protein_position2,bins=[0,0.25,0.50,0.75,1],labels=['Q1','Q2','Q3','Q4'])
    prots_df = prots_df.drop_duplicates(keep='first')
    prots_df = prots_df[['ID','Protein_position2','Protein_position_relative']]
    wt_df = ccle_df[['ID']][ccle_df['Phenotype']=='WT']
    wt_df['Protein_position_relative'] = 'WT'
    wt_df = wt_df.drop_duplicates(keep='first')
    prots_wt_mut_df = pd.concat([wt_df,prots_df],ignore_index=True)
    prots_wt_mut_df = prots_wt_mut_df.drop_duplicates(keep='first')
    ccle_df = pd.merge(ccle_df,prots_wt_mut_df,how='left')

    #######################
    #Creating cterm degron
    #######################
    print('Add creating cterm degron column')
    ccle_df['Creating_cterm_degron'] = ~ccle_df['Cterm_degron'].isnull()
    ccle_df = ccle_df.drop_duplicates(subset=['ID'])

    return ccle_df

