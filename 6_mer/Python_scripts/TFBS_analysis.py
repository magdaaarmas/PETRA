import os
import pickle
import pandas as pd
import numpy as np
import sys
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests

########################################################################################################################

gene=sys.argv[1]
FIMO_output_path=sys.argv[2]
FIMO_output_file=f"{FIMO_output_path}/fimo.tsv"
save_folder_general=sys.argv[3]
save_folder=f"{save_folder_general}TFBS_analysis/"
os.makedirs(save_folder, exist_ok=True)
genes_expressed_in_jurkat_ls=sys.argv[4] #genes_expressed_in_jurkat_ls='/Volumes/lab-findlayg/home/users/armasm/PETRA_large_screen_analysis/jurkat_expression/genes_expressed_in_jurkat_ls.pkl'

max_score_threshold=0.1 #splice AI maximum score for considering for analysis

########################################################################################################################
def WT_analysis(df):
    df=df.copy()
    WT_info=[]
    for gene in genes:
        df_gene = df[df['gene'] == gene].reset_index(drop=True)
        #analyse WT matches
        gene_WT_matches=df_gene[df_gene['ID']=='WT']
        gene_WT_matches_ls=gene_WT_matches['motif_alt_id'].tolist()
        #check if WT matches are interrupted by insertion
        condition=(gene_WT_matches['start']<=20)&(gene_WT_matches['stop']>20)
        gene_WT_matches_interrupted=gene_WT_matches[condition]
        gene_WT_matches_interrupted_ls=gene_WT_matches_interrupted['motif_alt_id'].tolist()
        WT_info.append([gene, gene_WT_matches_ls, gene_WT_matches_interrupted_ls])
    WT_analysis=pd.DataFrame(WT_info, columns=['gene', 'matches_in_WT', 'interrupted_BS'])
    WT_analysis.to_csv(f"{output_path}TFBS_in_WT.csv")
    return WT_analysis
def remove_non_insert(df):
    df=df.copy()
    #remove matches to WT
    df=df[df['ID']!='WT']
    #remove the matches that don't include the insert
    condition = (df['start'] <= 26) & (df['stop'] >= 21)
    df=df[condition]
    return df

def obtain_gene_hits(df):
    df = df.copy()
    genes_hits={}
    for gene in genes:
        genes_hits[gene]={}
        df_gene=df[df['gene']==gene]
        hits=df_gene['motif_alt_id'].unique().tolist()
        for hit in hits:
            sequences_to_hit=df_gene[df_gene['motif_alt_id']==hit]['ID'].tolist()
            genes_hits[gene][hit]=sequences_to_hit
    return genes_hits
def open_pickle(file):
    with open(file, 'rb') as fp:
        df=pickle.load(fp)
    return df
def merge_TFBS_ES(ES_df, dict_of_hits):
    new_dictionary={}
    df=ES_df.copy()
    for hit, sequences_ls in dict_of_hits[gene].items():
        df[hit]=np.where(df['sequence'].isin(sequences_ls), True, False)
    return df
def compute_TFBS_statistics(hits_list_dict, gene, score_df):
    dictionary = hits_list_dict.copy()
    # open data
    data = score_df.copy()
    # select hits
    gene_hits = dictionary[gene]

    # calculate statistics
    hits_statistics_KS = []
    for hit in gene_hits:
        positive = data[data[hit] == True]['filtered_score'].tolist()
        positive_median = np.median(positive)
        negative = data[data[hit] == False]['filtered_score'].tolist()
        negative_median = np.median(negative)
        log2FC = positive_median - negative_median
        if log2FC>0:
            change='TF_increase'
        else:
            change='TF_decrease'

        #statistical testing
        if len(positive) < 1 or len(negative) < 1:
            list = (hit, 'no_items','', '', '')
            list2 = list
        else:
            stat, p = ks_2samp(positive, negative)
            list = (hit, log2FC, change, f"{stat}_KS", p)
            list2 = list
        hits_statistics_KS.append(list2)
        print(f"done with {hit}")

    statistics_KS = pd.DataFrame(hits_statistics_KS, columns=['TF', 'log2FC','change_type', 't', 'p'])
    return statistics_KS


def p_value_corrections(df_with_stats):
    info=df_with_stats
    #filter out errors in statistical processing
    info=info[info['log2FC']!='no_items']
    #apply BH correction
    _, corrected_p_values, _, _ = multipletests(info['p'], method='fdr_bh')
    info['BH_q']=corrected_p_values
    info['BH_sig5']=np.where(info['BH_q']<=0.05, True, False)
    info['BH_sig1'] = np.where(info['BH_q'] <= 0.01, True, False)
    
    return info
########################################################################################################################
#------ merge TF hits to varinats ------ #

#obtain fimo scores
FIMO_scores=pd.read_csv(FIMO_output_file, sep='\t')

#re-order the dataframe for processing
#drop the last 3 rows
FIMO_scores = FIMO_scores.iloc[:-3]
FIMO_scores['gene']=FIMO_scores['sequence_name'].apply(lambda x: str(x).split('_')[0])
genes=list(set(FIMO_scores['gene']))
FIMO_scores['ID']=FIMO_scores['sequence_name'].apply(lambda x: str(x).split('_')[1])
FIMO_scores['TF_name']=FIMO_scores['motif_alt_id'].apply(lambda x: x.replace('(var.2)', ''))
FIMO_scores['TF_name']=FIMO_scores['TF_name'].apply(lambda x: x.split('.')[-1].upper())

#import list of genes expressed in jurkat cells
expressed_genes=open_pickle(genes_expressed_in_jurkat_ls)
expressed_genes=[s.upper() for s in expressed_genes]
FIMO_scores['expressed_in_jurkats']=FIMO_scores['TF_name'].isin(expressed_genes)
FIMO_scores=FIMO_scores[FIMO_scores['expressed_in_jurkats']==True]

#Remove the TFBS that match to regions that DONT include the insert
FIMO_insert_relevant=remove_non_insert(FIMO_scores)
#save the filtered FIMO scores
FIMO_insert_relevant.to_csv(f"{save_folder}{gene}_FIMO_scores_insertion_relevant.csv")
#After this point the information about scores and p-values is lost and only the TFBS that were hits based on the threshold
#set up in FIMO will remain.

#separate by gene and group into insert ID matches
hits_dict_ls=obtain_gene_hits(FIMO_insert_relevant)
with open(f"{save_folder}{gene}_FIMO_dict_hits_list.pkl", 'wb') as fp:
    pickle.dump(hits_dict_ls, fp)

#obtain ES
expression_scores_df=pd.read_csv(f"{save_folder_general}{gene}_scores_variant_features.csv")

#merge ES to TFBS detection
ES_TFBS_merged_df=merge_TFBS_ES(expression_scores_df, hits_dict_ls)

ES_TFBS_merged_df.to_csv(f"{save_folder}{gene}_ES_TFBS_merged.csv", index=False)


#------ STATISTICS ------ #

df=ES_TFBS_merged_df.copy()

#filter out vairants with high spliceAI score
df['filter_out'] = np.where(df['max_score'] >= max_score_threshold,True, False)
df = df[df['filter_out'] == False]
stat_KS_df=compute_TFBS_statistics(hits_dict_ls, gene, df)
stat_KS_df=p_value_corrections(stat_KS_df)

stat_KS_df.to_csv(f"{save_folder}{gene}_TFBS_analysis.csv")

