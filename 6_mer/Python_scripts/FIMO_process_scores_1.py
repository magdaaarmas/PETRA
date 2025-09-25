import os
import pickle
import pandas as pd
import numpy as np

########################################################################################################################

genes = ['IL2RA', 'VAV1', 'CD28', 'OTUD7B']
FIMO_threshold='0001'
path=f"/Volumes/lab-findlayg/home/users/armasm/PETRA_large_screen_analysis/FIMO_{FIMO_threshold}/"
FIMO_output_file=f"{path}/output/fimo.tsv"
expression_scores_path=f"/Volumes/lab-findlayg/home/users/armasm/PETRA_large_screen_analysis/filtered_scores/"
genes_expressed_in_jurkat_ls='/Volumes/lab-findlayg/home/users/armasm/PETRA_large_screen_analysis/jurkat_expression/genes_expressed_in_jurkat_ls.pkl'

output_path=f"{expression_scores_path}/FIMO_analysis_{FIMO_threshold}/"
os.makedirs(output_path, exist_ok=True)

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
def merge_TFBS_ES(ES_dict, dict_of_hits):
    dictionary=ES_dict.copy()
    new_dictionary={}
    for gene in genes:
        df=dictionary[gene]
        for hit, sequences_ls in dict_of_hits[gene].items():
            df[hit]=np.where(df['sequence'].isin(sequences_ls), True, False)
        new_dictionary[gene]=df
    return new_dictionary
########################################################################################################################

#obtain fimo scores
FIMO_scores=pd.read_csv(FIMO_output_file, sep='\t')

#re-order the dataframe for processing
#drop the last 3 rows
FIMO_scores = FIMO_scores.iloc[:-3]
FIMO_scores['gene']=FIMO_scores['sequence_name'].apply(lambda x: str(x).split('_')[0])
FIMO_scores['ID']=FIMO_scores['sequence_name'].apply(lambda x: str(x).split('_')[1])
FIMO_scores['TF_name']=FIMO_scores['motif_alt_id'].apply(lambda x: x.replace('(var.2)', ''))
FIMO_scores['TF_name']=FIMO_scores['TF_name'].apply(lambda x: x.split('.')[-1].upper())

#import list of genes expressed in jurkat cells
expressed_genes=open_pickle(genes_expressed_in_jurkat_ls)
expressed_genes=[s.upper() for s in expressed_genes]
FIMO_scores['expressed_in_jurkats']=FIMO_scores['TF_name'].isin(expressed_genes)
FIMO_scores=FIMO_scores[FIMO_scores['expressed_in_jurkats']==True]

#Analyse presence and interruption of TFBS in WT sequence
WT_df=WT_analysis(FIMO_scores)

#Remove the TFBS that match to regions that DONT include the insert
FIMO_insert_relevant=remove_non_insert(FIMO_scores)
#save the filtered FIMO scores
FIMO_insert_relevant.to_csv(f"{output_path}FIMO_scores_insertion_relevant.csv")
#After this point the information about scores and p-values is lost and only the TFBS that were hits based on the threshold
#set up in FIMO will remain.

#separate by gene and group into insert ID matches
hits_dict_ls=obtain_gene_hits(FIMO_insert_relevant)
with open(f"{output_path}FIMO_dict_hits_list.pkl", 'wb') as fp:
    pickle.dump(hits_dict_ls, fp)

#obtain ES
expression_scores_dict={}
for gene in genes:
    expression_scores_dict[gene]=open_pickle(f"{expression_scores_path}{gene}_filtered_scores_merged_filtered_ATG_SAI.pkl")

#merge ES to TFBS detection
ES_TFBS_merged_dict=merge_TFBS_ES(expression_scores_dict, hits_dict_ls)
#save the merged file in directory with reference to threshold
for gene in genes:
    ES_TFBS_merged_dict[gene].to_csv(f"{output_path}{gene}_ES_TFBS_merged.csv", index=False)
    with open(f"{output_path}{gene}_ES_TFBS_merged.pkl", 'wb') as fp:
        pickle.dump(ES_TFBS_merged_dict[gene], fp)