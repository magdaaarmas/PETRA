import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import sys
import os
from Bio.Seq import Seq
import pickle
import seaborn as sns
import scipy

########################################################################################################################

gene=sys.argv[1]
spliceAI_folder=sys.argv[2] #spliceai_file=f"{path}{gene}_output_{spliceAI_distance}.vcf"
coding_strand_info_file=sys.argv[3] #path+'coding_strand_info.csv'
save_folder=sys.argv[4]
filtered_scores_path=sys.argv[5]

splice_distance_used='300'
########################################################################################################################
def compute_reverse_complement_df(sequence):
    return str(Seq(sequence).reverse_complement())

########################################################################################################################

#import coding strand info
coding_strand_info=pd.read_csv(coding_strand_info_file)
coding_strand=coding_strand_info.loc[coding_strand_info['gene']==gene, 'strand'].iloc[0]

#import and format the .vcf file
column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
spliceai_file=f"{spliceAI_folder}/{gene}_output_{splice_distance_used}.vcf"
splice_scores = pd.read_csv(spliceai_file, sep='\t', comment='#', header=None, names=column_names)
info_df = splice_scores['INFO'].str.split('|', expand=True)
info_df.columns = ['ALLELE', 'SYMBOL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
splice_scores = pd.concat([splice_scores, info_df], axis=1)
splice_scores['ALLELE'] = splice_scores['ALLELE'].str[10:]
splice_scores.drop(columns=['INFO', 'ID', 'QUAL', 'FILTER'], inplace=True)
splice_scores['DS_AG'] = pd.to_numeric(splice_scores['DS_AG'], errors='coerce')
splice_scores['DS_DG'] = pd.to_numeric(splice_scores['DS_DG'], errors='coerce')
splice_scores['DS_DL'] = pd.to_numeric(splice_scores['DS_DL'], errors='coerce')
splice_scores['DS_AL'] = pd.to_numeric(splice_scores['DS_AL'], errors='coerce')
splice_scores.rename(columns={'ALLELE':'sequence'}, inplace=True)

#modify allele to match scoring direction
if coding_strand=='negative':
    splice_scores['sequence']=splice_scores['sequence'].apply(compute_reverse_complement_df)


list_of_scores=['DS_AG', 'DS_DG', 'DS_DL', 'DS_AL']
splice_scores['max_score']=splice_scores[list_of_scores].max(axis=1)


#incorporate splice scores
#import scores and extract replicate number
scores=pd.read_csv(f"{filtered_scores_path}{gene}_filtered_scores.csv")
scores.rename(columns={'ID':'sequence'}, inplace=True)

merged_scores=pd.merge(scores, splice_scores, on='sequence', how='left')
merged_scores.to_csv(f"{save_folder}{gene}_scores_variant_features.csv")
