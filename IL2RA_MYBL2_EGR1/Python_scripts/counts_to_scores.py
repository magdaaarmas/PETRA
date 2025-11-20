import sys
import os
import pickle
import numpy as np

########################################################################################################################
#to use locally
path=os.getcwd()
experiment_folder=sys.argv[1]
experiment_name=sys.argv[2]

replicates_to_consider=['1', '3']
########################################################################################################################
def obtain_pseudocounts(df):
    df=df.copy()
    col_names=df.columns.tolist()
    for col in col_names:
        if 'gDNA' in col or 'cDNA' in col:
            df[f"{col}_pseudocounts"]=df[col]+1
    return df

def calculate_relative_frequencies(df):
    df=df.copy()
    columns=df.columns.tolist()
    for column in columns:
        if 'pseudocounts' in column:
            total_pseudocounts=df[column].sum(axis=0)
            df[f"{column}_rel_freq"]=df[column]/total_pseudocounts
    return(df)

def average_across_replicates_general(df, columns):
    df=df.copy()
    columns_to_average=[f"{i}{columns}" for i in replicates_to_consider]
    df[columns+'_average_replicates']=df[columns_to_average].mean(axis=1)
    return df

def add_log2_column(df, column_name):
    df=df.copy()
    # Check if the specified column exists in the DataFrame
    if column_name not in df.columns:
        print(f"Column '{column_name}' not found in DataFrame.")
        return df

    # Add a new column with log2-transformed values
    log2_column_name = f"{column_name}_log2"
    df[log2_column_name] = np.log2(df[column_name])

    return df
def calculate_raw_scores(df, column_to_use):
    df=df.copy()
    columns = df.columns.tolist()
    for i in replicates_to_consider:
            df[f"score_{i}"]=df[f"{i}cDNA_{column_to_use}"]/df[f"{i}gDNA_{column_to_use}"]
            df=add_log2_column(df, f"score_{i}")
    return df

def normalise_to_WT(df):
    for i in replicates_to_consider:
        WT_score=df[df['ID']=='WT'][f"score_{i}_log2"].values[0]
        df[f"score_{i}_log2_norm"]=df[f"score_{i}_log2"]-WT_score
    return df

def average_scores_across_replicates(df):
    df=df.copy()
    columns_to_average=[f"score_{i}_log2_norm" for i in replicates_to_consider]
    df['score_log2_norm_average']=df[columns_to_average].mean(axis=1)
    return df
def counts_to_scores(raw_counts_df_input, save_path):
    raw_counts_df_input=raw_counts_df_input.copy()

    #remove the "total" row
    counts_df=raw_counts_df_input.drop('Total')

    #pseudocounts
    counts_df=obtain_pseudocounts(counts_df)

    #calculate relative frequencies
    counts_df=calculate_relative_frequencies(counts_df)

    #remove WT allele and unmatched alleles
    IDS_to_exclude=['no_extension_found']
    counts_df=counts_df[~counts_df['ID'].isin(IDS_to_exclude)]

    # calculate average gDNA rel_freq for variants
    # calculate average gDNA rel_freq for variants
    counts_df=average_across_replicates_general(counts_df, 'gDNA_count_pseudocounts_rel_freq')
    counts_df = average_across_replicates_general(counts_df, 'cDNA_count_pseudocounts_rel_freq')

    #calculate raw variant score and log2-transform
    counts_df=calculate_raw_scores(counts_df, 'count_pseudocounts_rel_freq')

    #normalise scores to WT
    counts_df=normalise_to_WT(counts_df)

    #average replicate scores
    counts_df=average_scores_across_replicates(counts_df)

    #turn the log2_norm_average score into the "filtered_score" to match with other codes.
    counts_df['filtered_score']=counts_df['score_log2_norm_average']

    #save dataframe with scores
    with open(save_path+gene+'_'+experiment_name+'_scores.pkl', 'wb') as fp:
        pickle.dump(counts_df, fp)
    counts_df.to_csv(save_path+gene+'_'+experiment_name+'_scores.csv')

########################################################################################################################
gene='IL2RA'

#create output directory
output_path=f"{path}/{experiment_folder}/scores/"
os.makedirs(output_path, exist_ok=True)

#obtain counts and number of replicates
with (open(f"{path}/{experiment_folder}/raw_counts/{gene}_raw_counts.pkl", 'rb') as fp):
    counts_dataframe=pickle.load(fp)

#score variants
counts_to_scores(counts_dataframe, output_path)
