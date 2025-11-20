import sys
import os
import pandas as pd
import pickle
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
########################################################################################################################
#to use locally
path=os.getcwd()
gene=sys.argv[1]
cell_type=sys.argv[2]
experiment_folder=sys.argv[3]
pegRNA_information_file=sys.argv[4]

########################################################################################################################
def obtain_replicate_number(df):
    col_names=df.columns.tolist()
    col_numbers=[]
    for col in col_names:
        if 'gDNA' in col and 'neg' not in col:
            col=col[0]
            col_numbers.append(int(col))
    replicates=max(col_numbers)
    return replicates
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
    columns_to_average=[f"{i}{columns}" for i in range(1,replicates+1)]
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
    for i in range(1, replicates+1):
            df[f"score_{i}"]=df[f"{i}cDNA_{column_to_use}"]/df[f"{i}gDNA_{column_to_use}"]
            df=add_log2_column(df, f"score_{i}")
    return df
def create_categories(df):
    # obtain variant information
    variant_information = pd.read_csv(f"{pegRNA_information_file}")
    variant_information.rename(columns={'sequence': 'ID'}, inplace=True)
    variant_information = variant_information[['ID', 'origin']]
    df = pd.merge(df, variant_information, on='ID', how='outer')
    df['origin'] = np.where(df['ID'] == 'WT', 'WT', df['origin'])

    df = df.copy()
    df['category'] = np.where(df['origin'].str.contains('WT'), 'WT',
                              np.where(df['origin'].str.contains('splice'), 'splice',
                              np.where(df['origin'].str.contains('negative'), df['origin'], 'variant')))
    df=df.dropna()
    return df
def normalize_to_neutral(df):
    replicate_columns=[f"score_{i}_log2" for i in range(1, replicates+1)]
    for column in replicate_columns:
        neutral_average=df[df['category']=='negative_neutral_sequence'][column].mean(axis=0)
        df[f'{column}_norm'] = df[column] - neutral_average
    return df
def average_scores_across_replicates(df):
    df=df.copy()
    columns_to_average=[f"score_{i}_log2_norm" for i in range(1,replicates+1)]
    df['score_log2_norm_average']=df[columns_to_average].mean(axis=1)
    return df
def normalise_average_neutral(df, column_name):
    df=df.copy()
    average_neutral=df[df['category']=='negative_neutral_sequence'][column_name].mean()
    df[column_name+'_norm']=df[column_name]-average_neutral
    return df
def statistics_welsh_t_test_to_neutral(df):
    df = df.copy()
    columns_to_use = [f"score_{i}_log2_norm" for i in range(1, replicates + 1)]
    df_melted = df.melt(id_vars=['ID', 'category'], value_vars=columns_to_use)

    # obtain neutral data
    neutral_data = df_melted[df_melted['category'] == 'negative_neutral_sequence']['value']
    results = []

    # Iterate over all non-neutral variants
    non_neutral_IDs = df_melted[df_melted['category'] != 'negative_neutral_sequence']['ID'].tolist()

    for var in df['ID'].unique():
        if var not in non_neutral_IDs:
            continue
        elif var in non_neutral_IDs:
            variant_data = df_melted[df_melted['ID'] == var]['value']

            # Welch's t-test
            t_stat, p_val = ttest_ind(variant_data, neutral_data, equal_var=False)

            results.append({
                'ID': var,
                't_statistic': t_stat,
                'p_value': p_val
            })
    results_df = pd.DataFrame(results)

    # perform BH-correction
    pvals = results_df['p_value'].values
    rejected, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')

    # Add to your results dataframe
    results_df['q_val'] = pvals_corrected

    # annotate significance
    results_df['sig_05'] = np.where(results_df['q_val'] <= 0.05, True, False)
    results_df['sig_01'] = np.where(results_df['q_val'] <= 0.01, True, False)

    # merge statistics with screen scores
    results_with_stats = pd.merge(df, results_df, on='ID', how='outer')

    return results_with_stats

def counts_to_scores(raw_counts_df_input):
    raw_counts_df_input = raw_counts_df_input.copy()

    # remove the "total" row
    counts_df = raw_counts_df_input.drop('Total')

    # pseudocounts
    counts_df = obtain_pseudocounts(counts_df)

    # calculate relative frequencies
    counts_df = calculate_relative_frequencies(counts_df)

    # remove WT allele and unmatched alleles
    IDS_to_exclude = ['WT', 'no_extension_found']
    counts_df = counts_df[~counts_df['ID'].isin(IDS_to_exclude)]

    # calculate average gDNA rel_freq for variants
    # calculate average gDNA rel_freq for variants
    counts_df = average_across_replicates_general(counts_df, 'gDNA_count_pseudocounts_rel_freq')
    counts_df = average_across_replicates_general(counts_df, 'cDNA_count_pseudocounts_rel_freq')

    # calculate raw variant score and log2-transform
    counts_df = calculate_raw_scores(counts_df, 'count_pseudocounts_rel_freq')

    #create the categories
    counts_df = create_categories(counts_df)

    #normalise to neutral variants
    counts_df = normalize_to_neutral(counts_df)

    #average scores across replicates
    counts_df = average_scores_across_replicates(counts_df)

    #normalise to the neutral average
    counts_df=normalise_average_neutral(counts_df, 'score_log2_norm_average')

    #statistics
    counts_df=statistics_welsh_t_test_to_neutral(counts_df)

    # save dataframe with scores_unfiltered
    with open(f"{output_path}scores.pkl", 'wb') as fp:
        pickle.dump(counts_df, fp)
    counts_df.to_csv(f"{output_path}scores.csv")

########################################################################################################################
sample_path=f"{path}/{experiment_folder}/{cell_type}/{gene}"
output_path=f"{sample_path}/scores/"
os.makedirs(output_path, exist_ok=True)

#obtain counts and number of replicates
with (open(f"{sample_path}/raw_counts/raw_counts.pkl", 'rb') as fp):
    counts_dataframe=pickle.load(fp)
replicates=obtain_replicate_number(counts_dataframe)

#score variants
scores=counts_to_scores(counts_dataframe)
