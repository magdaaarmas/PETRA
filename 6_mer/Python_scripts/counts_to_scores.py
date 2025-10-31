import sys
import os
import pandas as pd
from functools import reduce
from itertools import product
import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy
import matplotlib.ticker as ticker
import seaborn as sns

########################################################################################################################
#to use locally
path=os.getcwd()
gene=sys.argv[1]
alignment_filter=sys.argv[2]
filtering_info_file=sys.argv[3]

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
def normalize_to_median(df):
    replicate_columns=[f"score_{i}_log2" for i in range(1, replicates+1)]
    for column in replicate_columns:
        median_replicate = np.median(df[column])
        df[f'{column}_norm'] = df[column] - median_replicate
    return df
def average_scores_across_replicates(df):
    df=df.copy()
    columns_to_average=[f"score_{i}_log2_norm" for i in range(1,replicates+1)]
    df['score_log2_norm_average']=df[columns_to_average].mean(axis=1)
    return df
def normalise_average_to_median(df, column_name):
    df=df.copy()
    average_median=df[column_name].median()
    df[column_name+'_norm']=df[column_name]-average_median
    return df

def plot_correlations(df, save_path):
    df=df.copy()

    def plot_and_correlate(df, columns_list, save_path=None):
        df=df.copy()
        fig, axes = plt.subplots(3, 3, figsize=(15, 15),constrained_layout=True)

        # Define combinations of column indices
        if replicates==2:
            combinations = [(0, 1)]
        elif replicates==3:
            combinations = [(0, 1), (1, 2), (0, 2)]

        for row_index, columns in enumerate(columns_list):
            titles = [f"{columns[i]} vs {columns[j]}" for i, j in combinations]

            for ax, (i, j), title in zip(axes[row_index], combinations, titles):
                ax.scatter(df[columns[i]], df[columns[j]], alpha=0.5, color=gene_colour[gene], s=15)
                ax.set_xlabel(columns[i])
                ax.set_ylabel(columns[j])

                # Calculate Spearman correlation
                stats = scipy.stats.spearmanr(df[columns[i]], df[columns[j]])
                r_value_text = f"R={stats.correlation:.2f}"
                ax.text(0.95, 0.05, r_value_text, fontsize=14,
                        verticalalignment='bottom', horizontalalignment='right', color='darkgrey',
                        transform=ax.transAxes)

                # Remove top and right borders
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                # Format x-axis to scientific notation
                ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style='scientific', axis='x', scilimits=(0, 0))
            fig.subplots_adjust()
        if save_path:
            # Ensure the directory exists
            os.makedirs(save_path, exist_ok=True)
            # Replace "/" with "_" in the title for the filename
            plt.savefig(save_path+gene+"_correlations_replicates"+".png", format='png', dpi=300)

        # plt.show()
        plt.close()
    gDNA=[f"{i}gDNA_count_pseudocounts_rel_freq_log2" for i in range(1,replicates+1)]
    cDNA=[f"{i}cDNA_count_pseudocounts_rel_freq_log2" for i in range(1,replicates+1)]
    scores=[f"score_{i}_log2_norm" for i in range(1,replicates+1)]
    plot_variables=[gDNA, cDNA, scores]
    plot_and_correlate(df, plot_variables, save_path)

def filtering_check(df, save_path,thereshold=None):
    df=df.copy()
    df.sort_values(by='gDNA_count_pseudocounts_rel_freq_average_replicates', ascending=False, inplace=True, ignore_index=True)
    df['filter']=np.where(df['gDNA_count_pseudocounts_rel_freq_average_replicates']<thereshold, True, False)
    pelette={False:gene_colour[gene], True:'red'}
    fig, ax=plt.subplots(figsize=(8, 2))
    ax=sns.scatterplot(data=df, x=df.index, y='score_log2_norm_average_norm', hue='filter', palette=pelette, s=10)
    plt.title(gene+' gDNA relative frequency vs. variant score')
    plt.xlabel('variants ordered by decreasing gDNA relative frequency')
    plt.ylabel('Variant score')
    plt.grid(False)
    plt.tick_params(labelbottom=False, bottom=False)
    sns.despine(right=True, top=True)
    plt.tight_layout()
    threshold_name=str(f"{thereshold:4f}").split('.')[-1]
    plt.savefig(save_path+gene+'_filtering_visualisation'+threshold_name+'.png', format='png', dpi=300)
    # plt.show()
    plt.close()

def filtering_check_correlations(df, save_path, thereshold=None):
    df=df.copy()
    fig, ax=plt.subplots(figsize=(8, 4))
    ax.scatter(df['gDNA_count_pseudocounts_rel_freq_average_replicates'], df['score_log2_norm_average_norm'],
                color=gene_colour[gene], s=10, alpha=0.5)
    plt.title(gene+' gDNA relative frequency vs. variant score')
    plt.xlabel('gDNA rel freq replicate average')
    plt.ylabel('Variant score')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))  # Tick every 0.5 units
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.00001))
    plt.xticks(rotation=90)
    plt.grid(False)
    sns.despine(right=True, top=True)
    if thereshold:
        plt.vlines(thereshold, df['score_log2_norm_average_norm'].max(), df['score_log2_norm_average_norm'].min())
    plt.tight_layout()
    threshold_name=str(f"{thereshold:4f}").split('.')[-1]
    plt.savefig(save_path+gene+'_filtering_visualisation_correlation_'+threshold_name+'.png', format='png', dpi=300)
    # plt.show()
    plt.close()
def counts_to_scores(raw_counts_df_input, save_path):
    raw_counts_df_input=raw_counts_df_input.copy()

    #remove the "total" row
    counts_df=raw_counts_df_input.drop('Total')

    #pseudocounts
    counts_df=obtain_pseudocounts(counts_df)

    #calculate relative frequencies
    counts_df=calculate_relative_frequencies(counts_df)

    #remove WT allele and unmatched alleles
    IDS_to_exclude=['WT', 'no_extension_found']
    counts_df=counts_df[~counts_df['ID'].isin(IDS_to_exclude)]

    # calculate average gDNA rel_freq for variants
    # calculate average gDNA rel_freq for variants
    counts_df=average_across_replicates_general(counts_df, 'gDNA_count_pseudocounts_rel_freq')
    counts_df = average_across_replicates_general(counts_df, 'cDNA_count_pseudocounts_rel_freq')

    #calculate raw variant score and log2-transform
    counts_df=calculate_raw_scores(counts_df, 'count_pseudocounts_rel_freq')

    #intrareplicate normalization to median
    counts_df=normalize_to_median(counts_df)

    #average replicate scores
    counts_df=average_scores_across_replicates(counts_df)

    #normalise averaged score to median
    counts_df=normalise_average_to_median(counts_df, 'score_log2_norm_average')

    #save dataframe with scores_unfiltered
    with open(save_path+gene+'_unfiltered_scores.pkl', 'wb') as fp:
        pickle.dump(counts_df, fp)
    counts_df.to_csv(save_path+gene+'_unfiltered_scores.csv')

    ########################################################################################################################
    #make plots to decide on filtering
    counts_df_plotting=counts_df.copy()

    output_path_plots=save_path+gene+'_unfiltered_figures/'
    os.makedirs(output_path_plots, exist_ok=True)

    #calculate log2 of relative frequencies
    for i in range(1,replicates+1):
        column1=f"{i}gDNA_count_pseudocounts_rel_freq"
        column2 = f"{i}cDNA_count_pseudocounts_rel_freq"
        counts_df_plotting=add_log2_column(counts_df_plotting, column1)
        counts_df_plotting = add_log2_column(counts_df_plotting, column2)

    #plot replicate correlations
    plot_correlations(counts_df_plotting, output_path_plots)

    #plot gDNA frequency vs score
    filtering_check(counts_df_plotting, output_path_plots, tentative_filtering_thereshold)
    filtering_check_correlations(counts_df_plotting, output_path_plots, tentative_filtering_thereshold)


########################################################################################################################
#create output directory
output_path=f"{path}/unfiltered_scoring/"
os.makedirs(output_path, exist_ok=True)

#define colour scheme
set2_colors = plt.get_cmap('Set2').colors
gene_colour={'VAV1':set2_colors[0], 'CD28':set2_colors[1], 'IL2RA':set2_colors[2], 'OTUD7B':set2_colors[3]}


#obtain counts and number of replicates
with (open(f"{path}/raw_counts/{gene}_raw_counts.pkl", 'rb') as fp):
    counts_dataframe=pickle.load(fp)
#neg samples are named "eg" in IL2RA and the neg from IL2RA have sample number, replace column names to normalise for all screens
counts_dataframe.columns = [col.replace("eg_", "neg_") if "neg_" not in col else col for col in counts_dataframe.columns]
counts_dataframe.columns = [col.replace("neg_1gDNA", "neg_gDNA") if "neg_1gDNA" in col else col for col in counts_dataframe.columns]
counts_dataframe.columns = [col.replace("neg_1cDNA", "neg_cDNA") if "neg_1cDNA" in col else col for col in counts_dataframe.columns]

replicates=obtain_replicate_number(counts_dataframe)

parameters=pd.read_csv(f"{filtering_info_file}")
parameters=parameters[parameters['gene']==gene]
tentative_filtering_thereshold=parameters['gDNA_threshold'].values[0]
remove_replicate=parameters['replicate_removed'].values[0]
replicate_to_remove=parameters['removed_replicate'].values[0]

#score variants
counts_to_scores(counts_dataframe, output_path)

#score variants with replicate removed
if remove_replicate:
    output_path_no_rep = f"{output_path}/{gene}_unfiltered_scoring_rep_filtered/"
    os.makedirs(output_path_no_rep, exist_ok=True)
    columns=counts_dataframe.columns.tolist()
    columns_to_keep=[column for column in columns if str(replicate_to_remove) not in column]
    counts_dataframe_rep_removed=counts_dataframe[columns_to_keep]
    replicates = obtain_replicate_number(counts_dataframe_rep_removed)
    counts_to_scores(counts_dataframe_rep_removed, output_path_no_rep)