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
import re
from scipy.stats import gaussian_kde


########################################################################################################################
#to use locally
path=os.getcwd()
gene=sys.argv[1]
filtering_info_file=sys.argv[2]
experiment_folder=sys.argv[3]

unfiltered_path=f"{path}/{experiment_folder}/unfiltered_scoring/"
path_filtered=f"{path}/{experiment_folder}/filtered_scores/"
os.makedirs(path_filtered, exist_ok=True)

########################################################################################################################

# define colour scheme
set2_colors = plt.get_cmap('Set2').colors
gene_colour = {'VAV1': set2_colors[0], 'CD28': set2_colors[1], 'IL2RA': set2_colors[2], 'OTUD7B': set2_colors[3]}
def obtain_replicate_number(df):
    df.columns = [col.replace("eg_", "neg_") if "neg_" not in col else col for col in
                  df.columns]
    col_names = df.columns.tolist()
    col_numbers = []
    for col in col_names:
        if 'gDNA_count' in col and 'neg' not in col:
            match = re.match(r"(\d+)", col)  # Regex to match digits at the start of the string
            if match:
                col_number = int(match.group(1))  # Extract the matched number
                col_numbers.append(col_number)
    replicates = max(col_numbers)

    return replicates

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

def plot_correlations_density(df, save_path):
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
                xy = np.vstack([df[columns[i]], df[columns[j]]])
                z = gaussian_kde(xy)(xy)
                ax.scatter(df[columns[i]], df[columns[j]], c=z, s=10, cmap='plasma')
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
            plt.savefig(save_path+gene+"_correlations_replicates_density"+".png", format='png', dpi=300)

        # plt.show()
        plt.close()
    gDNA=[f"{i}gDNA_count_pseudocounts_rel_freq_log2" for i in range(1,replicates+1)]
    cDNA=[f"{i}cDNA_count_pseudocounts_rel_freq_log2" for i in range(1,replicates+1)]
    scores=[f"score_{i}_log2_norm" for i in range(1,replicates+1)]
    plot_variables=[gDNA, cDNA, scores]
    plot_and_correlate(df, plot_variables, save_path)

def filtering_check(df, save_path):
    df=df.copy()
    df.sort_values(by='gDNA_count_pseudocounts_rel_freq_average_replicates', ascending=False, inplace=True, ignore_index=True)
    plt.figure(figsize=(8, 2))
    plt.scatter(df.index, df['score_log2_norm_average_norm'], color=gene_colour[gene], s=10, alpha=0.5)
    plt.title(gene+' gDNA relative frequency vs. variant score')
    plt.xlabel('variants ordered by decreasing gDNA relative frequency')
    plt.ylabel('Variant score')
    plt.grid(False)
    plt.tick_params(labelbottom=False, bottom=False)
    sns.despine(right=True, top=True)
    plt.tight_layout()
    plt.savefig(save_path+gene+'_filtering_visualisation.png', format='png', dpi=300)
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
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.00005))
    plt.xticks(rotation=90)
    plt.grid(False)
    sns.despine(right=True, top=True)
    if thereshold:
        plt.vlines(thereshold, df['score_log2_norm_average_norm'].max(), df['score_log2_norm_average_norm'].min())
    plt.tight_layout()
    plt.savefig(save_path+gene+'_filtering_visualisation_correlation.png', format='png', dpi=300)
    # plt.show()
    plt.close()

def filter_scores(df, rep_number):
    columns_gDNA = [f"{i}gDNA_count_pseudocounts_rel_freq" for i in range(1, rep_number + 1)]
    df = df[(df[columns_gDNA] > filtering_threshold).all(axis=1)]
    df = normalise_average_to_median(df, 'score_log2_norm_average_norm')
    df.rename(columns={'score_log2_norm_average_norm_norm':'filtered_score'}, inplace=True)
    df.reset_index(drop=True, inplace=True)
    with open(output_path + gene + '_filtered_scores.pkl', 'wb') as fp:
        pickle.dump(df, fp)
    df.to_csv(output_path + gene + '_filtered_scores.csv')
    # make plots to decide on filtering

    counts_df_plotting = df.copy()
    output_path_plots = output_path + gene + '_filtered_figures/'
    os.makedirs(output_path_plots, exist_ok=True)
    # calculate log2 of relative frequencies
    for i in range(1, rep_number + 1):
        column1 = f"{i}gDNA_count_pseudocounts_rel_freq"
        column2 = f"{i}cDNA_count_pseudocounts_rel_freq"
        counts_df_plotting = add_log2_column(counts_df_plotting, column1)
        counts_df_plotting = add_log2_column(counts_df_plotting, column2)

    # plot replicate correlations
    plot_correlations(counts_df_plotting, output_path_plots)
    plot_correlations_density(counts_df_plotting, output_path_plots)

    # plot gDNA frequency vs score
    filtering_check(counts_df_plotting, output_path_plots)
    filtering_check_correlations(counts_df_plotting, output_path_plots)

########################################################################################################################
#obtain filtering parameters from .csv file in folder
parameters=pd.read_csv(f"{filtering_info_file}")
parameters=parameters[parameters['gene']==gene]
alignment_filter=parameters['alignment_filter'].values[0]
filtering_threshold=parameters['gDNA_threshold'].values[0]
replicate_removed=parameters['replicate_removed'].values[0]
removed_replicate=parameters['removed_replicate'].values[0]
output_path = (f"{path_filtered}")

#FILTER THE DATA and create output directory
#obtain scores df from unfiltered analysis to select variants to filter out
if replicate_removed:
    #create output_path
    filter_name = str(f"{filtering_threshold:4f}").split('.')[-1]
    os.makedirs(output_path, exist_ok=True)
    #obtain unfiltered scores
    with open(f"{unfiltered_path}/{gene}_unfiltered_scoring_rep_filtered/{gene}_unfiltered_scores.pkl",
              'rb') as fp:
        scores_unfiltered = pickle.load(fp)
    replicates = (obtain_replicate_number(scores_unfiltered))


else:
    # create output_path
    filter_name = str(f"{filtering_threshold:4f}").split('.')[-1]
    os.makedirs(output_path, exist_ok=True)
    #obtain unfiltered scores
    with open(f"{unfiltered_path}/{gene}_unfiltered_scores.pkl", 'rb') as fp:
        scores_unfiltered=pickle.load(fp)
    replicates = (obtain_replicate_number(scores_unfiltered))


#score the variants of the filtered data
filter_scores(scores_unfiltered, replicates)