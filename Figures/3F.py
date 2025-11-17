import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
from statsmodels.stats.multitest import multipletests
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl


cwd=os.getcwd()
figure_name='3F'
timepoints=['7', '15']
fig_width_mm=77
fig_height_mm=30
gene='IL2RA'

def mm_to_inches(w, h):
    return w / 25.4, h / 25.4
def compute_statistics(df, score_type):
    df = df.copy()

    # remove negative, WT and neutral samples for statistical comparison to neutral
    samples = list(set(df['name_for_stats'].tolist()))
    samples = [sample for sample in samples if 'neg' not in sample]
    samples = [sample for sample in samples if 'neutral' not in sample]
    samples = [sample for sample in samples if 'WT' not in sample]

    # determine score column
    score_column = f"{score_type}_log2_norm"

    # calculate statistics
    stats = {}
    neutral_data = df[df['name_for_stats'] == 'neutral'][score_column].tolist()
    for sample in samples:
        sample_data = df[df['name_for_stats'] == sample][score_column].tolist()
        stat, p_value = ttest_ind(neutral_data, sample_data, alternative='two-sided')
        stats[f"{sample}_vs_neutral"] = p_value

    # combine stats for all samples in a df
    stats_df = pd.DataFrame(list(stats.items()), columns=['Groups', 'ttest_p_value'])

    # corrections to control for FDR
    reject_bonferroni, pvals_corrected_bonferroni, _, _ = multipletests(stats_df["ttest_p_value"],
                                                                        alpha=0.05,
                                                                        method='bonferroni')
    reject_BH, pvals_corrected_BH, _, _ = multipletests(stats_df["ttest_p_value"],
                                                        alpha=0.05,
                                                        method='fdr_bh')

    # add corrections to df
    stats_df["p_bonferroni"] = pvals_corrected_bonferroni
    stats_df["sig_bonferroni"] = reject_bonferroni
    stats_df["p_BH"] = pvals_corrected_BH
    stats_df["sig_BH"] = reject_BH

    return stats_df
def open_colour_palette(file, item_to_search=False):
    palette_df=pd.read_csv(file)
    if item_to_search:
        color=palette_df[palette_df['use']==item_to_search]['color_code'].values[0]
        return color
    else:
        return palette_df
def obtain_font(file):
    palette_df=pd.read_csv(file)

    font=palette_df[palette_df['use']=='font']['color_name'].values[0]
    size=int(palette_df[palette_df['use']=='font_size']['color_name'].values[0])
    return font, size
def plot_flow_scores(df, column, y_col_name, width, height):
    df = df.copy()
    column_to_plot = [i for i in df.columns if column in i and 'log2_norm' in i]

    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))

    #define order and sort
    custom_order = ['neutral_AACTGT', 'neutral_CATTAG', 'MYBL2_homodimer', 'sequence_CGTTAT', 'sequence_GGTTTA',
                    'TFAP2A_homodimer', 'TFAP2A_homotrimer', 'sequence_ATGGAT']
    df["name_to_plot"] = pd.Categorical(df["name_to_plot"], categories=custom_order, ordered=True)
    df = df.sort_values('name_to_plot')

    #generate plot
    sns.barplot(data=df, y='name_to_plot', x=column_to_plot[0], ax=ax, hue='type',palette={'neutral':'lightgrey',
                                                                                           'positive':open_colour_palette(colour_palette, f"IL2RA")},
                errorbar=None, legend=None)
    sns.stripplot(data=df, y='name_to_plot', x=column_to_plot[0], ax=ax, color='dimgrey', s=3, legend=None)

    #rename variants for plotting
    relabel={
        'MYBL2_homodimer':'MYBL2 motif x2',
        'TFAP2A_homodimer':'TFAP2A motif x2',
        'TFAP2A_homotrimer': 'TFAP2A motif x3',
        'neutral_AACTGT':'Neutral - AACTGT',
        'neutral_CATTAG': 'Neutral - CATTAG',
        'sequence_CGTTAT':'Variant CGTTAT',
        'sequence_GGTTTA': 'Variant GGTTTA',
        'sequence_ATGGAT': 'Variant ATGGAT'
    }

        # Apply the mapping to existing tick labels
    current_labels = [tick.get_text() for tick in ax.get_yticklabels()]
    new_labels = [relabel.get(label, label) for label in current_labels]  # default to original if not in map

    ax.set_yticklabels(new_labels, rotation=0, ha='right')
    plt.ylabel('')
    plt.xlabel(y_col_name)

    sns.despine(top=True, right=True)
    sns.despine(top=True, right=True)
    plt.tight_layout()

    # plt.savefig(f"{figures_path}{figure_name}_day_{timepoint}_{y_col_name}.pdf", transparent=True, format='pdf')
    plt.show()
    plt.close()

figures_path=f"{cwd}/output_figures/"
colour_palette=(f"{cwd}/color_palette.csv")

#obtain information for plotting from master file
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)
color_to_use_negatives=greys

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################

for timepoint in timepoints:
    data_path = f"{cwd}/data/Flow_cytometry/flow_data_day_{timepoint}.csv"
    scores = pd.read_csv(data_path)

    #rename columns
    # rename columns
    total_MFI_column = 'Lymphocytes/Single Cells | Mean (FL2-A :: IL2RA: PE-A)'
    positive_MFI_column = 'Lymphocytes/Single Cells/IL2RA: PE-A, FITC-A subset | Mean (FL2-A :: IL2RA: PE-A)'
    frequency_column = 'Lymphocytes/Single Cells/IL2RA: PE-A, FITC-A subset | Freq. of Total (%)'
    scores.rename(columns={'Unnamed: 0': 'sample', total_MFI_column: 'MFI_total',
                           positive_MFI_column: 'MFI_positive',
                           frequency_column: 'frequency_positive'}, inplace=True)
    #remove unnecessary rows
    scores = scores[scores['sample'] != 'Mean']
    scores = scores[scores['sample'] != 'SD']

    # merge all neutral scores into a single "name_for_stats" category
    scores['name_for_stats'] = np.where(scores['name_to_plot'].str.contains('neutral'), 'neutral',
                                            scores['name_to_plot'])

    # calculate the scores
    score_columns = ['MFI_total', 'MFI_positive', 'frequency_positive']
    statistics = {}
    for score in score_columns:
        # log2 transform the scores
        scores[f"{score}_log2"] = scores[score].apply(lambda x: np.log2(x))
        # find the mean of all neutral measurements
        neutral_mean = scores[scores['type'] == 'neutral'][f"{score}_log2"].mean()
        # normalise to the mean from neutral samples
        scores[f"{score}_log2_norm"] = scores[f"{score}_log2"] - neutral_mean
        #calculate statistics
        statistics[score]=compute_statistics(scores, score)
        #remove WT and negatives for plotting
        scores = scores[~scores['type'].isin(['WT', 'neg'])]
        #make plot
        plot_flow_scores(scores, 'MFI_total', 'MFI (log2)', fig_width_mm, fig_height_mm)

