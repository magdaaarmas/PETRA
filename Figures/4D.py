import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
import os
import scipy
import matplotlib.ticker as ticker

cwd=os.getcwd()
figure_name='4D'
gene='IL2RA'
primary_data_directory=(f"{cwd}/data/TF_tandem/{gene}_CD3_scores.csv")
jurkat_data_directory=(f"{cwd}/data/TF_tandem/{gene}_jurkat_scores.csv")
fig_width_mm=50
fig_height_mm=50

########################################################################################################################
figures_path=f"{cwd}/output_figures/"
colour_palette=(f"{cwd}/color_palette.csv")

def mm_to_inches(w, h):
    return w / 25.4, h / 25.4
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
def merge_scores(name_1, type_1, name_2, type_2):
    def modify_df(df, cell_type):
        df.rename(columns={'ID':'sequence'}, inplace=True)
        df = df.copy()
        df = df[['sequence', 'origin', 'category', 'score_log2_norm_average_norm']]
        df.rename(columns={'score_log2_norm_average_norm': f"score_{cell_type}"}, inplace=True)
        return df
    type1=modify_df(type_1, name_1)
    type2=modify_df(type_2, name_2)

    remove_columns=['category', 'origin']
    type2=type2.drop(remove_columns, axis=1)
    df=pd.merge(type1, type2, on='sequence', how='inner')
    return df
def expand_categories(df):
    df = df.copy()
    df['category'] = np.where(df['origin'].str.contains('sequence_from'), '6_bp_TF_coding', df['category'])
    df['category'] = np.where(df['origin'].str.contains('motif_from'), 'TF_motif_x1', df['category'])
    df['category'] = np.where(df['origin'].str.contains('homodimer_from'), 'TF_motif_x2', df['category'])
    df['category'] = np.where(df['origin'].str.contains('homotrimer_from'), 'TF_motif_x3', df['category'])
    if gene == 'IL2RA':
        df['category'] = np.where(df['origin'].str.contains('top_sequence_screen_jurkat'), 'other_variant',
                                  df['category'])
    elif gene == 'VAV1':
        df['category'] = np.where(df['origin'].str.contains('top_sequence_screen'), 'other_variant',
                                  df['category'])
    df['category'] = np.where(df['category'] == 'variant', 'other_variant', df['category'])
    df['category'] = np.where(df['category'] == 'negative_neutral_SNV', 'other_variant', df['category'])
    df['category'] = np.where(df['category'] == 'splice', 'other_variant', df['category'])

    return df
def plot_correlation(dataframe_to_use, width, height):
    df = dataframe_to_use.copy()
    df=df[df['sequence']!='WT']

    palette = {
        'negative_neutral_sequence': 'orchid',
        'negative_neutral_SNV': 'tomato',
        'other_variant': 'lightgrey',
        'splice': 'red',
        '6_bp_TF_coding': 'lightcoral',
        'TF_motif_x1': 'deepskyblue',
        'TF_motif_x2': 'dodgerblue',
        'TF_motif_x3': 'darkblue',
        '6_bp_positive_control': 'limegreen'
    }

    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))

    # Scatterplot
    sns.scatterplot(
        data=df[df['category']=='variant'],
        y='score_jurkat',
        x='score_primary',
        s=10,
        ax=ax,
        hue='category',
        palette=palette,
        edgecolor=None
    )
    sns.scatterplot(
        data=df[df['category']!='variant'],
        y='score_jurkat',
        x='score_primary',
        s=10,
        ax=ax,
        hue='category',
        palette=palette,
        edgecolor=None
    )
    # Calculate Spearman correlation
    stats = scipy.stats.spearmanr(df['score_jurkat'], df['score_primary'])
    r_value_text = f"ρ = {stats.correlation:.2f}"
    mantissa, exponent = f"{stats.pvalue:.2e}".split("e")
    exponent = int(exponent)
    p_value_text = f"$p$ = {mantissa} × 10$^{{{exponent}}}$"
    ax.text(0.95, 0.1, f"{r_value_text}\n{p_value_text}", fontsize=font_size, fontname=font_to_use,
            verticalalignment='bottom', horizontalalignment='right', color=greys,
            transform=ax.transAxes)

    # Find global min and max to cover both axes equally
    min_val = min(df['score_jurkat'].min() - 2, df['score_primary'].min() - 2)
    max_val = max(df['score_jurkat'].max() + 2, df['score_primary'].max() + 2)

    # Plot the x=y line
    ax.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='lightgrey', linewidth=0.7)

    # Set limits
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    # Force same number format (1 decimal place) on both axes
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

    # Force equal aspect ratio
    ax.set_aspect('equal')

    # Labels and formatting
    plt.ylabel(f"Jurkat Expression score")
    plt.xlabel(f"Primary Expression score")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    sns.despine(top=True, right=True)
    plt.legend().set_visible(False)
    plt.title(r"$\mathit{" + gene + "}$")
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()

#obtain information for plotting from master file
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)
color_to_use_name = gene
color_to_use = open_colour_palette(colour_palette, color_to_use_name)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################
jurkat_data = pd.read_csv(jurkat_data_directory)
jurkat_data.drop(columns=['Unnamed: 0'], inplace=True)

primary_data = pd.read_csv(primary_data_directory)
primary_data.drop(columns=['Unnamed: 0'], inplace=True)

#merge jurkat and primary scores
df=merge_scores('jurkat', jurkat_data, 'primary', primary_data)

#expand variant categories for plotting
df = expand_categories(df)

#plot
plot_correlation(df, fig_width_mm, fig_height_mm)
