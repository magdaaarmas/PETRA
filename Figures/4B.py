import pandas as pd
import matplotlib.pyplot as plt
import pickle
import numpy as np
from scipy.stats import gaussian_kde
import scipy
import matplotlib.ticker as ticker
import seaborn as sns
import os
import matplotlib as mpl
from functools import reduce
import math

cwd=os.getcwd()
figure_name='4B'
data_paths={
    'primary':(f"{cwd}/data/6_mer/primary/IL2RA_scores_variant_features.csv"),
    'jurkat':(f"{cwd}/data/6_mer/jurkat/IL2RA_scores_variant_features.csv")
}
color_to_use_name=f"IL2RA_continuous"
fig_width_mm=60
fig_height_mm=60
reverse_gradient=True
gene='IL2RA'

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
def open_pickle(file):
    with open(file, 'rb') as fp:
        df=pickle.load(fp)
    return df
def reverse_scale(string):
    string=f"{string}_r"
    return string

#obtain information for plotting from master file
color_to_use=open_colour_palette(colour_palette, color_to_use_name)
if reverse_gradient==True:
    color_to_use=reverse_scale(color_to_use)
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)
color_to_use_negatives=greys


#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################
#merge jurkat and primary data with significance
data_dict={}
for i in ['primary', 'jurkat']:
    data=pd.read_csv(data_paths[i])
    data=data[['filtered_score', 'sequence']]
    data.rename(columns={'filtered_score':f"{i}_score"}, inplace=True)
    data_dict[i]=data
merged_df=reduce(lambda left, right: pd.merge(left, right, on='sequence', how='inner'), data_dict.values())

def plot_correlation(df, col1, col2, width, height, interval):

    df=df.copy()
    fig, ax =  plt.subplots(figsize=mm_to_inches(width, height))

    df_2=df[[col1, col2]]
    df_2 = df_2[np.isfinite(df_2).all(axis=1)]

    xy = np.vstack([df_2[col1], df_2[col2]])
    z = gaussian_kde(xy)(xy)

    # Scatterplot
    sns.scatterplot(
        data=df,
        y=col1,
        x=col2,
        s=10,
        ax=ax,
        c=z,
        cmap=color_to_use,
        edgecolor=None,
        zorder=1
    )

    # Calculate Spearman correlation
    stats = scipy.stats.spearmanr(df[col1], df[col2])
    if stats.pvalue <= 0.01:
        pvalue_text = 'p≤0.01'
    else:
        pvalue_text = f"p={stats.pvalue:.2f}"

    r_value_text = f"ρ = {stats.correlation:.2f}"
    pval = stats.pvalue
    exp = int(math.floor(math.log10(pval)))
    mantissa = pval / 10 ** exp

    p_value_text = rf"$p = {mantissa:.2f} \times 10^{{{exp}}}$"
    ax.text(0.95, 0.05, f"{r_value_text}\n{p_value_text}", fontsize=font_size, fontname=font_to_use,
            verticalalignment='bottom', horizontalalignment='right', color=greys,
            transform=ax.transAxes)

    # # Find global min and max to cover both axes equally
    min_val = min(df[col1].min() - 2, df[col2].min() - 2)
    max_val = max(df[col1].max() + 2, df[col2].max() + 2)

    #normalise the tick distance in both axes:
    tick_interval=interval
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))

    # Plot the x=y line
    ax.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='dimgrey', linewidth=0.5)
    #
    # Set limits
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)

    # Force equal aspect ratio
    ax.set_aspect('equal')

    # Labels and formatting
    plt.ylabel(f"{col1}")
    plt.xlabel(f"{col2}")
    sns.despine()
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", transparent=True, format='pdf')
    plt.show()
    plt.close()

plot_correlation(merged_df, f"jurkat_score", f"primary_score", fig_width_mm, fig_height_mm, 2)
