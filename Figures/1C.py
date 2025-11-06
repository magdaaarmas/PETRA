import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy
import matplotlib.ticker as ticker
import seaborn as sns
import matplotlib as mpl

cwd=os.getcwd()
figure_name='1C'
data_source=(f"{cwd}/data/6_mer/jurkat/VAV1_scores_variant_features.csv")
color_to_use_name='VAV1_continuous'
fig_width_mm=45
fig_height_mm=45

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
def reverse_scale(string):
    string=f"{string}_r"
    return string

#obtain information for plotting from master file
color_to_use=open_colour_palette(colour_palette, color_to_use_name)
color_to_use=reverse_scale(color_to_use)
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size
########################################################################################################################
#plot figures
#get data
data=pd.read_csv(data_source)

def create_log10_columns(df):
    df=df.copy()
    for i in range(1,3):
        df[f"{i}gDNA_plot"]=df[f"{i}gDNA_count_pseudocounts_rel_freq"].apply(lambda x: np.log10(x))
        df[f"{i}cDNA_plot"]=df[f"{i}cDNA_count_pseudocounts_rel_freq"].apply(lambda x: np.log10(x))
    return df

data=create_log10_columns(data)

def plot_and_correlate(df, width, height, rep1, rep2, column, axis_label, interval, save_path=None):
    df=df.copy()
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
    #define column to plot
    if column=='score':
        column1=f"score_{str(rep1)}_log2_norm"
        column2=f"score_{str(rep2)}_log2_norm"
    else:
        column1=f"{rep1}{column}"
        column2=f"{rep2}{column}"

    xy = np.vstack([df[column1], df[column2]])
    z = gaussian_kde(xy)(xy)

    #plot
    ax.scatter(df[column1], df[column2], c=z, s=3, cmap=color_to_use)
    ax.set_xlabel(f"{axis_label} - Rep. {rep1}")
    ax.set_ylabel(f"{axis_label} - Rep. {rep2}")

    # Calculate Spearman correlation
    stats = scipy.stats.spearmanr(df[column1], df[column2])

    r_value_text = f"ρ = {stats.correlation:.2f}"
    ax.text(0.95, 0.05, r_value_text, fontsize=font_size, fontname=font_to_use,
            verticalalignment='bottom', horizontalalignment='right', color=greys,
            transform=ax.transAxes)


    #normalise the tick distance in both axes:
    tick_interval=interval
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))


    # Remove top and right borders
    sns.despine()
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf')
    plt.show()
    plt.close()

plot_and_correlate(data, fig_width_mm, fig_height_mm, 1, 2, 'gDNA_plot',
                   'gDNA variant freq. (log10)', 0.5)

