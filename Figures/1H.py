import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import os
import matplotlib as mpl
import scipy


cwd=os.getcwd()
figure_name='1H'
data_source=(f"{cwd}/data/6_mer/jurkat/VAV1_scores_variant_features.csv")
color_to_use_name='VAV1'
fig_width_mm=60
fig_height_mm=60

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

#obtain information for plotting from master file
color_to_use=open_colour_palette(colour_palette, color_to_use_name)
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

def correlate_ES_SAI(df, width, height, SAI_column, interval_ES, interval_SAI, save_path=None):
    df=df.copy()
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))

    #plot
    ax.scatter(df[SAI_column],df['filtered_score'] , color=color_to_use, s=3)
    ax.set_xlabel(f"SpliceAI Donor Gain Score")
    ax.set_ylabel(f"Expression Score")

    #normalise the tick distance in both axes:
    ax.yaxis.set_major_locator(ticker.MultipleLocator(interval_ES))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(interval_SAI))

    #add correlation
    # Calculate Spearman correlation
    stats = scipy.stats.spearmanr(df[SAI_column], df['filtered_score'])
    r_value_text = f"ρ = {stats.correlation:.2f}"
    ax.text(0.95, 0.05, r_value_text, fontsize=font_size, fontname=font_to_use,
            verticalalignment='bottom', horizontalalignment='right', color=greys,
            transform=ax.transAxes)

    # Remove top and right borders
    sns.despine()
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf')
    plt.show()
    plt.close()

correlate_ES_SAI(data, fig_width_mm, fig_height_mm, 'DS_DG', 2, 0.2, save_path=figures_path)
