import pandas as pd
import matplotlib.pyplot as plt
import pickle
from scipy.stats import gaussian_kde
import scipy
import matplotlib.ticker as ticker
import seaborn as sns
import os
import matplotlib as mpl
from Bio.Seq import Seq

cwd=os.getcwd()
figure_name='1G'
data_source_screen=(f"{cwd}/data/6_mer/jurkat/VAV1_scores_variant_features.csv")
data_source_validation=(f"{cwd}/data/6_mer_validation/VAV1_filtered_scores.csv")
pegRNA_info_file=(f"{cwd}/data/6_mer_validation/VAV1_pegRNA_info.csv")
color_to_use_name_top='VAV1_high'
color_to_use_name_bottom='VAV1_low'
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
color_to_use_top=open_colour_palette(colour_palette, color_to_use_name_top)
color_to_use_bottom=open_colour_palette(colour_palette, color_to_use_name_bottom)
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size
########################################################################################################################
#plot figures
#get data and merge datasets
def obtain_data(validation_data, screen_data, pegRNAinfo):
    # get pegRNA information
    pegRNAs=pd.read_csv(pegRNAinfo)
    pegRNAs=pegRNAs[['sequence', 'category']]
    pegRNAs.rename(columns={'sequence': 'ID'}, inplace=True)

    #get validation data
    validation_df=pd.read_csv(validation_data)
    validation_df=validation_df[['ID', 'filtered_score']]
    validation_df.rename(columns={'filtered_score': f"validation_filtered_score"}, inplace=True)

    #get screen data
    screen_df=pd.read_csv(screen_data)
    screen_df.rename(columns={'sequence': 'ID'}, inplace=True)
    screen_df=screen_df[['ID', 'filtered_score']]
    screen_df.rename(columns={'filtered_score': f"screen_filtered_score"}, inplace=True)

    #merge data
    merged_data=pd.merge(validation_df, screen_df, on='ID', how='left')
    merged_data = pd.merge(merged_data, pegRNAs, on='ID', how='outer')

    return merged_data

data=obtain_data(data_source_validation, data_source_screen, pegRNA_info_file)

def correlation_minipool_screen(df, width, height):
    df=df.copy()
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
    colors={'neutral':greys,
            'top':color_to_use_top,
            'bottom':color_to_use_bottom
    }
    ax=sns.scatterplot(data=df, x='screen_filtered_score', y='validation_filtered_score', hue='category', palette=colors,s=6,
                       edgecolor=None)
    ax.set_xlabel(f"6-mer screen Expression Score")
    ax.set_ylabel(f"Validation library Expression Score")

    # Calculate Spearman correlation
    stats = scipy.stats.spearmanr(df['screen_filtered_score'], df['validation_filtered_score'])
    r_value_text = f"ρ = {stats.correlation:.2f}"
    ax.text(0.95, 0.05, r_value_text, fontsize=font_size, fontname=font_to_use,
            verticalalignment='bottom', horizontalalignment='right', color=greys,
            transform=ax.transAxes)

    # normalise the tick distance in both axes:
    tick_interval = 2
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))

    # Remove top and right borders
    sns.despine()
    plt.legend().remove()
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf')
    plt.show()
    plt.close()

correlation_minipool_screen(data, fig_width_mm, fig_height_mm)