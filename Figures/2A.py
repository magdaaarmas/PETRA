import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import seaborn as sns
import os
import matplotlib as mpl

cwd=os.getcwd()
figure_name='2A'
genes=['VAV1', 'CD28', 'IL2RA', 'OTUD7B']
data_directory=(f"{cwd}/data/6_mer/jurkat/")
files_name='_scores_variant_features.csv'
fig_width_mm=87
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
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

#define palette
# plot
gene_specific_palette = {'VAV1_low':open_colour_palette(colour_palette, 'VAV1_low'),
           'VAV1_high':open_colour_palette(colour_palette, 'VAV1_high'),
           'CD28_low': open_colour_palette(colour_palette, 'CD28_low'),
           'CD28_high': open_colour_palette(colour_palette, 'CD28_high'),
           'OTUD7B_low': open_colour_palette(colour_palette, 'OTUD7B_low'),
           'OTUD7B_high': open_colour_palette(colour_palette, 'OTUD7B_high'),
           'IL2RA_high': open_colour_palette(colour_palette, 'IL2RA_high'),
           'IL2RA_low': open_colour_palette(colour_palette, 'IL2RA_low')

}
########################################################################################################################
#plot figures
#get data
data_dict={}
for gene in genes:
    df=pd.read_csv(f"{data_directory}{gene}{files_name}")
    df=df[['filtered_score', 'p_val']]
    df['gene']=gene
    data_dict[gene]=df
#create df - gene, p_val, filtered_score
merged_data=pd.concat(list(data_dict.values()), axis=0).reset_index(drop=True)

def plot_strip_all_genes(df, width, height):
    df=df.copy()
    df['color']=np.where(
        (df["filtered_score"] < 0) & (df["p_val"] < 0.01),
        'low', 'neutral'
    )
    df['color'] = np.where(
        (df["filtered_score"] > 0) & (df["p_val"] < 0.01),
        'high', df['color']
    )
    df['gene_specific_color']=df['gene']+'_'+df['color']
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))

    sns.stripplot(data=df[df['color']!='neutral'], x='filtered_score', y='gene', hue='gene_specific_color',
                  palette=gene_specific_palette, s=2,
                  zorder=1, jitter=0.3)
    sns.stripplot(data=df[df['color'] == 'neutral'], x='filtered_score', y='gene', color='lightgrey', s=2,
                  zorder=1, jitter=0.3)
    sns.violinplot(data=df, y='gene', x='filtered_score', fill=False,
                   width=0.8, ax=ax,zorder=2, linewidth=0.8, color='grey', inner='quart')
    ax.set_xlabel(f"Expression Score")
    ax.set_ylabel(f"")
    plt.legend().remove()

    # normalise the tick distance in both axes:
    tick_interval = 2
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))

    # Remove top and right borders
    sns.despine()
    ax.set_yticklabels(
        [r"$\mathit{" + t.get_text() + "}$" for t in ax.get_yticklabels()]
    )
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}_strip.pdf", format='pdf')
    plt.show()
    plt.close()

plot_strip_all_genes(merged_data, fig_width_mm, fig_height_mm)