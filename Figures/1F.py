import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import seaborn as sns
import os
import matplotlib as mpl

cwd=os.getcwd()
figure_name='1F'
data_source=(f"{cwd}/data/6_mer/jurkat/VAV1_scores_variant_features.csv")
color_to_use_name='VAV1'
fig_width_mm=50
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
#get data
data=pd.read_csv(data_source)

def plot_strip(df, width, height):
    df=df.copy()
    df['color']=np.where(df['p_val']<0.01, 'p < 0.01', 'p > 0.01')
    df['gene']='VAV1'

    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))

    #plot
    palette={'p < 0.01':color_to_use, 'p > 0.01':'lightgrey'}
    sns.stripplot(data=df, x='filtered_score', y='gene', hue='color', palette=palette, s=3, zorder=1, jitter=0.3)
    sns.violinplot(data=df, y='gene', x='filtered_score', fill=False,
                   width=0.6, ax=ax,zorder=2, linewidth=1, color='grey', inner='quart')
    ax.set_xlabel(f"Expression Score")
    ax.set_ylabel(f"")
    plt.legend().remove()

    # normalise the tick distance in both axes:
    tick_interval = 2
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    plt.yticks([])  # Removes both major and minor y-ticks

    # Remove top and right borders
    sns.despine(top=True, right=True, left=True)
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf')
    plt.show()
    plt.close()

plot_strip(data, fig_width_mm, fig_height_mm)