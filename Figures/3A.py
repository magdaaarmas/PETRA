import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator


cwd=os.getcwd()
figure_name='3A'
data_directory=(f"{cwd}/data/IL2RA_MYBL2_EGR1/IL2RA_ME_scores.csv")
fig_width_mm=81
fig_height_mm=81
significance_threshold=0.05
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
def sort_df_ascending(df):
    df=df.copy()
    # df['names'] = df['ID'].apply(lambda x: x.replace(f"_", '\n'))
    df['names'] = df['ID']
    df_sorted = df.sort_values(by='filtered_score', ascending=False)
    return df_sorted
def modify_df_for_plotting(df):
    df=df.copy()
    df_long = df.melt(id_vars=['names'],
                      value_vars=[col for col in df.columns if 'log2_norm' in col and 'average' not in col],
                      var_name='Replicate',
                      value_name='expression_score')
    return df_long
def plot_scores(df,width, height, save_name):
    df=df.copy()
    x_order = df['names'].unique().tolist()

    color_palette={
        'MYB & EGR':combined_color,
        'MYB only':MYB_color,
        'EGR only':EGR_color
    }
    # PLOT
    df['color_code'] = np.where(
        df['names'].str.contains('EGR'),
        np.where(df['names'].str.contains('MYB'), 'MYB & EGR', 'EGR only'),
        'MYB only'
    )

    # colum plot
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
    sns.stripplot(data=df,order=x_order, x='expression_score', y='names', color='black', alpha=0.6, s=2)

    # Plot mean values as bars
    barplot = sns.barplot(x='expression_score', y='names', hue='color_code', palette=color_palette, data=df,
                          errorbar=None, order=x_order,
                          width=1, linewidth=1, edgecolor='white')
    plt.xlabel('Expression score')
    plt.ylabel('')

    ax.legend(frameon=False, loc='lower right')

    sns.despine(top=True, right=True, left=True)
    plt.tight_layout()
    # plt.savefig(f"{figures_path}_{figure_name}.pdf", transparent=True, format='pdf')
    plt.show()
    plt.close()

#obtain information for plotting from master file
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################
data=pd.read_csv(data_directory)
data.drop(columns=['Unnamed: 0'], inplace=True)

MYB_color=open_colour_palette(colour_palette, 'MYBL2')
EGR_color=open_colour_palette(colour_palette, 'EGR')
combined_color=open_colour_palette(colour_palette, 'MYBL2EGR')

data_sorted_ascending=sort_df_ascending(data)
data_sorted_ascending=modify_df_for_plotting(data_sorted_ascending)

plot_scores(data_sorted_ascending, fig_width_mm, fig_height_mm, 'AAA')
