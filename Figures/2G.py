import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib as mpl
import numpy as np

cwd=os.getcwd()
figure_name='2G'
genes=['OTUD7B']
data_directory=(f"{cwd}/data/6_mer/jurkat/")
files_name='_scores_variant_features.csv'
fig_width_mm=34.25
fig_height_mm=60
splice_AI_filter=0.1

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

########################################################################################################################
for gene in genes:
    df=pd.read_csv(f"{data_directory}{gene}{files_name}")

    #filter spliceai scores over 0.1
    df['filter_splice']=np.where(df['max_score']>splice_AI_filter, False, True)
    df = df[df['filter_splice'] == True]

    #define palette
    color_to_use_high_name=f"{gene}_high"
    color_to_use_high = open_colour_palette(colour_palette, color_to_use_high_name)

    color_to_use_low_name= f"{gene}_low"
    color_to_use_low = open_colour_palette(colour_palette, color_to_use_low_name)

    color_to_use_medium_name=f"{gene}"
    color_to_use_medium = open_colour_palette(colour_palette, color_to_use_medium_name)

    palette={'strong': color_to_use_high, 'moderate': color_to_use_medium, 'weak': color_to_use_low}
    def ATG_violins_kozak(df, width, height, gene):
        df = df.copy()
        df["uORF"] = df["uORF"].fillna("NaN")
        if gene=='OTUD7B':
            fig, ax = plt.subplots(figsize=mm_to_inches((width), height*0.95))
        else:
            fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
        df.sort_values(by='uORF', ascending=True, inplace=True)

        hue_order = ['strong', 'moderate', 'weak']
        white_palette = {k: "white" for k in hue_order}
        palette_background = {'NaN': 'lightgrey', 'oORF_inframe': 'white', 'oORF_outframe': 'white', 'uORF': 'white'}
        sns.violinplot(data=df, x='uORF', y='filtered_score', hue='uORF', palette=palette_background, inner=None,
                       linewidth=0)
        sns.violinplot(data=df, x='uORF', y='filtered_score', hue='kozak_strength',
                       hue_order=hue_order, palette=white_palette, inner=None, linewidth=0.2, linecolor='dimgrey',
                       legend=False)
        sns.stripplot(data=df, x='uORF', y='filtered_score', hue='kozak_strength',
                       hue_order=hue_order, palette=palette, s=2, dodge=True)
        # plot median of each group as horizontal line
        medians = df.groupby('uORF')['filtered_score'].median()
        category_positions = {category: pos for pos, category in enumerate(df['uORF'].unique())}
        for category, median in medians.items():
            if category in category_positions:
                if category != 'NaN':
                    ax.hlines(y=median, xmin=category_positions[category] - 0.5,
                              xmax=category_positions[category] + 0.5, color='dimgrey', linestyle='dashed', linewidth=0.7)

        sns.despine(top=True, right=True)
        plt.xlabel('')

        #rename x ticks for plotting
        tick_rename_map = {
            'NaN': 'None',
            'oORF_inframe': 'oORF in',
            'oORF_outframe': 'oORF out',
            'uORF':'uORF'
        }
        xticks = ax.get_xticks()
        xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
        new_labels = [tick_rename_map.get(label, label) for label in xticklabels]

        # Set labels with right alignment
        ax.set_xticklabels(new_labels, rotation=45, ha='right')

        plt.legend().set_visible(False)
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        fig.legend(
            by_label.values(),
            by_label.keys(),
            title="Kozak strength",
            frameon=False,
            loc="lower center",
            ncol=3,  # horizontal legend
            columnspacing=0.5,  # space between columns
            handletextpad=0.3,  # ← reduce space between marker and text
            handlelength=1.2,  # ← control the length of legend handles
            labelspacing=0.4  # ← reduce vertical spacing between labels
        )
        # Reserve space at the bottom for the legend
        plt.title(r"$\mathit{" + gene + "}$")
        plt.ylabel('Expression Score')
        plt.tight_layout(rect=[0, 0.12, 1, 1])
        # plt.savefig(f"{figures_path}{figure_name}_{gene}.pdf", format='pdf', dpi=300)
        plt.show()
        plt.close()

    ATG_violins_kozak(df, fig_width_mm, fig_height_mm, gene)