import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator


cwd=os.getcwd()
figure_name='2I'
genes=['CD28']
data_directory=(f"{cwd}/data/6_mer/jurkat/")
files_name='_TFBS_analysis.csv'
fig_width_mm=45
fig_height_mm=60
significance_threshold=0.05

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
# color_to_use=open_colour_palette(colour_palette, color_to_use_name)
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################


for gene in genes:
    data = pd.read_csv(f"{data_directory}{gene}{files_name}")
    def plot_volcano_unlabelled(df, gene, TF_to_label_ls, width, height):
        df=df.copy()
        # make plot
        # define palette
        color_to_use_high_name = f"{gene}_high"
        color_to_use_high = open_colour_palette(colour_palette, color_to_use_high_name)

        color_to_use_low_name = f"{gene}_low"
        color_to_use_low = open_colour_palette(colour_palette, color_to_use_low_name)

        volcano_palette = {'increase': color_to_use_high, 'decrease':color_to_use_low,
                           'neutral': 'lightgrey'}

        #calculate -log10 q value
        df[f"-log10p_BH"]=df['BH_q'].apply(lambda x: -np.log10(x))

        #assign a significance column
        df['BH_sig']=np.where(df['BH_q'] < significance_threshold, True, False)

        #asssign volcano classification for colouring
        conditions_increase = (df['BH_sig'] == True) & (df['log2FC'] > 0)
        conditions_decrease = (df['BH_sig'] == True) & (df['log2FC'] < 0)
        df['volcano_classification'] = np.where(conditions_increase, 'increase',
                                                np.where(conditions_decrease, 'decrease', 'neutral'))

        #plot
        fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
        ax = sns.scatterplot(data=df, x='log2FC', y=f"-log10p_BH", hue='volcano_classification',
                             edgecolor=None,
                             s=6, palette=volcano_palette)

        # label TFs in IL2RA and VAV1
        if gene in ['IL2RA', 'VAV1']:

            # Start with everything unlabelled
            df['label'] = False

            # For each TF in the list, label only the row with the highest log2FC
            for tf in TF_to_label_ls:
                subset = df[df['TF'].str.contains(tf)]
                if not subset.empty:
                    top_idx = subset['log2FC'].idxmax()
                    df.loc[top_idx, 'label'] = True

            for i, row in df[df['label']].iterrows():
                x = row['log2FC']
                y = row[f"-log10p_BH"]

                # Increase randomness for label placement
                offset_x = -5 + np.random.uniform(-15, 0)  # in points
                offset_y = -5 + np.random.uniform(-15, 0)  # in points

                ax.annotate(
                    str(row['TF']).split('.')[-1],
                    xy=(x, y),
                    xytext=(offset_x, offset_y),
                    textcoords="offset points",
                    fontsize=5,
                    ha='left',
                    va='bottom',
                    color='dimgrey',
                    arrowprops=dict(
                        arrowstyle="-",
                        color='dimgrey',
                        lw=0.4,
                        shrinkA=0,  # start exactly at the point
                        shrinkB=0,  # end exactly at the label
                        connectionstyle="arc3,rad=0"  # ensures straight line connection
                    )
                )

        plt.xlabel('log2 FC')
        plt.ylabel('-log10 q value')
        # ax.set_xlim(-5, 5)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.legend().remove()
        sns.despine(right=True, top=True)
        plt.title(r"$\mathit{" + gene + "}$")
        plt.tight_layout()
        # plt.savefig(f"{figures_path}{figure_name}_{gene}.pdf", format='pdf', dpi=300)
        plt.show()
        plt.close()

    plot_volcano_unlabelled(data, gene, ['MYBL2', 'SRF', 'RELA', 'CTCF', 'TEAD1', 'MA0003.3.TFAP2A', 'ZNF707'], fig_width_mm, fig_height_mm)