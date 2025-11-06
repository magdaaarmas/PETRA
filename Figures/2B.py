import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib as mpl
from functools import reduce

cwd=os.getcwd()
figure_name='2B'
genes=['VAV1', 'CD28', 'IL2RA', 'OTUD7B']
data_directory=(f"{cwd}/data/6_mer/jurkat/")
files_name='_scores_variant_features.csv'
fig_width_mm=90
fig_height_mm=60
variants_selected=50


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


gene_specific_palette = {
    'VAV1':{'low':open_colour_palette(colour_palette, 'VAV1_low'),
           'high':open_colour_palette(colour_palette, 'VAV1_high')},
    'CD28':{'low': open_colour_palette(colour_palette, 'CD28_low'),
           'high': open_colour_palette(colour_palette, 'CD28_high')},
    'OTUD7B':{'low': open_colour_palette(colour_palette, 'OTUD7B_low'),
           'high': open_colour_palette(colour_palette, 'OTUD7B_high')},
    'IL2RA':{'high': open_colour_palette(colour_palette, 'IL2RA_high'),
           'low': open_colour_palette(colour_palette, 'IL2RA_low')}
}

genes_df_strip={}
to_merge_dict={}
for gene in genes:
    df=pd.read_csv(f"{data_directory}{gene}{files_name}")
    df = df[['sequence', 'filtered_score']]
    # mark the top and bottom hits
    df[f"category_{gene}"] = 'NaN'
    df.loc[df['filtered_score'].nlargest(variants_selected).index, f"category_{gene}"] = f"{gene}_high"
    df.loc[df['filtered_score'].nsmallest(variants_selected).index, f"category_{gene}"] = f"{gene}_low"
    df['gene'] = gene
    genes_df_strip[gene] = df
    df = df[['sequence', f"category_{gene}"]]
    to_merge_dict[gene] = df

categories_df_general = reduce(lambda df1, df2: pd.merge(df1, df2, on='sequence', how='outer'), to_merge_dict.values())

reversed_violin_dict={}
for gene in genes:
    df=genes_df_strip[gene]
    df=df[['sequence', 'filtered_score']]
    df=pd.merge(df, categories_df_general, on='sequence', how='left')
    reversed_violin_dict[gene]=df


def plot_violin_per_gene_reversed_modified(dictionary, gene, second_gene, ax):
    dictionary=dictionary.copy()
    #select the scores for gene
    df=dictionary[gene]
    df.rename(columns={'filtered_score':f"score_{gene}"}, inplace=True)
    #only keep categories for gene2
    df=df[[f"score_{gene}", f"category_{second_gene}"]]
    df['category_gene']=second_gene
    #remove the variants that are not top or bottom in gene2
    df=df[df[f"category_{second_gene}"]!='NaN']
    df=df.dropna()
    #categorise based on top or bottom
    df['hit_type']=df[f"category_{second_gene}"].apply(lambda x:x.split('_')[1])

    #plot
    sns.violinplot(
        data=df,
        x='category_gene',
        y=f"score_{gene}",
        hue=f"hit_type",
        hue_order=['high', 'low'],
        split=True,
        gap=.1,
        linewidth=0.2,
        inner=None,
        palette=gene_specific_palette[gene],
        ax=ax  # This makes the plot go to the correct subplot
    )

    # Customize each subplot
    ax.set_title(f"Expression\nscore\nin " + r"$\mathit{" + gene + "}$")
    ax.set_xlabel('')
    ax.set_ylabel('Expression score')
    ax.get_legend().remove()  # Remove legend to avoid duplication


#plot violin plot reversed
fig, axes = plt.subplots(1, 4, figsize=mm_to_inches(fig_width_mm, fig_height_mm), sharey=True)
for i, gene in enumerate(genes):
    plot_violin_per_gene_reversed_modified(reversed_violin_dict, gene, gene, ax=axes[i])
    for gene2 in genes:
        if gene2 != gene:
            plot_violin_per_gene_reversed_modified(reversed_violin_dict, gene, gene2, ax=axes[i])
            if i!=0:
                axes[i].tick_params(axis='y', which='both', left=False, right=False)
for i, ax in enumerate(axes):
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if i!=0:
        ax.spines['left'].set_visible(False)
# Adjust layout and show the plot
plt.tight_layout()
# plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf', dpi=300)
plt.show()
plt.close()