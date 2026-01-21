import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
import os
import scipy

cwd=os.getcwd()
figure_name='5C'
genes=['VAV1', 'CD28', 'IL2RA', 'OTUD7B']
data_directory=(f"{cwd}/data/6_mer/")
fig_width_mm=80
fig_height_mm=70
min_heatmap=-0.5
max_heatmap=0.7

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
def modify_column_name(text, gene_name, cell_type):
    parts = text.split("_")
    parts = [p for p in parts if p != gene_name and p != cell_type]
    return "_".join(parts)

#obtain information for plotting from master file
font_to_use, font_size=obtain_font(colour_palette)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################
predictions_df={}
for gene in genes:
    gene_predictions={}

    #obtain predictions
    predictions=pd.read_csv(f"{data_directory}jurkat/6_mer_{gene}_PETRA_and_predictions_jurkat.csv")
    predictions.drop(columns=['Unnamed: 0'], inplace=True)
    predictions=predictions[[col for col in predictions.columns if 'filtered_score' not in col]]

    #obtain scores from validation screen
    validation_data=pd.read_csv(f"{cwd}/data/6_mer_validation/{gene}_filtered_scores.csv")
    validation_data.rename(columns={'ID':'sequence', 'filtered_score':f"{gene}_filtered_score"}, inplace=True)

    #select only variants scored in screen
    data=pd.merge(validation_data[['sequence', f"{gene}_filtered_score"]], predictions, on='sequence', how='inner')

    predictions = [prediction for prediction in data.columns if "filtered_score" not in prediction and "sequence" not
                   in prediction]
    for column in predictions:
        stats = scipy.stats.spearmanr(data[f"{gene}_filtered_score"], data[column])

        modified_column_name=modify_column_name(column, gene, 'jurkat')
        gene_predictions[modified_column_name] = stats
    predictions_df[gene]=gene_predictions

def predictions_to_dfs(jurkat_predictions, predictor_order=None):
    rho_dict = {}
    pval_dict = {}

    for gene, predictor_dict in jurkat_predictions.items():
        rho_dict[gene] = {}
        pval_dict[gene] = {}

        for predictor, stat in predictor_dict.items():
            rho_dict[gene][predictor] = stat.correlation
            pval_dict[gene][predictor] = stat.pvalue

    rho_df = pd.DataFrame(rho_dict)
    pval_df = pd.DataFrame(pval_dict)

    if predictor_order is not None:
        rho_df = rho_df.reindex(predictor_order)
        pval_df = pval_df.reindex(predictor_order)

    return rho_df, pval_df
predictor_order=['enformer_CAGE', 'borzoi_CAGE', 'borzoi_RNA',
                 'AG_CAGE_positive','AG_CAGE_negative' , 'AG_RNA_seq']
heatmap_R_df, heatmap_p_df = predictions_to_dfs(predictions_df, predictor_order)

def plot_heatmap(df_correlation, df_pvalues, width, height):
    # Mask for significance
    signif_mask = df_pvalues < 0.05
    # Create figure and axis
    fig, ax = plt.subplots(figsize=mm_to_inches(width, height))
    heatmap = sns.heatmap(
        df_correlation.where(signif_mask),
        cmap=plt.get_cmap('RdYlBu').reversed(),
        mask=~signif_mask,
        cbar=True,
        linewidths=0.5,
        linecolor='lightgray',
        ax=ax,
        vmin=min_heatmap,  # Force symmetric scale
        vmax=max_heatmap,
        center=0  # Ensure 0 is the midpoint
    )
    # set up color bar
    colorbar = heatmap.collections[0].colorbar
    colorbar.set_label("Spearman's ρ")  # or "Correlation", "LogFC", etc.

    # Now overlay gray cells for non-significant values and annotate all cells
    for l in range(df_correlation.shape[0]):
        for m in range(df_correlation.shape[1]):
            value = df_correlation.iloc[l, m]
            is_signif = signif_mask.iloc[l, m]

            # If not significant, color cell gray manually
            if not is_signif:
                ax.add_patch(plt.Rectangle((m, l), 1, 1, color='lightgray', ec='lightgray'))

            # Annotate the cell (regardless of significance)
            ax.text(m + 0.5, l + 0.5, f'{value:.2f}',
                    ha='center', va='center', color='black')

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    # plt.savefig(f"{figures_path}{figure_name}.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()
plot_heatmap(heatmap_R_df, heatmap_p_df, fig_width_mm, fig_height_mm)
