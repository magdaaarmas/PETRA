import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
import os

cwd=os.getcwd()
figure_name='3C'
gene='IL2RA'
data_directory=(f"{cwd}/data/TF_tandem/{gene}_jurkat_scores.csv")
fig_width_mm=57
fig_height_mm=100

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
def obtain_TF_name(string):
    TF_count = string.count('_MA')

    if TF_count == 0:
        TF_name = ''
    elif TF_count == 1:
        index = string.find("_MA")
        substring1 = string[index:]
        first_dot = substring1.find('.')
        second_dot = substring1.find('.', first_dot + 1)
        TF_name_extended = substring1[second_dot + 1:]
        TF_name = TF_name_extended.split('_')[0]
    elif TF_count == 2 and 'sequence' in string:
        index = string.find("_MA")
        substring1 = string[index:]
        first_dot = substring1.find('.')
        second_dot = substring1.find('.', first_dot + 1)
        first_TF_name_extended = substring1[second_dot + 1:]
        first_TF_name = first_TF_name_extended.split('_')[0]

        third_dot = substring1.find('.', second_dot + 1)
        fourth_dot = substring1.find('.', third_dot + 1)
        second_TF_name_extended = substring1[fourth_dot + 1:]
        second_TF_name = second_TF_name_extended.split('_')[0]

        if first_TF_name == second_TF_name:
            TF_name = first_TF_name
        elif first_TF_name != second_TF_name:
            TF_name = f"{first_TF_name} & {second_TF_name}"
    else:
        TF_name = ''
    return TF_name
def obtain_TF_sequence_type(string):
        TF_sequence_type = string.split("_")[0]
        return TF_sequence_type
def plot_scores(df, width, height):
    df = df.copy()

    #define TF names
    df['TF'] = df['origin'].apply(lambda x: obtain_TF_name(x))

    #duplicate rows for sequences that bind 2 TFs:
    df['TF'] = df['TF'].str.split(r'\s*&\s*')  # split on "&" and remove surrounding spaces
    df = df.explode('TF', ignore_index=True)

    #define TF motif insertion type
    df['TF_sequence_type'] = np.where(df['TF'] != '',
                                      df['origin'].apply(lambda x: obtain_TF_sequence_type(x)), '')
    df['TF_detail'] = np.where(df['TF_sequence_type'] == 'sequence',
                                      df['TF_sequence_type']+'_'+df['ID'], df['TF_sequence_type'])

    # make plot
    df = df[df['TF'] != '']
    custom_order=['sequence', 'motif', 'homodimer', 'homotrimer']
    custom_tf_order = ['MYBL2', 'RELA', 'SRF', 'TFAP2A', 'ZNF707', 'TEAD1', 'CTCF']
    df_long = pd.melt(
        df,
        id_vars=['TF', 'TF_sequence_type', 'TF_detail'],
        value_vars=[f"score_{i}_log2_norm" for i in range(1, 4)],
        value_name='score',
        var_name='replicate'
    )

    # Set TF_sequence_type order
    df_long['TF_sequence_type'] = pd.Categorical(
        df_long['TF_sequence_type'],
        categories=custom_order,
        ordered=True
    )

    # Set TF order
    df_long['TF'] = pd.Categorical(
        df_long['TF'],
        categories=custom_tf_order,
        ordered=True
    )

    # Find unique TF_sequence_type values from df
    unique_types = df['TF_detail'].unique()

    # Now sort by both in your preferred order
    df_long = df_long.sort_values(['TF', 'TF_sequence_type'])
    TFs = df_long['TF'].dropna().unique()
    TFs=[Tf for Tf in TFs if '&' not in Tf]

    fig, ax = plt.subplots(nrows=len(TFs), ncols=1,figsize=mm_to_inches(width, height),sharex=True)
    # Reduce horizontal spacing between subplots

    for i, TF in enumerate(TFs):
        if TF in ['MYBL2', 'RELA', 'SRF']:
            palette = {
                'sequence': open_colour_palette(colour_palette, 'IL2RA_sequence'),
                'motif': open_colour_palette(colour_palette, 'IL2RA_motif'),
                'homodimer': open_colour_palette(colour_palette, 'IL2RA_homodimer'),
                'homotrimer': open_colour_palette(colour_palette, 'IL2RA_homotrimer')
            }
        elif TF in ['ZNF707', 'TFAP2A']:
            palette = {
                'sequence': open_colour_palette(colour_palette, 'VAV1_sequence'),
                'motif': open_colour_palette(colour_palette, 'VAV1_motif'),
                'homodimer': open_colour_palette(colour_palette, 'VAV1_homodimer'),
                'homotrimer': open_colour_palette(colour_palette, 'VAV1_homotrimer')
            }
        else:
            palette = {
                'sequence': 'mediumpurple',
                'motif': 'mediumslateblue',
                'homodimer': 'slateblue',
                'homotrimer': 'darkslateblue'
            }
        # Extend palette for any "sequence*" value missing from it
        for val in unique_types:
            if val.startswith('sequence') and val not in palette:
                palette[val] = palette['sequence']

        sns.barplot(
            data = df_long[df_long['TF'].str.contains(TF, na=False)],  # plot only this TF's data
            y='TF_detail',
            x='score',
            hue='TF_detail',
            palette=palette,
            errorbar='se',
            err_kws={"linewidth": 0.5, "color": "dimgrey"},
            ax=ax[i],
            legend=None,
        )
        # Despine this subplot only
        sns.despine(top=True, right=True, ax=ax[i])

        # No individual x-axis labels
        ax[i].set_xlabel("")

        # Set y-axis label to TF name for each subplot
        ax[i].set_ylabel(TF)

        # Get original x-tick labels from Seaborn
        yticklabels = [t.get_text() for t in ax[i].get_yticklabels()]

        # Transform them according to your rules
        new_labels = []
        for label in yticklabels:
            if "sequence" in label:
                # Keep only the part after the underscore
                if "_" in label:
                    new_labels.append(label.split("_", 1)[1])
                else:
                    new_labels.append(label)
            elif label == "motif":
                new_labels.append("motif")
            elif label == "homodimer":
                new_labels.append("motif x2")
            elif label == "homotrimer":
                new_labels.append("motif x3")
            else:
                # Keep unexpected labels unchanged
                new_labels.append(label)

        # Apply the custom x-tick labels
        ax[i].set_yticklabels(new_labels, rotation=0)

        # Hide x-axis completely for all but the last subplot
        if i != len(TFs) - 1:
            ax[i].tick_params(axis='x', bottom=False, labelbottom=False)  # no ticks, no labels
            ax[i].spines['bottom'].set_visible(False)  # hide x-axis line
            ax[i].set_xlabel("")  # no x-label
        else:
            ax[i].set_xlabel("Expression Score")  # add x-axis label only on the last subplot
            ax[i].spines['bottom'].set_visible(True)  # make sure last subplot keeps x-axis line
            ax[i].tick_params(axis='x', bottom=True, labelbottom=True)  # show ticks + labels

    # Rotate all x-tick labels for readability
    for axis in ax:
        for label in axis.get_xticklabels():
            label.set_rotation(0)

    # plt.tight_layout(rect=[0.05, 0, 1, 1])  # leave room for global y-label
    # plt.tight_layout()  # leave room for global y-label
    fig.subplots_adjust(
        left=0.30,  # room for global y-label
        right=0.995,
        top=0.99,  # room for figure title
        bottom=0.08,  # room for rotated x-ticks
        hspace=0.05# gaps between panels

    )
    # plt.savefig(f"{figures_path}.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()

#obtain information for plotting from master file
greys=open_colour_palette(colour_palette, 'greys')
font_to_use, font_size=obtain_font(colour_palette)
color_to_use_name = gene
color_to_use = open_colour_palette(colour_palette, color_to_use_name)

#set plotting parameters
mpl.rcParams['font.size'] = font_size
mpl.rcParams['font.family'] = font_to_use  # or 'Arial', 'Times New Roman', etc.
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

########################################################################################################################
data = pd.read_csv(data_directory)
data.drop(columns=['Unnamed: 0'], inplace=True)

plot_scores(data, fig_width_mm, fig_height_mm)