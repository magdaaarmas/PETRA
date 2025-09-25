import os
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests
import seaborn as sns
from scipy import stats

########################################################################################################################

genes = ['IL2RA', 'VAV1', 'CD28', 'OTUD7B']
FIMO_threshold='001'
#define the maximum score to be included in the analysis - to filter out the variants with high splice effects
# based on max_score
#if none, use 1
max_score_threshold=0.1
#define the minimum ES to be considered (filter to be used in combination with the max_score_threshold).
#if none, use False
ES_threshold=False
#input path to data processed by FIMO_procees_score_1 code
input_path=f"/Volumes/lab-findlayg/home/users/armasm/PETRA_large_screen_analysis/filtered_scores/FIMO_analysis_{FIMO_threshold}/"

if ES_threshold==False:
    ES_name=''
else:
    ES_name=str(ES_threshold)
if max_score_threshold==1:
    max_score_threshold_name='None'
else:
    max_score_threshold_name=str(max_score_threshold).split('.')[-1]
output_path=f"{input_path}/statistics_SAI_{max_score_threshold_name}_ES_{ES_name}/"
os.makedirs(output_path, exist_ok=True)

########################################################################################################################
def open_pickle(file):
    with open(file, 'rb') as fp:
        df=pickle.load(fp)
    return df

def plot_histogram(list, gene, name):
    list=list.copy()
    fig, ax = plt.subplots(figsize=(5,5))
    ax=sns.histplot(list, bins=50)
    y_min, y_max = ax.get_ylim()
    plt.vlines(x=8, ymin=y_min, ymax=y_max, colors='red')
    plt.title(f"{gene}_{name}")
    sns.despine(top=True, right=True)
    plt.savefig(f"{output_path}{gene}_{name}.png", dpi=300)
    # plt.show()
    plt.close()
def compute_TFBS_statistics(hits_list_dict, gene, score_df):
    dictionary = hits_list_dict.copy()
    # open data
    data = score_df

    # select hits
    gene_hits = dictionary[gene]

    # calculate statistics
    hits_statistics_mixed = []
    hits_statistics_KS = []
    number_negative=[]
    number_positive=[]
    normal_distribution_count=0
    for hit in gene_hits:
        positive = data[data[hit] == True]['filtered_score'].tolist()
        positive_median = np.median(positive)
        negative = data[data[hit] == False]['filtered_score'].tolist()
        negative_median = np.median(negative)
        log2FC = positive_median - negative_median
        if np.abs(positive_median)>np.abs(negative_median):
            closer_to_zero='negative'
        else:
            closer_to_zero = 'positive'
        if log2FC>0:
            if closer_to_zero=='negative':
                change='TF_increase'
            else:
                change = 'TF_zero'
        else:
            if closer_to_zero == 'negative':
                change = 'TF_decrease'
            else:
                change = 'TF_zero'

        #test normality
        if len(positive)>3 and len(negative)>3:
            norm_stat_positive, norm_p_positive=stats.shapiro(positive)
            norm_stat_negative, norm_p_negative = stats.shapiro(negative)
            if norm_p_negative>0.05 and norm_p_positive>0.05:
                normal_distribution_count+=1

        #statistical testing
        if len(positive) > 8 and len(negative) > 8:
            stat, p = mannwhitneyu(positive, negative)
            list = (hit, log2FC, change, stat, p)
            stat2, p2 = ks_2samp(positive, negative)
            list2 = (hit, log2FC, change, f"{stat2}_KS", p2)
        elif len(positive) < 1 or len(negative) < 1:
            list = (hit, 'no_items','', '', '')
            list2 = list
        else:
            stat, p = ks_2samp(positive, negative)
            list = (hit, log2FC, change, f"{stat}_KS", p)
            list2 = list
        hits_statistics_mixed.append(list)
        hits_statistics_KS.append(list2)
        number_positive.extend([len(positive)])
        number_negative.extend([len(negative)])
        print(f"done with {hit}")

    statistics_mixed = pd.DataFrame(hits_statistics_mixed, columns=['TF', 'log2FC', 'change_type','t', 'p'])
    statistics_KS = pd.DataFrame(hits_statistics_KS, columns=['TF', 'log2FC','change_type', 't', 'p'])
    plot_histogram(number_positive, gene,f"hist_positive_TF")
    plot_histogram(number_negative, gene, f"hist_negative_TF")
    data_normality=f"{normal_distribution_count}/{len(gene_hits)}"
    with open(f"{output_path}{gene}_TF_normality.txt", 'w') as fp:
        fp.write(data_normality)
    return statistics_mixed, statistics_KS
def p_value_corrections(df_with_stats, save_name):
    info=df_with_stats
    #filter out errors in statistical processing
    info=info[info['log2FC']!='no_items']
    #apply BH correction
    _, corrected_p_values, _, _ = multipletests(info['p'], method='fdr_bh')
    info['BH_q']=corrected_p_values
    info['BH_sig5']=np.where(info['BH_q']<=0.05, True, False)
    info['BH_sig1'] = np.where(info['BH_q'] <= 0.01, True, False)
    #apply Bonferroni correction
    total_comparisons=len(info)
    Bonferroni_p1=0.01/total_comparisons
    Bonferroni_p5=0.05/total_comparisons
    info['Bonferroni_sig1']=np.where(info['p']<Bonferroni_p1, True, False)
    info['Bonferroni_sig5'] = np.where(info['p'] < Bonferroni_p5, True, False)
    info.to_csv(f"{output_path}{save_name}.csv")
    with open(f"{output_path}{save_name}.pkl", 'wb') as fp:
        pickle.dump(info, fp)

    significant=info[info['BH_sig5']==True]
    significant_top=len(significant[significant['log2FC']>0])
    significant_bottom=len(significant[significant['log2FC']<0])
    return significant_top, significant_bottom

########################################################################################################################
full_df_dict={}
variants_remaining=[]
for gene in genes:
    df=open_pickle(f"{input_path}{gene}_ES_TFBS_merged.pkl")
    if ES_threshold:
        df['filter_out']=np.where((df['filtered_score']<=ES_threshold)&(df['max_score']>=max_score_threshold), True, False)
        df=df[df['filter_out']==False]
    else:
        df['filter_out'] = np.where(df['max_score'] >= max_score_threshold,True, False)
        df = df[df['filter_out'] == False]
    full_df_dict[gene]=df
    variants=[gene, len(df)]
    variants_remaining.append(variants)
variants_remaining=pd.DataFrame(variants_remaining, columns=['gene', 'variants_post_filter'])
variants_remaining.to_csv(f"{output_path}post_filter_variant_numbers.csv")

hits_ls_dict=open_pickle(f"{input_path}FIMO_dict_hits_list.pkl")

stat_mixed_dict={}
stat_KS_dict={}
for gene in genes:
    print(f"starting with {gene}")
    stat_mixed_dict[gene], stat_KS_dict[gene]=compute_TFBS_statistics(hits_ls_dict, gene, full_df_dict[gene])


#find BH and Bonferroni significances
volcano_data={}
hits_numbers=[]
for gene in genes:
    info_KS=stat_KS_dict[gene]
    info_mixed=stat_mixed_dict[gene]
    hits_KS_top, hits_KS_bottom =p_value_corrections(info_KS, f"{gene}_stats_KS")
    hits_mixed_top, hits_mixed_bottom= p_value_corrections(info_mixed, f"{gene}_stats_mixed")
    list=[gene, hits_KS_top, hits_KS_bottom, hits_mixed_top, hits_mixed_bottom]
    hits_numbers.append(list)
hits_numbers=pd.DataFrame(hits_numbers, columns=['gene', 'top_KS', 'bottom_KS', 'top_mixed', 'bottom_mixed'])
hits_numbers.to_csv(f"{output_path}hits_numbers.csv")
