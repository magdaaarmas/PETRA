import sys
import os
import pandas as pd
from functools import reduce
from itertools import product
import pickle

########################################################################################################################
#for use in NEMO
path=os.getcwd()
gene=sys.argv[1]
cell_type=sys.argv[2]
alignment_filter=sys.argv[3]
loci_context_file=sys.argv[4]
variant_information_file=sys.argv[5]
experiment_folder=sys.argv[6]

########################################################################################################################

def obtain_files_and_replicates(input_path):
    files={}
    for filename in os.listdir(input_path):
        index=filename.find('DNA')
        if index != -1 and 'neg' not in filename:
            name=filename[index-2:index+3]
            files[name]=input_path+filename
        elif  index != -1 and 'neg' in filename:
            name = filename[index - 5:index + 3]
            files[name] = input_path + filename
    gDNA_files={key:value for key,value in files.items() if 'gDNA' in key and 'neg' not in key}
    cDNA_files={key:value for key,value in files.items() if 'cDNA' in key and 'neg' not in key}
    neg_files={key:value for key,value in files.items() if 'neg' in key}
    number_of_replicates=len(gDNA_files)
    return files, gDNA_files, cDNA_files, neg_files, number_of_replicates
def analyse_fastq_reads(sequence_file, dictionary, sample):
    # Prepare a dictionary to store counts for each extension
    results = {key: 0 for key, value in dictionary.items()}
    # Counter for exclusions
    nothing_found = 0
    read_number=0
    # Iterate over the reads in the FASTQ file
    with open(sequence_file, 'r') as file:
        for line in file:
            read = line.strip()
            variant_found=False
            read_number+=1
            for key, value in dictionary.items():
                variant_to_search = value
                if variant_to_search in read:
                    results[key] += 1
                    variant_found = True
                    break  # Exit the loop once a match is found
            # Increment no_extension_count only if no extensions were found
            if not variant_found:
                nothing_found += 1
            print(f"done with {read_number}")
    # Create a results DataFrame with counts
    result_df = pd.DataFrame(list(results.items()), columns=['ID', f"{sample}_count"])

    # Add the exclusion count for no extensions found
    exclusion_df = pd.DataFrame({
        'ID': ['no_extension_found'],
        f"{sample}_count": [nothing_found]
    })

    # Combine the two DataFrames
    final_result_df = pd.concat([result_df, exclusion_df], ignore_index=True)

    return final_result_df, read_number
def obtain_genomic_context(file):
    genomic_context=pd.read_csv(file)
    pre_sequence_value=genomic_context.loc[genomic_context['gene']==gene, 'presequence'].values
    pre_sequence=pre_sequence_value[0] if len(pre_sequence_value) > 0 else None
    post_sequence_value=genomic_context.loc[genomic_context['gene']==gene, 'postsequence'].values
    post_sequence=post_sequence_value[0] if len(post_sequence_value) > 0 else None
    return pre_sequence, post_sequence
def generate_GREP_variants_from_file(peg_information_file, pre_extension, post_extension):
    peg_info_df=pd.read_csv(f"{path}/{peg_information_file}")
    inserts=peg_info_df['sequence'].tolist()
    inserts_in_context={}
    for insert in inserts:
        insert_with_context=pre_extension+insert+post_extension
        inserts_in_context[insert]=insert_with_context
    inserts_in_context['WT']=pre_extension+post_extension
    return inserts_in_context

########################################################################################################################
sample_path=f"{path}/{experiment_folder}/{cell_type}/{gene}"
reads_path=f"{sample_path}/aligned_reads/{gene}_aligned/{alignment_filter}_aligned_sequences/"

#obtain files and name them appropriately, obtain number of replicates
files, gDNA_files, cDNA_files, neg_files, number_of_replicates=obtain_files_and_replicates(reads_path)

#obtain locus genomic context information
pre_sequence, post_sequence=obtain_genomic_context(loci_context_file)

variants_for_detection=generate_GREP_variants_from_file(variant_information_file, pre_sequence, post_sequence)

########################################################################################################################
#obtain count matrix
total_reads_dict={}
counts_dict={}
for sample, document in files.items():
    count_df, total_reads=analyse_fastq_reads(document, variants_for_detection, sample)
    total_reads_dict[sample]=total_reads
    counts_dict[sample]=count_df

sample_read_counts = pd.DataFrame(list(total_reads_dict.items()), columns=["sample", "total_reads_for_analysis"])
raw_counts=reduce(lambda left, right:pd.merge(left, right, on='ID', how='outer'), counts_dict.values())
raw_counts.loc['Total'] = raw_counts.sum()

########################################################################################################################
output_path=f"{sample_path}/raw_counts/"
os.makedirs(output_path, exist_ok=True)

raw_counts.to_csv(f"{output_path}raw_counts.csv")

with open(f"{output_path}raw_counts.pkl", 'wb') as fp:
    pickle.dump(raw_counts, fp)





