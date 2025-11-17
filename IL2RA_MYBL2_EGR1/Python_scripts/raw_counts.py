import sys
import os
import pandas as pd
from functools import reduce
from itertools import product
import pickle

# ########################################################################################################################
#for use in NEMO
path=os.getcwd()
alignment_filter=sys.argv[1]
loci_context_file=sys.argv[2]
insertion_sequences_file=sys.argv[3]
experiment_folder=sys.argv[4]
experiment=sys.argv[5] #Use ME if MYBL2+EGR1 insertions, else use "EGR" or "MYB"

def obtain_files_and_replicates(input_path):
    files={}
    for filename in os.listdir(input_path):
        if 'R1' in filename and experiment in filename:
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
def obtain_genomic_context_multiple_insertions(file):
    genomic_context=pd.read_csv(file)
    dictionary=genomic_context.iloc[0].to_dict()
    return dictionary
def generate_GREP_variants_multiple_insertions(context_dict, insertions_file):
    #obtain insertion sequence info
    insertions_df=pd.read_csv(f"{path}/{insertions_file}")
    insertions_dict=insertions_df.iloc[0].to_dict()
    number_of_inserts=len(context_dict)/2

    def create_greppable(insert_dict):
        dictionary = {}

        # Always include the wild-type sequence
        dictionary["WT"] = (
                context_dict["pre_insert_1"] + context_dict["inter_1_2"] +
                context_dict["inter_2_3"] + context_dict["post_3"]
        )

        insertion_sites = ["inter_1_2", "inter_2_3", "post_3"]  # Valid insert positions
        site_numbers = {"inter_1_2": 1, "inter_2_3": 2, "post_3": 3}  # Map positions to numbers
        insert_keys = list(insert_dict.keys())  # List of insert types (e.g., ['A', 'B'])

        # Generate all possible insert assignments (each site can have None, 'A', or 'B')
        for insert_assignment in product([None] + insert_keys, repeat=len(insertion_sites)):
            if all(i is None for i in insert_assignment):  # Skip WT (already included)
                continue

            mod_seq = [context_dict["pre_insert_1"]]  # Always starts with pre_insert_1
            insert_positions = {}  # Track where we put inserts

            for site, insert in zip(insertion_sites, insert_assignment):
                if insert is not None:
                    mod_seq.append(insert_dict[insert])  # Add the insert sequence
                    insert_positions[site] = insert  # Store assigned insert
                mod_seq.append(context_dict[site])  # Always add the original region

            # Create a dictionary key based on inserted positions
            key_parts = [f"{site_numbers[site]}{insert}" for site, insert in insert_positions.items()]
            dict_key = "_".join(key_parts)

            dictionary[dict_key] = "".join(mod_seq)

        return dictionary
    def add_combinations(dictionary):
        dictionary=dictionary.copy()
        #todo: figure this out.
        return dictionary

    sequences_to_grep={}
    if experiment=='MYB':
        insert={'MYB':insertions_dict['MYB']}
        sequences_to_grep=create_greppable(insert)
    elif experiment=='EGR':
        insert={'EGR':insertions_dict['EGR']}
        sequences_to_grep = create_greppable(insert)
    elif experiment=='ME':
        insert={'EGR':insertions_dict['EGR'], 'MYB':insertions_dict['MYB']}
        sequences_to_grep = create_greppable(insert)
        sequences_to_grep=add_combinations(sequences_to_grep)

    return sequences_to_grep

gene='IL2RA'
#define gene path
gene_path=f"{path}/{experiment_folder}/aligned_reads/{gene}_aligned/{alignment_filter}_aligned_sequences/"

#obtain files and name them appropriately, obtain number of replicates
files, gDNA_files, cDNA_files, neg_files, number_of_replicates=obtain_files_and_replicates(gene_path)

#obtain locus genomic context information
loci_dictionary=obtain_genomic_context_multiple_insertions(f"{path}/{loci_context_file}")

#generate variants to detect
variants_for_detection=generate_GREP_variants_multiple_insertions(loci_dictionary, insertion_sequences_file)

#######################################################################################################################
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
output_path=f"{path}/{experiment_folder}/{experiment}/raw_counts/"
os.makedirs(output_path, exist_ok=True)

raw_counts.to_csv(f"{output_path}{gene}_raw_counts.csv")
sample_read_counts.to_csv(f"{output_path}{gene}_read_counts.csv")

with open(f"{output_path}{gene}_raw_counts.pkl", 'wb') as fp:
    pickle.dump(raw_counts, fp)

with open(f"{output_path}{gene}_read_counts.pkl", 'wb') as fp:
    pickle.dump(sample_read_counts, fp)
