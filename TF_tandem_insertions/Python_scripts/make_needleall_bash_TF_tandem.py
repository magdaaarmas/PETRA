import os
import sys
import pandas as pd

file_locations=sys.argv[1]
gene=sys.argv[2]
cell_type=sys.argv[3]
read=sys.argv[4]
reference_directory=sys.argv[5]
experiment_name=sys.argv[6]
output_dir=sys.argv[7]

working_dir = os.getcwd()
shell_file = open(f"{output_dir}/run_needle.sh", 'w')
shell_file.write("#!/usr/bin/env bash\n\n#SBATCH --partition=ncpu\n#SBATCH --time='20:00:00'\n#SBATCH --ntasks=15\n#SBATCH --mem=64GB\n#SBATCH -o align.log\n#SBATCH -e align.stderr\n")
shell_file.write('ml purge\nmodule load EMBOSS/6.6.0-foss-2016b\n') #ideally would implement version control here...
shell_file.write(f"mkdir -p aligned_reads/{gene}_aligned/\n")


file_paths=pd.read_csv(file_locations)
path_value = file_paths.loc[file_paths['gene'] == gene, 'path'].values
fastq_directory=path_value[0] if len(path_value) > 0 else None

for filename in os.listdir(reference_directory):
    searchable=filename.upper()
    if gene in searchable and 'GDNA' in searchable:
        reference_gDNA=reference_directory+filename
    elif gene in searchable and 'CDNA' in searchable:
        reference_cDNA=reference_directory+filename

for i in os.listdir(fastq_directory):
    searchable=i.upper()
    if i.endswith(".fastq") and gene in searchable and cell_type in searchable and read in searchable:
        index_of_first_dot = i.find('.')
        sample_name = i[:index_of_first_dot]
        if "gDNA" in i:
            reference = reference_gDNA
        elif "cDNA" in i:
            reference = reference_cDNA
        else:
            continue  # Skip files that don't contain either "gDNA" or "cDNA"
        shell_file.write(f"needleall -asequence {working_dir}/{reference} -bsequence {fastq_directory}/{i} -gapopen 10 -gapextend 0.5 -outfile aligned_reads/{gene}_aligned/{sample_name}.sam -aformat sam &\n")

shell_file.write('wait\necho "finished running all alignments"\n')
shell_file.close()

