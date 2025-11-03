#!/usr/bin/env bash
#SBATCH --job-name=PETRA_scoring
#SBATCH --partition=ncpu
#SBATCH --time='20:00:00'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100GB
#SBATCH -o PETRA_fastq_to_scores.log
#SBATCH -e PETRA_fastq_to_scores.stderr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=magdalena.armas@crick.ac.uk
set -e  # stop on error

experiment_name=$1
min_alignment_score=300 # minimum alignment score needed to consider read
fastq_path_info="${experiment_name}/file_locations.csv" # contains gene names and path to folder with fastq files
genomic_context_file="${experiment_name}/loci_genomic_context.csv" # gene names + context sequences
filtering_information_file="${experiment_name}/filtering_parameters.csv" # file with filtering parameters
alignment_files_folder="${experiment_name}/alignment_reference_files/" #folder containing .fasta files with amplicon sequences for alignment

# should be run in directory containing "alignment_reference_files" with .fasta files
# align to reference amplicon
#python make_needleall_bash.py "$fastq_path_info" "$alignment_files_folder" "$experiment_name"
#bash run_needle.sh

# filter based on alignment score
#python alignment_filtering.py "$fastq_path_info" "$min_alignment_score" "$experiment_name"


IFS=',' read -r -a genes <<< "$2"
for gene in "${genes[@]}"; do
  echo "Processing gene: $gene"

	# obtain variant counts
	python raw_counts.py "$gene" "$min_alignment_score" "$genomic_context_file" "$experiment_name"

	# score all variants to decide filtering threshold
	python counts_to_scores.py "$gene" "$min_alignment_score" "$filtering_information_file" "$experiment_name"

	# filter scores
	python filter_scores.py "$gene" "$filtering_information_file" "$experiment_name"

	#statistics
	python stats_pnorm.py "$gene" "$experiment_name"
done


