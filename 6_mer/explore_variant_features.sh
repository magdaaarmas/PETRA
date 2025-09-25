#!/usr/bin/env bash

#SBATCH --job-name=PETRA_explore_variant_features
#SBATCH --partition=ncpu
#SBATCH --time='20:00:00'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100GB
#SBATCH -o PETRA_explore_variant_features.log
#SBATCH -e PETRA_explore_variant_features.stderr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=magdalena.armas@crick.ac.uk

experiment_name=$1
splice_ai_folder=$2
FIMO_path=$3
gene_expression_jurkat_list=$4

coding_strand_file="${experiment_name}/coding_strand_info.csv"
UTR_info_file="${experiment_name}/5UTR_information.csv"

filtered_scores="${experiment_name}/filtered_scores/"
output_folder="${experiment_name}/explore_variant_features/"
mkdir -p "$output_folder"

IFS=',' read -r -a genes <<< "$5"
for raw in "${genes[@]}"; do
  gene="${raw//[[:space:]]/}"              # remove spaces/tabs/newlines
	gene="$(printf '%s' "$gene" | tr -d '[:space:]"'\''')"
  echo "Processing gene: '$gene'"

  #incorporate SpliceAI scores
	python SpliceAI_score_processing.py "$gene" "$splice_ai_folder" "$coding_strand_file" "$output_folder" "$filtered_scores"

	# incorporate uORF effects
	python uORF_detection.py "$gene" "$UTR_info_file" "$output_folder"

	#TFBS analysis
	python TFBS_analysis.py "$gene" "$FIMO_path" "$output_folder" "$gene_expression_jurkat_list"

done

