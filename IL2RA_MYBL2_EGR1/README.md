# Scoring MYBL2 and EGR1 insertions in IL2RA

Folder contains the necessary codes for analysis of libraries of multiple MYBL2 and EGR1 insertions. Each step consists of a bash script with dependent Python scripts. 
Examples for all input information files can be found in the Example_files directory

## 1) Fastq to expression scores 
Bash script fastq.sh is used to score variants from fastq files from single-end Ilumina reads (using dependent python scripts contained in folder). The code can be run once to score different experiments. 
Outputs final expression scores.

Code line: 
```bash
sbatch fastq.sh experiment_name experiment_type
```
**experiment_type:** use "ME" if editing with MYBL2 and EGR1 sites simultaneously, else use "MYB" or "EGR1"


**experiment_name:** name of folder containing the following experiment information data:	

file_locations.csv - contains gene names and path to folder with fastq files (single-end Illumina reads). Fastq files for all samples from a single experiment (gDNA, cDNA, replicates, negatives) should be stored in a single folder. File names have to start with the experiment type and need to include replicate number and sample type in the format "{experiment_type}_{rep_number}{sample_type}". For example, gDNA from replicate 1 in MYBL3+EGR1 experiment = ME_1gDNA; cDNA from replicate 3 in only EGR1 experiment= EGR_3cDNA. Files from negative controls should be names "neg_cDNA" or "neg_gDNA".

loci_genomic_context.csv -  contains gene name and sequence context for all segments between inserts plus 10bp up- and downstream from first and last insertion, repectively (see example file) 

alignment_reference_files (folder): folder containing .fasta files of reference amplicons for alignment. Files for all genes should be included in this folder. Each filename should include sequence type "gDNA" or "cDNA" and gene name.

**Outputs:** 
For each experiment (gene name), the code outputs: 
- aligned_reads directory: contains alignment information, .sam files and .fastq and cigar strings for reads with alignment score > 300. 
- raw_counts.csv - file with number of reads matching each variant
- unfiltered_scores.csv - file with PETRA expression scores prior to gDNA frequency filtering. 
- filtered_scores.csv - file with final PETRA expression scores after gDNA filtering under the "filtered_score" column.  