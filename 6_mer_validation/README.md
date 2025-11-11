# Scoring 6-mer validation library

Folder contains the necessary codes for analysis of 6-mer insertion libraries. Each step consists of a bash script with dependent Python scripts. 
Examples for all input information files can be found in the Example_files directory

## 1) Fastq to expression scores 
Bash script fastq.sh is used to score 6-mer insertion variants from fastq files from single-end Ilumina reads (using dependent python scripts contained in folder). The code can be run once to score different experiments. 
Outputs final expression scores.

Code line: 
```bash
sbatch fastq.sh experiment_name gene_list
```

experiment_name: name of folder containing the following experiment information data:	
file_locations.csv - contains gene names and path to folder with fastq files (single-end Illumina reads). Fastq files for all samples from a single experiment (gDNA, cDNA, replicates, negatives) should be stored in a single folder. File names need to include replicate number and sample type in the format "{rep_number}{sample_type}". For example, gDNA from replicate 1 = 1gDNA; cDNA from replicate 3 = 3cDNA. Files from negative controls should be names "neg_cDNA" or "neg_gDNA". 

filtering_parameters.csv - contains gene name, minimim gDNA frequency

loci_genomic_context.csv - contains gene name and sequence context. presequence: 10-bp upstream from insertion site; postsequence: 10-bp downstream from insertion site. Sequences in R1 direction (CDS direction). 

pegRNA_info.csv - contains gene names and path to folder with to file containing pegRNA information (see example_pegRNA_file in Example_files directory). Contains "peg_target_strand" colum specifying "coding" or "non-coding" based on the strand targeted by the pegRNA.  

alignment_reference_files (folder): folder containing .fasta files of reference amplicons for alignment. Files for all genes should be included in this folder. Each filename should include sequence type "gDNA" or "cDNA" and gene name.

gene_list: list of genes to score in the format "gene1,gene2,gene3,gene4" (no spaces).

**Outputs:** 
For each experiment (gene name), the code outputs: 
- aligned_reads directory: contains alignment information, .sam files and .fastq and cigar strings for reads with alignment score > 300. 
- raw_counts.csv - file with number of reads matching each variant
- unfiltered_scores.csv - file with PETRA expression scores prior to gDNA frequency filtering. 
- filtered_scores.csv - file with final PETRA expression scores after gDNA filtering under the "filtered_score" column.  