# Scoring non-coding variants with PETRA

Code for scoring non-coding variants with PETRA.

## 1) 6-mer library:
Folder contains the necessary codes for analysis of 6-mer insertion libraries. Each step consists of a bash script with dependent Python scripts. 

### a) Fastq to expression scores 
Bash script fastq.sh is used to score 6-mer insertion variants from fastq files from single-end Ilumina reads (using dependent python scripts contained in folder). The code can be run once to score different experiments. 

Code line: 
```bash
sbatch fastq.sh experiment_name gene_list
```

experiment_name: name of folder containing the following experiment information data:	
file_locations.csv - contains gene names and path to folder with fastq files (single-end Illumina reads). Fastq files for all samples from a single experiment (gDNA, cDNA, replicates, negatives) should be stored in a single folder. File names need to include replicate number and sample type in the format "{rep_number}{sample_type}". For example, gDNA from replicate 1 = 1gDNA; cDNA from replicate 3 = 3cDNA. Files from negative controls should be names "neg_cDNA" or "neg_gDNA". 

filtering_parameters.csv - contains gene name, minimim gDNA frequency, condition for replicate removal, number of replicate to remove, maximum spliceAI score for downstream analysis, SpliceAI score category ("max_score" for maximum). 

loci_genomic_context.csv - contains gene name and sequence context. presequence: 10-bp upstream from insertion site; postsequence: 10-bp downstream from insertion site. Sequences in R1 direction (CDS direction). 

alignment_reference_files (folder): folder containing .fasta files of reference amplicons for alignment. Files for all genes should be included in this folder. Each filename should include sequence type "gDNA" or "cDNA" and gene name. 

gene_list: list of genes to score in the format "gene1,gene2,gene3,gene4" (no spaces).  

### b) Explore variant features
Code line: 
```bash
  sbatch explore_variant_features.sh experiment_name splice_ai_data_folder path_to_FIMO_files genes_jurkat_ls gene_list
```

Arguments:
experiment_name - same as above (see 1a). Should additionally include: 
coding_strand_info.csv - contains gene name and coding strand information ("positive" for VAV1 & CD28; "negative" for IL2RA & OTDU7B). 

5UTR_information.csv - contains gene name, the full 5'UTR sequence, the insert position (relative to the 5'TSS), and the 2bp flanking the insertion on either side.  

splice_ai_data_folder - folder containing spliceAI scores. Files should be named "{gene}_output_{splice_distance}.vcf"

path_to_FIMO_files: path to folder containing FIMO output file "fimo.tsv"

genes_jurkat_ls: list of genes with TPM+1 >0 in Jurkat cells

gene_list: list of genes to score in the format "gene1,gene2,gene3,gene4" (no spaces).  


