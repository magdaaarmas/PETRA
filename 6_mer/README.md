# Scoring 6-mer insertion variants

Folder contains the necessary codes for analysis of 6-mer insertion libraries. Each step consists of a bash script with dependent Python scripts (included in the `Python_scripts/` directory). 
Examples for all input information files can be found in the `Example_files/` directory

## 1) Fastq to PETRA expression scores 
Bash script fastq.sh is used to score 6-mer insertion variants from fastq files from single-end Ilumina reads (using dependent python scripts contained in folder). The code can be run once to score different experiments. 
Should be run from a directory containing all bash and python codes and containing a subdirectory `{experiment_name}/` containing the necessary files (examples of which can be found in the `Example_files/` directory).

Code line: 
```bash
sbatch fastq.sh experiment_name gene_list
```

experiment_name: name of folder containing the following experiment information data (examples of which can be found in the `Example_files/` directory):	
file_locations.csv - contains gene names and path to folder with fastq files (single-end Illumina reads). Fastq files for all samples from a single experiment (gDNA, cDNA, replicates, negatives) should be stored in a single folder. File names need to include replicate number and sample type in the format "{rep_number}{sample_type}". For example, gDNA from replicate 1 = 1gDNA; cDNA from replicate 3 = 3cDNA. Files from negative controls should be names "neg_cDNA" or "neg_gDNA". 

filtering_parameters.csv - contains gene name, minimim gDNA frequency, condition for replicate removal, number of replicate to remove, maximum spliceAI score for downstream analysis, SpliceAI score category ("max_score" for maximum). 

loci_genomic_context.csv - contains gene name and sequence context. presequence: 10-bp upstream from insertion site; postsequence: 10-bp downstream from insertion site. Sequences in R1 direction (CDS direction). 

alignment_reference_files (folder): folder containing .fasta files of reference amplicons for alignment. Files for all genes should be included in this folder. Each filename should include sequence type "gDNA" or "cDNA" and gene name. 

gene_list: list of genes to score in the format "gene1,gene2,gene3,gene4" (no spaces).  

**Outputs:** 
For each experiment (gene name), the code outputs: 
- aligned_reads directory: contains alignment information, .sam files and .fastq and cigar strings for reads with alignment score > 300. 
- raw_counts.csv - file with number of reads matching each variant
- unfiltered_scores.csv - file with PETRA expression scores prior to gDNA frequency filtering. 
- filtered_scores.csv - file with final PETRA expression scores after gDNA filtering under the "filtered_score" column, and variant p values (normal distribution test).  

## 2) Explore variant features
Bash script explore_variant_features.sh is used to assign variants with a RNA splicing score (previously computed by SpliceAI) and with a ORF category and Kozak sequence strength(using dependent python scripts contained in folder). Taked pre-computed data for TFBS scanning from FIMO to assign variants to TF binding, and compares effects of different TFBS. The code can be run once to score different experiments. 
Outputs variant expression score file containing SpliceAI score and uORF classification, plus statistical data for TFBS analysis. 
Should be run after running the fastq.sh script, in the same directory. 

Code line: 
```bash
  sbatch explore_variant_features.sh experiment_name splice_ai_data_folder path_to_FIMO_files genes_jurkat_ls gene_list
```

Arguments:
experiment_name - same as above (see 1). Should additionally include (examples in the `Example_files/` directory): 
coding_strand_info.csv - contains gene name and coding strand information ("positive" for VAV1 & CD28; "negative" for IL2RA & OTDU7B). 

5UTR_information.csv - contains gene name, the full 5'UTR sequence, the insert position (relative to the 5'TSS), and the 2bp flanking the insertion on either side.  

splice_ai_data_folder - folder containing spliceAI scores. Files should be named "{gene}_output_{splice_distance}.vcf"

path_to_FIMO_files: path to folder containing FIMO output file "fimo.tsv"

genes_jurkat_ls: list of genes with TPM+1 >0 in Jurkat cells

gene_list: list of genes to score in the format "gene1,gene2,gene3,gene4" (no spaces).  

**Outputs:** 
For each experiment (gene name), the code outputs:  
- A `explore_variant_features/` directory containing a {gene}_scores_variant_features.csv file with spliceAI scores and ORF classification for each variant. 
- A `TFBS_analysis/` subdirectory containing a {gene}_TFBS_analysis.csv file with a TF score (log2TF) and q-value (2-sample KS test + Bh correction - BH_q) for each considered TF motif. 
For example of these output files, refer to `Figures/data/`.