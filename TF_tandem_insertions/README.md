# Scoring TF motif tandem insertions

Folder contains the necessary codes for analysis of TF motif tandem insertion libraries. Each step consists of a bash script with dependent Python scripts. 
Examples for all input information files can be found in the Example_files directory

## 1) Fastq to expression scores 
Bash script fastq.sh is used to score variants from fastq files from single-end Ilumina reads (using dependent python scripts contained in folder). The code needs to be run once for each experiment.  
Outputs final expression scores and q-value (Welsch's t-test and BH correction against neutral values). 

Code line: 
```bash
sbatch fastq.sh experiment_name gene_name cell_type read
```
**gene_name:** VAV1 or IL2RA

**cell_type:** Jurkat or CD3

**read:** R1 or R2

**experiment_name:** name of folder containing the following experiment information data:	

file_locations.csv - contains gene names and path to folder with fastq files (single-end Illumina reads). Fastq files for all samples from a single experiment (gDNA, cDNA, replicates, negatives) should be stored in a single folder. File names need to contain the gene name and cell type and need to include replicate number and sample type in the format "{experiment_type}_{rep_number}{sample_type}". For example, gDNA from replicate 1 in experiment targeting IL2RA in Jurkat cells = Jurkat_IL2RA_1gDNA. Files from negative controls should be named "neg_cDNA" or "neg_gDNA".

loci_genomic_context.csv - contains gene name and sequence context. presequence: 10-bp upstream from insertion site; postsequence: 10-bp downstream from insertion site. Sequences in R1 direction (CDS direction). 

{gene}_variant_info.csv - contains variant information (insert sequence, variant type, etc - see example files)

alignment_reference_files (folder): folder containing .fasta files of reference amplicons for alignment. Files for all genes should be included in this folder. Each filename should include sequence type "gDNA" or "cDNA" and gene name.

**Outputs:** 
For each experiment a dedicated directory is created as: experiment_name/cell_type/gene, containing: 
- aligned_reads directory: contains alignment information, .sam files and .fastq and cigar strings for reads with alignment score > 300. 
- raw_counts.csv - file with number of reads matching each variant
- scores.csv - file with final PETRA expression scores and statistical parameters.  
  