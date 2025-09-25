# Example files
The Experiment1 folder contains the necessary files to run the full analysis for 6-mer insertions targeting the 5' UTR of VAV1, IL2RA, OTUD7B and CD28, as published. 

For this, the bash codes need to be run (from the Examples directory) as follows: 

```bash
sbatch fastq.sh Experiment_1 "VAV1,IL2RA,CD28,OTUD7B"

sbatch explore_variant_features.sh Experiment_1 splice_AI_outputs/  FIMO_outputs/ genes_expressed_in_jurkat_ls.pkl "VAV1,IL2RA,CD28,OTUD7B"
```

