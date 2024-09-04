This directory contains the snakemake files used to produce the annotations used in this manuscript.

The analysis was run using LALB_v3.1.smk. 

We have also provided two python scripts which are used by the snakemake files:
process_braker.py: processes the BRAKER3 results to produce a different naming schema.
softmask.py: creates a softmasked fasta using both the original and hard-masked assembly files

If you wish to use this pipeline on your data, please consider visiting: https://github.com/kocherlab/pipemake