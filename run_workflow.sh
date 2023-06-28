#!/bin/bash
# Run scRNAseq organoids calprotectin analysis pipeline

nextflow run ./main.nf \
-profile cluster \
-params-file params.yaml \
-w /data/scratch/marteau/nf-work-dir/calpro/work \
-resume