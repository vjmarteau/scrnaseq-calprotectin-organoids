#!/bin/bash
# Run scRNAseq organoids calprotectin analysis pipeline

nextflow run ./main.nf -profile cluster \
-w /data/scratch/marteau/nf-work-dir/nf-scRNAseq-organoids-calprotectin/work \
-resume