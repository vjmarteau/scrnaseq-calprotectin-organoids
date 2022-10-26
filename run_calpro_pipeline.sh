#!/bin/bash
# Run scRNAseq organoids calprotectin analysis pipeline

nextflow run ./main.nf -profile cluster \
-resume