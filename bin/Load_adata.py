#!/usr/bin/env python3

"""
Usage:
  Load_adata.py --samplesheet=<samplesheet> [options]

Mandatory arguments:
  --samplesheet=<samplesheet>       Samplesheet from nfcore/rnaseq pipeline

Optional arguments:
  --resDir=<resDir>                 Output directory [default: ./]
"""

from docopt import docopt
import pandas as pd

args = docopt(__doc__)
samplesheet = args['--samplesheet']
#path_dir = args['--cellrangerDir']   --cellrangerDir=<cellrangerDir>   Cellranger filtered_feature_bc_matrix.h5 directory
resDir = args['--resDir']

# New run
# Get metadata from samplesheet
meta = pd.read_csv(samplesheet)
meta.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# Reorder by Sample ID, drop double columns, and update index
sample_idx = []
for s in meta['sample'].values:
    sample_idx.append(int(s[2:]))

meta['sample_idx'] = sample_idx
meta = meta.sort_values('sample_idx')
meta.index = meta.sample_idx
meta=meta.drop(columns=['sample_idx'])

meta['sample_idx']=meta['sample']
meta=meta.set_index('sample_idx')
meta.index.name=None

# Drop patient P163
meta = meta[meta['patient'] != 'P163']

# Save metadata as csv
meta.to_csv(f'{resDir}/metadata.csv')
