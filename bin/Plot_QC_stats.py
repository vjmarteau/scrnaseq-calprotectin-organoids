#!/usr/bin/env python3

"""
Usage:
  Plot_QC_stats.py --adata=<adata> --adata_filtered=<adata_filtered> --adata_nodoublet=<adata_nodoublet> [options]

Mandatory arguments:
  --adata=<adata>                     adata
  --adata_filtered=<adata_filtered>   adata_filtered
  --adata_nodoublet=<adata_nodoublet> adata_nodoublet

Optional arguments:
  --resDir=<resDir>                    Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sc.settings.set_figure_params(figsize=(5,5))
sns.set(font_scale=2)

args = docopt(__doc__)
adata = args["--adata"]
adata_filtered = args["--adata_filtered"]
adata_nodoublet = args["--adata_nodoublet"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)
filtered = sc.read_h5ad(adata_filtered)
nodoublet = sc.read_h5ad(adata_nodoublet)

# Make stats table filter steps n_counts per sample
# Create new label: patient + group
adata.obs['label'] = adata.obs['patient'].astype(str) + '-' + adata.obs['group'].astype(str)
# Get label factor levels as list in right order
label_l = []
for p in adata.obs.patient.cat.categories:
    for g in adata.obs.group.cat.categories:
        label_l.append(f'{p}-{g}')
# Use created list to relevel label factor
adata.obs["label"] = pd.Categorical(adata.obs["label"],
                                     categories=label_l)

stats_adata = pd.DataFrame(adata.obs[['sample']].value_counts(), columns=['adata'])
stats_filtered = pd.DataFrame(filtered.obs[['sample']].value_counts(), columns=['filtered'])
stats_nodoublet = pd.DataFrame(nodoublet.obs[['sample']].value_counts(), columns=['nodoublet'])

stats = pd.concat([stats_adata, stats_filtered, stats_nodoublet], axis = 'columns')
stats = stats.loc[list(adata.obs['sample'].unique())]
stats['label'] = list(adata.obs['label'].unique())
stats = stats.reset_index()
stats.index = stats['sample']
stats = stats.T.iloc[:-1]

# Save stats as csv
stats.to_csv(f"{resDir}/stats.csv")

# Plot Violin Plots and _highest_expr_genes
p_name_l = ['adata','filtered','nodoublet']
for _adata, p_name in zip([adata, filtered, nodoublet], p_name_l):

    plot_adata = _adata.copy()
    adata_p = adata[plot_adata.obs.index, plot_adata.var.index].copy()

    # Append "sample_counts" as column
    # First initialize with NAN and loop to fill values
    adata_p.obs['sample_counts'] = np.NaN
    for s in adata_p.obs["sample"].cat.categories:
        index_w_sel_sample = adata_p.obs.where(adata_p.obs['sample']==s).dropna(how='all').index
        adata_p.obs.loc[index_w_sel_sample,"sample_counts"] = adata_p.obs['sample'].value_counts().loc[s]
    adata_p.obs['sample_counts'] = adata_p.obs['sample_counts'].astype(int)

    # Create new label column for plots: patient-batch + group + sample_counts
    adata_p.obs["label"] = adata_p.obs['patient'].astype(str) + '-B' + adata_p.obs['batch'].astype(str) + '\n' + adata_p.obs['group'].astype(str) + '\nn=' + adata_p.obs['sample_counts'].astype(str)

    # Get label factor levels as list in right order
    label_l = []
    for p in adata_p.obs.patient.cat.categories:
        for g in adata_p.obs.group.cat.categories:
            _sel_adata = adata_p.obs[(adata_p.obs.patient==p)&(adata_p.obs.group==g)]
            c = _sel_adata['sample_counts'].unique()[0]
            b = _sel_adata['batch'].unique()[0]
            label_l.append(f'{p}-B{b}\n{g}\nn={c}')

    # Use created list to relevel label factor
    adata_p.obs["label"] = pd.Categorical(adata_p.obs["label"],
                                          categories=label_l)

    # Calculate QC metrics for all samples
    adata_p.var['mito'] = adata_p.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata_p, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)

    fig, ax = plt.subplots(1,1,figsize=(10,10))
    sc.pl.highest_expr_genes(adata_p, n_top=20, ax=ax)
    fig.savefig(f"{resDir}/{p_name}_highest_expr_genes.png", bbox_inches="tight")

    fig, ax = plt.subplots(1,1,figsize=(30,10))
    sc.pl.violin(adata_p, 'total_counts', groupby='label', size=2, log=True, cut=0, ax=ax)
    fig.savefig(f"{resDir}/{p_name}_Violin_total_counts.png", bbox_inches="tight")

    fig, ax = plt.subplots(1,1,figsize=(30,10))
    sc.pl.violin(adata_p, 'pct_counts_mito', groupby='label', ax=ax)
    fig.savefig(f"{resDir}/{p_name}_Violin_pct_counts_mito.png", bbox_inches="tight")