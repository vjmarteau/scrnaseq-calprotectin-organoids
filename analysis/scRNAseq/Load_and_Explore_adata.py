# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python [conda env:CRCA-2022-crca-scanpy]
#     language: python
#     name: conda-env-CRCA-2022-crca-scanpy-py
# ---

# %% [markdown]
# # Calprotectin organoids:

# %% [markdown]
# ## 1. Load scRNAseq data
# Load cellranger output matrices and concatenate all samples to single adata

# %%
import os
import gzip
import anndata as ad
import scanpy as sc
import scipy as sp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sc.settings.set_figure_params(figsize=(5,5))
sns.set(font_scale=2)
fig_dir="../../results/figures/scRNAseq"

# %% [markdown]
# ### Get metadata
# This project had to be repeated due to technical issues with 10x. I will load both runs anyway.

# %%
# New run
# Get metadata from samplesheet
meta = pd.read_csv('../../tables/samplesheet.csv')
meta.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# Reorder by Sample ID, drop double columns, and update index
sample_idx = []
for s in meta['sample'].values:
    sample_idx.append(int(s[2:]))

meta['sample_idx'] = sample_idx
meta = meta.sort_values('sample_idx')
meta.index = meta.sample_idx
meta=meta.drop(columns=['sample_idx'])

# Make column type categorical and specify factor levels and order
# For some reason this is not taken over when compiling the anndata object ...
meta["sample"] = pd.Categorical(meta["sample"],
                                categories=meta["sample"].values)
meta["group"] = pd.Categorical(meta["group"],
                               categories=['Ctrl', 'A8', 'A9', 'A8A9'])

# %%
# Old run
# Get metadata from samplesheet
meta_old = pd.read_csv('../../tables/samplesheet_old.csv')
meta_old.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# Reorder by Sample ID, drop double columns, and update index
sample_idx = []
for s in meta_old['sample'].values:
    sample_idx.append(int(s[2:]))

meta_old['sample_idx'] = sample_idx
meta_old = meta_old.sort_values('sample_idx')
meta_old.index = meta_old.sample_idx
meta_old=meta_old.drop(columns=['sample_idx'])

# %% [markdown]
# ### Save/Load meta table as/from csv

# %%
pd.concat([meta,meta_old])

# %%
# Make single table with all samples old/new
# Save to path
pd.concat([meta,meta_old]).to_csv('../../results/produced_data/scRNAseq/metadata.csv')

# Load from path
pd.read_csv('../../results/produced_data/scRNAseq/metadata.csv', index_col=0)

# %% [markdown]
# ### Load features.tsv.gz from single sample

# %%
# This should be equivalent for every sample
features_path = f"/data/projects/2022/Adolph-scRNA-organoids/01_nfcore_scrnaseq/cellranger/sample-AJ10/outs/filtered_feature_bc_matrix/features.tsv.gz"
genes = pd.read_csv(gzip.open(features_path, mode="rt"), delimiter="\t", header=None)
genes.drop(genes.columns[[2]], axis=1, inplace=True)
genes.rename(columns={0:'gene_id', 1: 'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)

# %% [markdown]
# ### Get marker genes for cell types

# %%
marker = pd.read_csv('../../tables/marker_genes.csv', sep=";")

# %%
# subset by cell type marker and store as list
dict_marker={}
for m in marker.marker.unique():
    dict_marker[m] = list(marker.loc[(marker.marker==m),'Human'].dropna().values)

# %%
marker.marker.unique()

# %%
dict_marker['M_cells']

# %%
# Not sure how to handle genes that were modified with var_names_make_unique()
# When using the genes from the features.tsv.gz I should get the actual interaction - before any filtering though
intersect_genes = list(set(marker.Human.dropna().values).intersection(genes.index))

# %% [markdown]
# ## Load cellranger output to single adata

# %% [markdown]
# ### 1. adata.obs `gene_ids as index`, symbols as column
# This way there should be no identical indices.
# <br>
# **Note:** Mitochondrial genes do not work as only symbols have `MT-` string in name!
# <br>
# ### **-> Skip this block and use symbols as indices below!!!**

# %%
# Version 1: Use ensembl gene_ids as index
# Iterate over sample h5ads and concatenate
p_dir="/data/projects/2022/Adolph-scRNA-organoids"
adatas = dict()
key_save_l = []
# Use conditional statement to discriminate between new and old samples
for sample in pd.concat([meta,meta_old]).to_dict(orient="records"):
    if 'AJ' in sample['sample']:
        p_h5ads = f"{p_dir}/01_nfcore_scrnaseq/cellranger/sample-{sample['sample']}/outs/filtered_feature_bc_matrix.h5"
    elif 'TA' in sample['sample']:
        p_h5ads = f"{p_dir}/01_cellranger/cellranger/sample-{sample['sample']}/outs/filtered_feature_bc_matrix.h5"
    else:
         print(sample['sample'])
    tmp_adata = sc.read_10x_h5(p_h5ads)
                                             
    # save gene conversion key and switch index to ensembl ids before making unique
    key_save_l.append(tmp_adata.var.copy())
    tmp_adata.var['gene_symbols'] = tmp_adata.var.index
    tmp_adata.var.index = tmp_adata.var.gene_ids
    tmp_adata.var = tmp_adata.var.drop(columns=['gene_symbols','feature_types','genome'])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas[sample['sample']] = tmp_adata # assign sample_id to barcodes
    
# when concatenating all, columns in .var are somehow dropped
# index_unique in .concat appends sample ids to barcodes
adata = ad.concat(adatas, index_unique="_")

for k in key_save_l[1:]:
    assert np.all(k==key_save_l[0])

key = key_save_l[-1]
key= key.reset_index()
key.index=key.gene_ids

# Use conversion key to re-assign symbols to ensembl ids
adata.var['gene_symbols'] = key.loc[adata.var.index]['index']

# %% [markdown]
# ### 2. `Symbols as index`, gene_ids as column
# Dup symbols are appended with -1, -2, etc. using `var_names_make_unique()`
# <br>
# ### **-> Use this for QC!!!**

# %%
# Iterate over sample h5ads and concatenate
p_dir="/data/projects/2022/Adolph-scRNA-organoids"
adatas = dict()
key_save_l = []
# Use conditional statement to discriminate between new and old samples
for sample in pd.concat([meta,meta_old]).to_dict(orient="records"):
    if 'AJ' in sample['sample']:
        p_h5ads = f"{p_dir}/01_nfcore_scrnaseq/cellranger/sample-{sample['sample']}/outs/filtered_feature_bc_matrix.h5"
    elif 'TA' in sample['sample']:
        p_h5ads = f"{p_dir}/01_cellranger/cellranger/sample-{sample['sample']}/outs/filtered_feature_bc_matrix.h5"
    else:
         print(sample['sample'])
    tmp_adata = sc.read_10x_h5(p_h5ads)
                                             
    # save gene conversion key and switch index to ensembl ids before making unique
    key_save_l.append(tmp_adata.var.copy())
    tmp_adata.var = tmp_adata.var.drop(columns=['feature_types','genome'])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas[sample['sample']] = tmp_adata # assign sample_id to barcodes
    
# when concatenating all, columns in .var are somehow dropped
# index_unique in .concat appends sample ids to barcodes
adata = ad.concat(adatas, index_unique="_")

for k in key_save_l[1:]:
    assert np.all(k==key_save_l[0])

key = key_save_l[-1]

# Use conversion key to re-assign symbols to ensembl ids
adata.var.loc[key.index,'gene_ids'] = key.gene_ids

# %% [markdown]
# ### Filter min requirements
# -> before any refined filtering only remove barcodes with less than 200 genes and genes found in less than 3 cells

# %%
# Basic filter thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# %% [markdown]
# ### Specify `adata.obs` factor levels/order

# %%
adata.obs["group"] = pd.Categorical(adata.obs["group"],
                                     categories=['Ctrl', 'A8', 'A9', 'A8A9'])
adata.obs["sample"] = pd.Categorical(adata.obs["sample"],
                                     categories=adata.obs["sample"].unique())
adata.obs["patient"] = pd.Categorical(adata.obs["patient"])
adata.obs["batch"] = pd.Categorical(adata.obs["batch"])

# %%
# Append "sample_counts" as column
# First initialize with NAN and loop to fill values
adata.obs['sample_counts'] = np.NaN
for s in adata.obs["sample"].cat.categories:
    index_w_sel_sample = adata.obs.where(adata.obs['sample']==s).dropna(how='all').index
    adata.obs.loc[index_w_sel_sample,"sample_counts"] = adata.obs['sample'].value_counts().loc[s]
adata.obs['sample_counts'] = adata.obs['sample_counts'].astype(int)

# %%
# Create new label column for plots: patient-batch + group + sample_counts
adata.obs["label"] = adata.obs['patient'].astype(str) + '-' + adata.obs['batch'].astype(str) + '\n' + adata.obs['group'].astype(str) + '\nn=' + adata.obs['sample_counts'].astype(str)

# Get label factor levels as list in right order
label_l = []
for p in adata.obs.patient.cat.categories:
    for g in adata.obs.group.cat.categories:
        _sel_adata = adata.obs[(adata.obs.patient==p)&(adata.obs.group==g)]
        c = _sel_adata['sample_counts'].unique()[0]
        b = _sel_adata['batch'].unique()[0]
        label_l.append(f'{p}-{b}\n{g}\nn={c}')
        
# Use created list to relevel label factor
adata.obs["label"] = pd.Categorical(adata.obs["label"],
                                        categories=label_l)

# %%
# Create new label: patient + group
adata.obs['label1'] = adata.obs['patient'].astype(str) + '\n' + adata.obs['group'].astype(str)

# Get label factor levels as list in right order
label_l = []
for p in adata.obs.patient.cat.categories:
    for g in adata.obs.group.cat.categories:
        label_l.append(f'{p}\n{g}')
        
# Use created list to relevel label factor
adata.obs["label1"] = pd.Categorical(adata.obs["label1"],
                                     categories=label_l)

# %%
# Create new label: patient-batch + group
adata.obs["label2"] = adata.obs['patient'].astype(str) + '-' + adata.obs['batch'].astype(str) + '\n' + adata.obs['group'].astype(str)

# Get label factor levels as list in right order
label_l = []
for p in adata.obs.patient.cat.categories:
    for g in adata.obs.group.cat.categories:
        for b in adata.obs.batch.cat.categories:
            label_l.append(f'{p}-{b}\n{g}')
        
# Use created list to relevel label factor
adata.obs["label2"] = pd.Categorical(adata.obs["label2"],
                                     categories=label_l)

# %% [markdown] tags=[]
# ### Relevel factor of existing adata

# %%
# Add .copy() to mute warning
adata.obs["group"] = adata.obs.group.cat.reorder_categories(['Ctrl', 'A8', 'A9', 'A8A9'])#.copy()

# %% [markdown]
# ### Look at dup indices

# %% [markdown]
# #### Check unique gene ids/symbols(indices)

# %%
# Number of unique gene ids
len(adata.var['gene_ids'].unique())

# %%
# Number of unique gene symbols/indices
len(adata.var.index.unique())

# %%
# Apparently some gene symbols are not unique
adata.shape[1] - len(adata.var.index.unique())

# %%
# check if any indices contain a certain string (look for "-" as .var_names_make_unique appends -1, -2, etc. to duplicated gene symbols)
for g in adata.var.index:
    if '-' in g:
        print(g)

# %% [markdown]
# ### Summary stats raw adata

# %%
# Dimensions of adata - barcodes X Genes
adata.shape

# %%
# Orders highest to lowest
print(adata.obs['sample'].value_counts())
print('')
print(adata.obs['patient'].value_counts())
print('')
print(adata.obs['group'].value_counts())
print('')
print(adata.obs['sex'].value_counts())
print('')
print(adata.obs['batch'].value_counts())

# %%
# Look at .obs subset 
adata.obs[['label', 'sample_counts', 'patient']]

# %% [markdown]
# ### Plots: All samples - Basic summary stats
# - highest_expr_genes
# - violin plots:
#     - genes by counts
#     - total counts
#     - mitochondrial counts [per]

# %%
sc.pl.highest_expr_genes(adata, n_top=20)

# %%
# Calculate QC metrics for all samples
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %% [markdown]
# ## Subset adata by multiple samples

# %%
# Subset adata for samples 1 to 16

# Make copy, save index (barcodes) to column, set samples as index
_adata_subset = adata.copy()
_adata_subset.obs['barcodes'] = _adata_subset.obs.index
_adata_subset.obs.index = _adata_subset.obs['sample']

start=0
end=16
barcodes_sel = _adata_subset.obs.loc[_adata_subset.obs['sample'].cat.categories[start:end]].barcodes

adata.obs['col_sel'] = False
adata.obs.loc[barcodes_sel,'col_sel'] = True

adata_subset = adata[adata.obs.col_sel]

# %%
# Subset adata for only Ctrls - much easier as not selcting multiple
adata_subset_ctrl = adata[adata.obs.group == 'Ctrl']

# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# %%
fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata_subset, 'total_counts', groupby='label', size=2, log=True, cut=0, ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_total_counts_raw_sub.pdf", bbox_inches="tight")

# %%
fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata_subset, 'pct_counts_mt', groupby='label', ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_pct_counts_mt_raw_sub.pdf", bbox_inches="tight")

# %% [markdown]
# ### Plot weird samples

# %%
# Subset adata for samples 1 to 16

# Make copy, save index (barcodes) to column, set samples as index
_adata_subset = adata.copy()
_adata_subset.obs['barcodes'] = _adata_subset.obs.index
_adata_subset.obs.index = _adata_subset.obs['sample']

start=16
end=37
barcodes_sel = _adata_subset.obs.loc[_adata_subset.obs['sample'].cat.categories[start:end]].barcodes

adata.obs['col_sel'] = False
adata.obs.loc[barcodes_sel,'col_sel'] = True

adata_subset2 = adata[adata.obs.col_sel]

# %%
fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata_subset2, 'total_counts', groupby='label', size=2, log=True, cut=0, ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_total_counts_raw_sub2.pdf", bbox_inches="tight")

# %%
fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata_subset2, 'pct_counts_mt', groupby='label', ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_pct_counts_mt_raw_sub2.pdf", bbox_inches="tight")

# %% [markdown]
# ### Write AnnData object to disk

# %%
# Save h5ad
adata.write('../../results/produced_data/scRNAseq/raw_calprotectin.h5ad', compression='gzip')
# !h5ls '../../results/produced_data/scRNAseq/raw_calprotectin.h5ad'

# %%
# Save h5ad
adata_subset.write('../../results/produced_data/scRNAseq/raw_sub_calprotectin.h5ad', compression='gzip')
# !h5ls '../../results/produced_data/scRNAseq/raw_sub_calprotectin.h5ad'

# %% [markdown]
# ## Subset adata for ith sample and return adata

# %%
# Subset for any factor level in obs
sample_AJ1 = adata[adata.obs["sample"] == "AJ1", :]
sample_AJ1.obs["value"] = 0  # This makes AJ1 a “real” AnnData object

# %%
sample_AJ1

# %%
# sample_Aj1 adata can now be used in the same way as the main adata
sc.pl.violin(sample_AJ1, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# %% [markdown]
# ## Automate for every sample
# -> plot summary figures for every iundividual sample in adata seperately.

# %%
sample_d = dict()
for s in adata.obs['sample'].values.unique():
    _sampli = adata[adata.obs["sample"] == s, :]
    _sampli.obs["value"] = 0
    sample_d[s] = _sampli

# %%
sample_d[s]

# %%
for s in adata.obs['sample'].values.unique():
    print(s)
    sc.pl.violin(sample_d[s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# %%
# Filter thresholds for individual samples
sc.pp.filter_cells(sample_d['AJ1'], min_genes=200)
sc.pp.filter_cells(sample_d['AJ2'], min_genes=1000)

# %% [markdown]
# ### Concat all samples back to one adata

# %%
# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata.obs['sample'].values.unique():
    adatas_new_l.append(sample_d[s])
adata_filtered = anndata.concat(adatas_new_l, index_unique="_")

# %%
adata_filtered

# %% [markdown]
# #### Look at table for single obs

# %%
adata[adata.obs.sample == "AJ1"].var
