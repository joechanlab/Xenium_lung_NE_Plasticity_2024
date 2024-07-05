#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
from anndata import AnnData
import anndata

def find_key_by_value(dictionary, target_value):
    for key, value_list in dictionary.items():
        if target_value in value_list:
            return key
    return None

def find_keys_for_all_values(annotation_dict, search_dict):
    result = []
    for group, genes in search_dict.items():
        annotations = [find_key_by_value(annotation_dict, gene) for gene in genes]
        for gene, annotation in zip(genes, annotations):
            result.append({'group': group, 'gene': gene, 'annotation': annotation})
    return result

basedir = "/data/chanjlab/Xenium.lung.NE_plasticity.010124/preprocessing/20231013__183137__JC_PILOT_RU151_RU263_RESECTION/"
sample_name = "output-XETG00065__0005465__Region_1__20231013__183151" # output-XETG00065__0005504__RU151__20231013__183151

transcripts = pd.read_parquet(os.path.join(basedir, sample_name, "transcripts.parquet"))
cells = pd.read_parquet(os.path.join(basedir, sample_name, "cells.parquet"))

# Filter for transcripts with QV score >= 20
# https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/xoa-algorithms#qvs
QV_CUTOFF = 20
pass_qv = transcripts.qv >= QV_CUTOFF
transcripts = transcripts.loc[pass_qv,:]

overlap_nucleus = transcripts.overlaps_nucleus == 1
zero_nucleus_distance = transcripts.nucleus_distance == 0
transcripts = transcripts.loc[zero_nucleus_distance,:]

cell_gene_matrix = pd.crosstab(transcripts['cell_id'], transcripts['feature_name'])

obs = pd.DataFrame({'cell_id': cell_gene_matrix.index})
obs = obs.merge(cells.loc[:,['cell_id','x_centroid','y_centroid','nucleus_area']])
obs.index = cell_gene_matrix.index

# Convert to sparse matrix format
sparse_matrix = csr_matrix(cell_gene_matrix.values)

# Create AnnData object
adata = AnnData(X=sparse_matrix, obs=obs, var=cell_gene_matrix.columns.to_frame())

adata.obsm['spatial'] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)

sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],
    wspace=0.4,
)

adata.write('/scratch/wangm10/RU1124.h5ad')

# Get markers
gene_info = pd.read_csv("/lila/data/chanjlab/Xenium.lung.NE_plasticity.010124/ref/Final_Gene_Panel_Design.081523/Xenium_hLung_v1_metadata.csv")
gene_info_2 = pd.read_table("/lila/data/chanjlab/Xenium.lung.NE_plasticity.010124/ref/Final_Gene_Panel_Design.081523/Xenium.custom_gene_panel.v5.txt")
annotation_dict = gene_info.groupby('Annotation')['Gene'].apply(list).to_dict()
annotation_dict_2 = gene_info_2.groupby('Info')['Gene'].apply(list).to_dict()
annotation_dict.update(annotation_dict_2)

basedir = "/scratch/wangm10/RU1124/"
adata = anndata.read_h5ad(os.path.join(basedir, 'RU1124_celltypist.h5ad'))

ranked_genes = adata.uns['rank_genes_groups']
gene_names = ranked_genes['names']
gene_list = gene_names.tolist()
genes_dict = {}
for group, genes in zip(ranked_genes['names'].dtype.names, gene_list):
    genes_dict[group] = genes

# Get the result as a list of dictionaries
table_data = find_keys_for_all_values(annotation_dict, genes_dict)

# Convert to pandas DataFrame
df = pd.DataFrame(table_data)

df.to_csv('/scratch/wangm10/RU1124.h5ad')
