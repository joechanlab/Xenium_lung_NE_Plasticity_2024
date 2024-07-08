#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc
from anndata import AnnData


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
            result.append({"group": group, "gene": gene, "annotation": annotation})
    return result


parser = argparse.ArgumentParser(description="Preprocessing Xenium to h5ad.")
parser.add_argument("input", help="Paths to directory containing Xenium files.")
parser.add_argument("--output", required=True, help="The output h5ad file path.")

args = parser.parse_args()

transcripts = pd.read_parquet(os.path.join(args.input, "transcripts.parquet"))
cells = pd.read_parquet(os.path.join(args.input, "cells.parquet"))

# Filter for transcripts with QV score >= 20
# https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/xoa-algorithms#qvs
QV_CUTOFF = 20
pass_qv = transcripts.qv >= QV_CUTOFF
transcripts = transcripts.loc[pass_qv, :]

# filter the transcripts within nucleus
# overlap_nucleus = transcripts.overlaps_nucleus == 1
zero_nucleus_distance = transcripts.nucleus_distance == 0
transcripts = transcripts.loc[zero_nucleus_distance, :]

# Create AnnData object
cell_gene_matrix = pd.crosstab(transcripts["cell_id"], transcripts["feature_name"])
obs = pd.DataFrame({"cell_id": cell_gene_matrix.index})
obs = obs.merge(cells.loc[:, ["cell_id", "x_centroid", "y_centroid", "nucleus_area"]])
obs.index = cell_gene_matrix.index
sparse_matrix = csr_matrix(cell_gene_matrix.values)
adata = AnnData(X=sparse_matrix, obs=obs, var=cell_gene_matrix.columns.to_frame())
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
print("Created the AnnData object.")

# Preprocessing
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True, target_sum=10**4)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
print("Preprocessed the data.")

adata.write(args.output)
