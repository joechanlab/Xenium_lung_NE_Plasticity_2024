#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix
import statsmodels.formula.api as smf
import scanpy as sc
from anndata import AnnData

from typing import Optional, Callable

parser = argparse.ArgumentParser(description="Preprocessing Xenium to h5ad.")
parser.add_argument("input", help="Paths to directory containing Xenium files.")
parser.add_argument("output", help="The output h5ad file path.")
parser.add_argument(
    "--sample_name", default="sample_name", help="The sample column name."
)
parser.add_argument("--n_pca", default=200, help="The numbers of PC.")
args = parser.parse_args()

# Definitions of helper functions
def calculate_featcount_dist(
    adata: AnnData, *, key_added: Optional[str] = "featcount_dist"
) -> None:
    """
    Calculates feature-count distances and stores them in adata.obs[key_added].
    Feature-count distance is the difference between the observed and expected
    log no. of features given log total counts.
    From Germain et al. 2020 (DOI: 10.1186/s13059-020-02136-7).
    """
    mod = smf.ols("log1p_n_genes_by_counts ~ log1p_total_counts", adata.obs).fit()
    pred = mod.predict()
    adata.obs[key_added] = adata.obs["log1p_n_genes_by_counts"] - pred


def calculate_group_featcount_dist(
    adata: AnnData, *, group_key: str, key_added: Optional[str] = None
) -> None:
    """
    Calculates group-wise feature-count distances for groups in
    adata.obs[group_key] and stores them in adata.obs[key_added]
    """
    if key_added is None:
        key_added = f"featcount_dist_{group_key}"
    groups = adata.obs[group_key].unique().tolist()
    for group in groups:
        tmp_adata = adata[adata.obs[group_key] == group].copy()
        calculate_featcount_dist(tmp_adata, key_added="tmp_fcd")
        adata.obs.loc[adata.obs[group_key] == group, key_added] = tmp_adata.obs[
            "tmp_fcd"
        ]


def mad_outlier(
    adata: AnnData, *, metric: str, nmads_upper: float, nmads_lower: float
) -> pd.Series:
    """
    Marks cell as outlier (i.e. returns True) if it is more than
    nmads_upper/nmads_lower median absolute deviations above/below
    the median of the metric.
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads_lower * sp.stats.median_abs_deviation(M)) | (
        M > np.median(M) + nmads_upper * sp.stats.median_abs_deviation(M)
    )
    return outlier


def default_filter(adata: AnnData, group_key: str) -> np.ndarray:
    """
    Default filter adapted from Germain et al. 2020 (DOI: 10.1186/s13059-020-02136-7)
    and Heumos et al. 2023 (DOI: 10.1038/s41576-023-00586-w). Outliers must meet at
    least two of the five given criteria:
    1. log1p_total_counts < median - 3 MADs OR log1p_total_counts > median + 5 MADs
    2. log1p_n_genes_by_counts < median - 3 MADs OR log1p_n_genes_by_counts > median + 5 MADs
    3. pct_counts_in_top_20_genes < median - 5 MADs OR pct_counts_in_top_20_genes > median + 5 MADs
    4. featcount_dist_<group> < median - 5 MADs OR featcount_dist_<group> > median + 5 MADs
    5. pct_counts_mito > median + 2.5 MADs AND pct_counts_mito > 8%
    """
    condition = (
        np.sum(
            [
                mad_outlier(
                    adata, metric="log1p_total_counts", nmads_lower=3, nmads_upper=5
                ),
                mad_outlier(
                    adata,
                    metric="log1p_n_genes_by_counts",
                    nmads_lower=3,
                    nmads_upper=5,
                ),
                mad_outlier(
                    adata,
                    metric="pct_counts_in_top_20_genes",
                    nmads_lower=5,
                    nmads_upper=5,
                ),
                mad_outlier(
                    adata,
                    metric=f"featcount_dist_{group_key}",
                    nmads_lower=5,
                    nmads_upper=5,
                )
            ],
            axis=0,
        )
        >= 2
    )
    return condition


def designate_outliers(
    adata: AnnData,
    *,
    condition: Callable[[AnnData, str], np.ndarray],
    group_key: str,
    key_added: Optional[str] = "outlier",
) -> None:
    """
    Marks group-wise outliers according to condition given as a function
    of adata and group_key, and saves them to adata.obs[key_added]
    """
    groups = adata.obs[group_key].unique().tolist()
    for group in groups:
        adata.obs.loc[adata.obs[group_key] == group, key_added] = condition(
            adata[adata.obs[group_key] == group], group_key
        )
    adata.obs[key_added] = adata.obs[key_added].astype(bool)


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
adata.obs['sample_name'] = args.sample_name
print("Created the AnnData object.")

# Filtering
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True, target_sum=10**4)
sc.pp.log1p(adata)

# Remove outliers
calculate_group_featcount_dist(adata, group_key='sample_name')
designate_outliers(adata, condition=default_filter, group_key='sample_name')

# Processing
sc.pp.pca(adata, n_comps=args.n_pca)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

print("Preprocessed the data.")

adata.write(args.output)
