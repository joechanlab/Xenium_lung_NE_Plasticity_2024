#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scenvi
import pickle
import jax

parser = argparse.ArgumentParser(description="ENVI to integrate paired scRNAseq and spatial data.")
parser.add_argument("st_path", help="Paths to single-cell RNAseq data h5ad.")
parser.add_argument("sc_path", help="Path to spatial data h5ad.")
parser.add_argument("st_outpath", help="Output path to single-cell RNAseq data h5ad.")
parser.add_argument("sc_outpath", help="Output path to spatial data h5ad.")
parser.add_argument("imputation_outpath", help="Output path to imputed spatial data.")
parser.add_argument("model_outpath", help="Output path to model data checkpoint.")
parser.add_argument("--patient", help="Patient to subset.", default="None")
parser.add_argument("--downsample", help="Downsample data to this number of cells.", type = int, default=50000)
parser.add_argument("--HVG", help="Number of highly variable genes to use.", type = int, default=5000)
args = parser.parse_args()

# Downsample data
def downsample(adata, num_cells = 10000):
    if num_cells > adata.n_obs:
        return(adata)
    np.random.seed(42)
    random_indices = np.random.choice(adata.n_obs, num_cells, replace=False)
    adata_downsampled = adata[random_indices,:].copy()
    return(adata_downsampled)

# Load data
print('Loading data...')
st_data = sc.read_h5ad(args.st_path)
sc_data = sc.read_h5ad(args.sc_path)

if args.patient != "None":
    print(f'Filtering for the patient {args.patient}')
    sc_data = sc_data[sc_data.obs['patient']== args.patient]

# Remove outliers
st_data = st_data[~st_data.obs["outlier"],:]

# Preprocess data
print('Preprocessing data...')
sc_data.layers['norm'] = sc_data.X.copy()
st_data.layers['norm'] = st_data.X.copy()
# Input matrices are normalized but not log-transformed
sc_data.X = sc_data.raw.X.copy() 
sc.pp.normalize_total(sc_data, inplace=True, target_sum=10**4)
st_data.X = st_data.layers['counts']
sc.pp.normalize_total(st_data, inplace=True, target_sum=10**4)

if isinstance(sc_data.X, np.ndarray):
    sc_data.layers['log'] = np.log(sc_data.X + 1)
else:
    sc_data.layers['log'] = np.log(sc_data.X.toarray() + 1)

if isinstance(st_data.X, np.ndarray):
    st_data.layers['log'] = np.log(st_data.X + 1)
else:
    st_data.layers['log'] = np.log(st_data.X.toarray() + 1)

sc.pp.highly_variable_genes(sc_data, n_top_genes=args.HVG, layer="log")
marker_genes_df = sc.get.rank_genes_groups_df(sc_data, group = None, log2fc_min = 1, pval_cutoff = 0.05)
marker_genes_df = marker_genes_df.sort_values(['group', 'scores'], ascending=False).groupby('group').head(100)
marker_genes = marker_genes_df.names.unique()
sc_data.var['highly_variable'][marker_genes] = True

# Downsample data
downsampled = False
if st_data.shape[0] > args.downsample:
    print(f'Downsampled data of {args.downsample} cells...')
    downsampled = True
    st_data_raw = st_data.copy()
    st_data = downsample(st_data, int(args.downsample))

# ENVI Integration
print('Training ENVI...')
envi_model = scenvi.ENVI(st_data, sc_data)
envi_model.train()
with open(args.model_outpath, 'wb') as pickle_file:
    print("Saved model to ", args.model_outpath)
    pickle.dump(envi_model.params, pickle_file)

print("Saving single-cell RNAseq data to ", args.sc_outpath)
sc_data.obsm['envi_latent'] = envi_model.sc_data.obsm['envi_latent']
sc_data.X = sc_data.layers['norm']
del sc_data.layers['norm']
sc_data.write(args.sc_outpath)
del sc_data

print("Saving spatial data to ", args.sc_outpath)
if downsampled:
    print(f'Encoding on all {st_data.shape[0]} cells...')
    st_data = st_data_raw
    st_X = st_data[:,list(envi_model.overlap_genes)].X
    if not isinstance(st_X, np.ndarray):
        st_X = st_X.toarray()
    st_data.obsm['envi_latent'] = envi_model.encode(st_X, mode = 'spatial')
    envi_model.spatial_data = st_data.copy()
else:
    st_data.obsm['envi_latent'] = envi_model.spatial_data.obsm['envi_latent']

st_data.X = st_data.layers['norm']
del st_data.layers['norm']
st_data.write(args.st_outpath)
del st_data

print('Imputing genes...')
st_envi_latent = envi_model.spatial_data.obsm["envi_latent"].copy()
var_names = envi_model.sc_data.var_names
obs_names = envi_model.spatial_data.obs_names
del envi_model.spatial_data
imputation = pd.DataFrame(envi_model.decode_exp(st_envi_latent, 
                                                mode="sc", 
                                                max_batch = 128),
                                                columns=var_names, 
                                                index=obs_names)
with open(args.imputation_outpath, 'wb') as pickle_file:
    pickle.dump(imputation, pickle_file)
    print("Saved imputed spatial data to ", args.imputation_outpath)
