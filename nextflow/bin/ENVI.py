#!/usr/bin/env python3
import argparse
import scanpy as sc
import scenvi
import numpy as np

parser = argparse.ArgumentParser(description="ENVI to integrate paired scRNAseq and spatial data.")
parser.add_argument("st_path", help="Paths to single-cell RNAseq data h5ad.")
parser.add_argument("sc_path", help="Path to spatial data h5ad.")
parser.add_argument("st_outpath", help="Output path to single-cell RNAseq data h5ad.")
parser.add_argument("sc_outpath", help="Output path to spatial data h5ad.")
parser.add_argument("--downsample", help="Downsample data to this number of cells.", default=None, type=int)
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

# Downsample data
if args.downsample is not None:
    print('Downsampling data...')
    st_data = downsample(st_data, args.downsample)
    sc_data = downsample(sc_data, args.downsample)

# Preprocess data
print('Preprocessing data...')
if 'counts' not in sc_data.layers:
    sc_data.layers['counts'] = sc_data.raw.X.copy()
if 'counts' not in st_data.layers:
    st_data.layers['counts'] = st_data.raw.X.copy()
sc_data.layers['log'] = np.log(sc_data.layers['counts'].toarray() + 1)
st_data.layers['log'] = np.log(st_data.layers['counts'].toarray() + 1)

# ENVI Integration
print('Training ENVI...')
envi_model = scenvi.ENVI(st_data, sc_data)
envi_model.train()

print('Imputing genes...')
envi_model.impute_genes()

print('Inferring niche covariates...')
envi_model.infer_niche_covet()

# Save results
print('Saving results...')
st_data.obsm['envi_latent'] = envi_model.spatial_data.obsm['envi_latent']
st_data.obsm['COVET'] = envi_model.spatial_data.obsm['COVET']
st_data.obsm['COVET_SQRT'] = envi_model.spatial_data.obsm['COVET_SQRT']
st_data.uns['COVET_genes'] =  envi_model.CovGenes
st_data.obsm['imputation'] = envi_model.spatial_data.obsm['imputation']

sc_data.obsm['envi_latent'] = envi_model.sc_data.obsm['envi_latent']
sc_data.obsm['COVET'] = envi_model.sc_data.obsm['COVET']
sc_data.obsm['COVET_SQRT'] = envi_model.sc_data.obsm['COVET_SQRT']
sc_data.uns['COVET_genes'] =  envi_model.CovGenes

st_data.write(args.st_outpath)
sc_data.write(args.sc_outpath)
print(f'Results saved at {args.st_outpath} and {args.sc_outpath}')
