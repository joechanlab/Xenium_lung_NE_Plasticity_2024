#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scenvi
import pickle
# from scenvi.utils import compute_covet
from sklearn.neighbors import KNeighborsClassifier

parser = argparse.ArgumentParser(description="ENVI to integrate paired scRNAseq and spatial data.")
parser.add_argument("st_path", help="Paths to single-cell RNAseq data h5ad.")
parser.add_argument("sc_path", help="Path to spatial data h5ad.")
parser.add_argument("st_outpath", help="Output path to single-cell RNAseq data h5ad.")
parser.add_argument("sc_outpath", help="Output path to spatial data h5ad.")
parser.add_argument("model_outpath", help="Output path to model data checkpoint.")
parser.add_argument("--patient", help="Patient to subset.", default="None")
parser.add_argument("--downsample", help="Downsample data to this number of cells.", type = int, default=50000)
parser.add_argument("--celltype", help="Cell type column name in scRNAseq.", default="None")
parser.add_argument("--k", help="Number of neighbors for KNN.", default=30)
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
st_data = st_data[~st_data.obs["outlier"],:].copy()

# Preprocess data
print('Preprocessing data...')
sc_data.layers['norm'] = sc_data.X.copy()
st_data.layers['norm'] = st_data.X.copy()
sc_data.X = sc_data.raw.X.copy()
st_data.X = st_data.layers['counts']

if isinstance(sc_data.X, np.ndarray):
    sc_data.layers['log'] = np.log(sc_data.X + 1)
else:
    sc_data.layers['log'] = np.log(sc_data.X.toarray() + 1)

if isinstance(st_data.X, np.ndarray):
    st_data.layers['log'] = np.log(st_data.X + 1)
else:
    st_data.layers['log'] = np.log(st_data.X.toarray() + 1)

# Downsample data
downsampled = False
if st_data.shape[0] > args.downsample:
    print('Downsampling data...')
    downsampled = True
    st_data_raw = st_data.copy()
    st_data = downsample(st_data, int(args.downsample))
    sc_data = downsample(sc_data, int(args.downsample))

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
sc_data.obsm['envi_latent'] = envi_model.sc_data.obsm['envi_latent']
sc_data.obsm['COVET'] = envi_model.sc_data.obsm['COVET']
sc_data.obsm['COVET_SQRT'] = envi_model.sc_data.obsm['COVET_SQRT']
sc_data.uns['COVET_genes'] =  envi_model.CovGenes
sc_data.X = sc_data.layers['norm']
del sc_data.layers['norm']
sc_data.write(args.sc_outpath)

if downsampled:
    st_data_raw.obsm['envi_latent'] = envi_model.encode(st_data_raw[:,list(envi_model.overlap_genes)].layers['log'], 
                                                        mode = 'spatial')
    # st_data_raw.uns['COVET_genes'] =  envi_model.CovGenes
    # st_data.obsm['imputation'] = pd.DataFrame(
    #         envi_model.decode_exp(st_data_raw.obsm["envi_latent"], mode="sc"),
    #         columns=envi_model.sc_data.var_names,
    #         index=envi_model.spatial_data.obs_names,
    #     )
    st_data = st_data_raw
else:
    st_data.obsm['envi_latent'] = envi_model.spatial_data.obsm['envi_latent']
    st_data.obsm['COVET'] = envi_model.spatial_data.obsm['COVET']
    st_data.obsm['COVET_SQRT'] = envi_model.spatial_data.obsm['COVET_SQRT']
    st_data.uns['COVET_genes'] =  envi_model.CovGenes
    st_data.obsm['imputation'] = envi_model.spatial_data.obsm['imputation']

if args.celltype != "None":
    print('Predicting cell types...')
    knn = KNeighborsClassifier(n_neighbors=args.k)
    celltype = args.celltype.split(',')
    celltype_match = None
    for ct in celltype:
        if ct in sc_data.obs.columns:
            celltype_match = ct
            break
    knn.fit(sc_data.obsm['envi_latent'], sc_data.obs[celltype_match])
    st_data_cell_type = knn.predict(st_data.obsm['envi_latent'])
    st_data.obs['cell_type_predicted'] = st_data_cell_type
st_data.X = st_data.layers['norm']
del st_data.layers['norm']
st_data.write(args.st_outpath)

# Save model parameters
with open(args.model_outpath, 'wb') as pickle_file:
    pickle.dump(envi_model.params, pickle_file)

print(f'Results saved at {args.st_outpath}, {args.sc_outpath} and {args.model_outpath}.')
