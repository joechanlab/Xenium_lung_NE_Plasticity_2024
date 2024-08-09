#!/usr/bin/env python3
import argparse
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier

parser = argparse.ArgumentParser(description="ENVI to integrate paired scRNAseq and spatial data.")
parser.add_argument("st_path", help="Paths to single-cell RNAseq data h5ad.")
parser.add_argument("sc_path", help="Path to spatial data h5ad.")
parser.add_argument("celltype_path", help="Cell type column name in scRNAseq.", default="None")
parser.add_argument("st_outpath", help="Output path to single-cell RNAseq data h5ad.")
parser.add_argument("--k", help="Number of neighbors for KNN.", default=30)
parser.add_argument("--celltype", help="Cell type column name in scRNAseq.", default="None")
args = parser.parse_args()

# Load data
print('Loading data...')
st_data = sc.read_h5ad(args.st_path)
sc_data = sc.read_h5ad(args.sc_path)
cell_type_data = pd.read_table(args.celltype_path, sep='\t', header=0, index_col=0)
celltype = args.celltype.split(',')
celltype_match = None
for ct in celltype:
    if ct in cell_type_data.columns:
        celltype_match = ct
        break
sc_data.obs[celltype_match] = cell_type_data.loc[sc_data.obs_names, celltype_match]

# Preprocess data
print(f'Predicting cell types from the `{celltype_match}` column...')
knn = KNeighborsClassifier(n_neighbors=args.k)
knn.fit(sc_data.obsm['envi_latent'], sc_data.obs[celltype_match])
st_data_cell_type = knn.predict(st_data.obsm['envi_latent'])
st_data.obs['cell_type_predicted'] = st_data_cell_type
st_data.write(args.st_outpath)

print(f'Results saved at {args.st_outpath}')
