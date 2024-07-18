#!/usr/bin/env python3
import argparse
import squidpy as sq
import scanpy as sc

parser = argparse.ArgumentParser(description="Preprocessing Xenium to h5ad.")
parser.add_argument("input", help="Paths to directory containing Xenium files.")
parser.add_argument("--output", required=True, help="The output h5ad file path.")

args = parser.parse_args()

adata = sc.read_h5ad(args.input)
sq.gr.spatial_neighbors(adata)
genes = adata.var_names.values
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
)
print("Computed the SVG.")

adata.write(args.output)
