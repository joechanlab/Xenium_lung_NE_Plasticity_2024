import numpy as np
import pandas as pd
import tifffile
import dask.array as da
from matplotlib import pyplot as plt
import sopa.io
import sopa.segmentation
from xarray import DataArray
from spatialdata.transformations import Identity, set_transformation, get_transformation
from spatialdata_io.readers.xenium import xenium
import spatialdata_plot
from pathlib import Path
import shutil

parser = argparse.ArgumentParser(
    description="Cellpose segmentation."
)
parser.add_argument("input", help="The input Xenium path")
parser.add_argument("output", help="The output spatialdata file")
parser.add_argument("workdir", help="Working directory.")

args = parser.parse_args()

sdata = xenium(args.input, 
               cells_boundaries = True,
               nucleus_boundaries = True,
               cells_as_circles = False,
               cells_labels = False,
               nucleus_labels = False,
               transcripts = True,
               morphology_mip = False,
               morphology_focus = True,
               aligned_images = True,
               cells_table = False
              )

# nucleus segmentation
patches = sopa.segmentation.Patches2D(sdata, 
                                      'if_image', 
                                      patch_width=1000,
                                      patch_overlap=100)
out = patches.write()

# cell boundaries
method_cyto = sopa.segmentation.methods.cellpose_patch(diameter=45, 
                                                  channels=[0, 2],
                                                  model_type = 'cyto3',
                                                  cellpose_model_kwargs = {'gpu': True})
segmentation_cyto = sopa.segmentation.StainingSegmentation(sdata, 
                                                           method_cyto, 
                                                           channels = [0, 2], 
                                                      image_key = 'if_image')
for patch_index in range(len(cropped_sdata['sopa_patches'])):
    segmentation_cyto.write_patch_cells(args.workdir, patch_index)
cells = sopa.segmentation.StainingSegmentation.read_patches_cells(args.workdir)
cells = sopa.segmentation.shapes.solve_conflicts(cells)
sopa.segmentation.StainingSegmentation.add_shapes(sdata, cells, 
                                                  image_key = 'if_image', 
                                                  shapes_key = 'cellpose_cell_boundaries')

# nucleus segmentation
method_nuc = sopa.segmentation.methods.cellpose_patch(diameter=35, 
                                                      channels=[2],
                                                      model_type = 'nuclei',
                                                      cellpose_model_kwargs = {'gpu': True})
segmentation_nuc = sopa.segmentation.StainingSegmentation(sdata, 
                                                          method_nuc,
                                                          channels = [2],
                                                          image_key = 'if_image')
for patch_index in range(len(cropped_sdata['sopa_patches'])):
    segmentation_nuc.write_patch_cells(args.workdir, patch_index)
nuc = sopa.segmentation.StainingSegmentation.read_patches_cells(args.workdir)
nuc = sopa.segmentation.shapes.solve_conflicts(nuc)

sopa.segmentation.StainingSegmentation.add_shapes(sdata, nuc, 
                                                  image_key = 'if_image', 
                                                  shapes_key = 'cellpose_nucleus_boundaries')

# write results
path_write = Path(args.output)
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
