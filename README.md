# Xenium Lung NE Plasticity 

This repository contains the nextflow pipeline and scripts to analyze Xenium data for lung NE plasticity project. 

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows, where each row contains the sample name, path to Xenium Ranger output.

`samplesheet.csv`:
```csv
sample, xenium
CONTROL_REP1, xenium_output_folder
```

Next, prepare a parameter YAML file that looks as follows:

`params.yml`:
```yaml
samplesheet: "../results/samplesheet.csv"             # path to the sample sheet
outdir: "../results/"                                 # output folder
experiment:
  name: "Xenium"
QC:
  gene_information_file: "../data/annotation.csv"     # gene annotation file
  extra_markers: "PHOX2B"                             # custom markers to plot
max_memory: "36.GB"                                   # memory
max_cpus: 6                                           # cpu
```

Now, you can run the pipeline. For local run on HPC with singularity installed, execute the following command

```bash
nextflow run ./main.nf \
   -profile singularity \
   -params-file ./params.yml \
   -w ./work/
```

If you are using MSKCC lilac, you can use the pre-defined `lilac` profile that uses the LSF executor.
```bash
nextflow run ./main.nf \
   -profile lilac \
   -params-file ./params.yml \
   -w ./work/
```
