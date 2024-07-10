module load singularity/3.7.1
module load gcc/10.2.0
module load cuda/11.7
export NXF_SINGULARITY_CACHEDIR="/lila/data/chanjlab/wangm10/work-nf-scrnaseq/singularity/"

nextflow run ./main.nf -profile lilac -params-file ../data/params.yml -w "../work"
