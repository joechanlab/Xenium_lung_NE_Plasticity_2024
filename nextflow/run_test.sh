module load singularity/3.7.1
module load gcc/10.2.0
module load cuda/11.7
export NXF_SINGULARITY_CACHEDIR="/lila/data/chanjlab/wangm10/work-nf-scrnaseq/singularity/"

nextflow run ./main.nf -resume -profile lilac -params-file ./params.yml -w "../work"
