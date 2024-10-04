# module load singularity
# module load cuda
export NXF_SINGULARITY_CACHEDIR="/usersoftware/chanj3/singularity"

nextflow run ./main.nf -resume -profile iris \
    -params-file /data1/chanj3/Xenium.lung.NE_plasticity.010124/nf-xenium/inputs/params_run3.yml \
    -w /data1/chanj3/Xenium.lung.NE_plasticity.010124/nf-xenium/work
