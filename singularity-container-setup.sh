#!/bin/bash

# save tutorial docker image as local singularity container
module load singularity
singularity build rnaseq-nf.simg docker://nextflow/rnaseq-nf
