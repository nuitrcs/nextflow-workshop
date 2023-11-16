#!/bin/bash
#SBATCH -A a9009
#SBATCH -p a9009
#SBATCH -t 01:00:00
#SBATCH -n 1 
#SBATCH -N 1
#SBATCH --mem=4GB
#SBATCH --job-name="nf-parent-RNAseq"
#SBATCH -o %j-%x.out

module purge all
module load nextflow/23.04.3
nextflow run \
  -c nextflow-with-profiles.config \
  -profile quest_slurm \
  example-rnaseq-mixed-executors.nf

