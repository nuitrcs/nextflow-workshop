# nextflow-workshop
Reference materials for RCDS' "Introduction to Nextflow" workshop.

### Setup
Before doing anything, be sure to run the following, which will create the Singularity container used in this directory's Nextflow processes.

`bash singularity-container-setup.sh`

### Simple example

This submits four jobs to SLURM:

```
module load nextflow/23.04.3
nextflow run -c nextflow.config example-rnaseq.nf
```

### `resume` example

This will use the cached data from the previous run:

```
module load nextflow/23.04.3
nextflow run -c nextflow.config example-rnaseq.nf -resume
```

### Mixed executors example

This parent job will run two processes as Slurm jobs and two processes 'locally' (i.e., on the compute node)

```
sbatch rnaseq-parent.sh
```
