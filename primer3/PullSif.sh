export SINGULARITY_CACHEDIR=/tmp/$USER/apptainer/cache
export SINGULARITY_TMPDIR=/tmp/$USER/apptainer/cache
export APPTAINER_CACHEDIR=/tmp/$USER/apptainer/cache
export APPTAINER_TMPDIR=//tmp/$USER/apptainer/cache
mkdir -p $APPTAINER_TMPDIR
singularity pull primer3v2.5.sif docker://quay.io/biocontainers/primer3:v2.5.0--pl526he1b5a44_0
