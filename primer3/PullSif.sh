export SINGULARITY_CACHEDIR=/tmp/$USER/apptainer/cache
export SINGULARITY_TMPDIR=/tmp/$USER/apptainer/cache
export APPTAINER_CACHEDIR=/tmp/$USER/apptainer/cache
export APPTAINER_TMPDIR=//tmp/$USER/apptainer/cache
mkdir -p $APPTAINER_TMPDIR
singularity build primer3v2.6.1.sif docker://quay.io/biocontainers/primer3:2.6.1--pl5321h503566f_7
