# Activate environment
source activate hypoxia

# Get arguments from main script
export OUTDIR=${1}
export FQDIR=${2}

# Set locale values
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# MultiQC
multiqc --title 'T84_IL22_IFNL_Data_BoulantLab' --outdir $OUTDIR $FQDIR
