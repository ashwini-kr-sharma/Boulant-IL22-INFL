#!/bin/bash

#-------------------------------------------------------------------------------

# Scripts for generating counts from the paired-end bulk T84 cells trteated with 
# IL22, INF-L and both
# Data genereted from Boulant Lab

# Analysis done by - Ashwini Kumar Sharma, PhD
# Last update - 3rd March, 2020

# Starting folder structure requirement inside the root directory
# NOTE: the names should be exact as shown below and dont change any script names etc

# <ROOTDIR>/data/fastq/   > containing all the fastq files (can be downloaded from GEO)
# <ROOTDIR>/data/anno/ > containing the annotation information (can be downloaded from GEO/Supplementary data)
# <ROOTDIR>/src/       > containing all the analysis scripts (can be downloaded from GitHub)

# We assume that the user has conda installed and has created a conda environment -
# conda create --name hypoxia -c bioconda r fastqc multiqc bioconductor-rsubread r-base
# A detailed environment .yaml file is also provided in GitHub

# This analysis was run on the LSF - IBM HPC system, if using a different HPC system,
# like qsub, SLURM etc the user has to modify the cluster submission code chunk accordingly

#-------------------------------------------------------------------------------

ROOTDIR='/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/'

#-------------------------------------------------------------------------------
# FastQC
#-------------------------------------------------------------------------------

# SCRIPT=$ROOTDIR$'src/01_fastQC.sh'
# mkdir -p $ROOTDIR$'analysis/FastQC'
# mkdir -p $ROOTDIR'logs/FastQC'
# 
# for i in $ROOTDIR$'data/fastq/'CG01*/*.fastq.gz ; do
# k=$(basename "$i" .fastq.gz)
#  bsub -n 1  -R "rusage[mem=2G]" -W 2:00 -J T84_fastqc \
#  -o $ROOTDIR'logs/FastQC/'$k'.out' -e $ROOTDIR'logs/FastQC/'$k'.err' \
#  "$SCRIPT $ROOTDIR$'analysis/FastQC' $i"
# done

#-------------------------------------------------------------------------------
# Download genomes
#-------------------------------------------------------------------------------

# Using the genomes downloaded for Carmon's Hypoxia project

#-------------------------------------------------------------------------------
# Rsubread indexing
#-------------------------------------------------------------------------------

# Using the Rsubread index generated for Carmon's Hypoxia project

#-------------------------------------------------------------------------------
# Rsubread alignment
#-------------------------------------------------------------------------------

# mkdir -p $ROOTDIR'data/bam/alignStats/RDS'
# indexdir='/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/data/hg38/index/'
# 
# for i in $ROOTDIR$'data/fastq/'* ; do
#   
#   for j in $i; do
#     id=$(basename $i)
#     mkdir -p $ROOTDIR'logs/Align/'$id
#     
#     f=($j/*)
# 
#     f1=${f[0]}
#     f2=${f[1]}
# 
#     echo "source activate hypoxia; Rscript $ROOTDIR$'src/02_alignmentRsubread.R' $f1 $f2 $indexdir $ROOTDIR$'data/bam/'" | \
#     bsub -n 15 -R "rusage[mem=20G]" -W 15:00 -J T84_align \
#     -o $ROOTDIR'logs/Align/'$id'/stdout.out' -e $ROOTDIR'logs/Align/'$id$'/stderr.err'
# 
#   done
# done

#-------------------------------------------------------------------------------
# Rsubread featurecounts
#-------------------------------------------------------------------------------

# mkdir -p $ROOTDIR'data/counts/'
# mkdir -p $ROOTDIR'logs/Count'
# hg38gtf='/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/data/hg38/Homo_sapiens.GRCh38.102.gtf.gz'
# 
# echo "source activate hypoxia; Rscript $ROOTDIR$'src/03_countsRsubread.R' $ROOTDIR$'data/bam' $hg38gtf $ROOTDIR'data/counts/'" | \
# bsub -n 15 -R "rusage[mem=5G]" -W 5:00 -J T84_count \
# -o $ROOTDIR'logs/Count/featureCount.out' -e $ROOTDIR'logs/Count/featureCount.err'

#-------------------------------------------------------------------------------
# MultiQC
#-------------------------------------------------------------------------------

SCRIPT=$ROOTDIR$'src/04_multiQC.sh'
mkdir -p $ROOTDIR$'analysis/MultiQC'
mkdir -p $ROOTDIR'logs/MultiQC'

bsub -n 1  -R "rusage[mem=5G]" -W 2:00 -J T84_multiqc \
-o $ROOTDIR'logs/MultiQC/MultiQC.out' -e $ROOTDIR'logs/MultiQC/MultiQC.err' \
"$SCRIPT $ROOTDIR$'analysis/MultiQC' $ROOTDIR"
