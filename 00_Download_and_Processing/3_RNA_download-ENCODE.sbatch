#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/GROcap_download-ENCODE.log.out
#SBATCH -e logs/GROcap_download-ENCODE.log.err

# Download capped GRO-cap data from ENCODE

# data/BAM
#   |--ENCFF112MXH_Grocap_hg38.bam
#   |--ENCFF663UAN_CoPRO_hg38.bam

### CHANGE ME
WRK=/Path/to/Title/00_Download_and_Preprocessing/
THREADS=8
###

# Dependencies
# - java
# - samtools
# - wget

set -exo
module load samtools
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Download BAM files
wget -c -O $BAMDIR/ENCFF112MXH_Grocap_hg38.bam https://www.encodeproject.org/files/ENCFF112MXH/@@download/ENCFF112MXH.bam
wget -c -O $BAMDIR/ENCFF663UAN_CoPRO_hg38.bam https://www.encodeproject.org/files/ENCFF663UAN/@@download/ENCFF663UAN.bam

# Index
samtools index $BAMDIR/ENCFF112MXH_Grocap_hg38.bam
samtools index $BAMDIR/ENCFF663UAN_CoPRO_hg38.bam
