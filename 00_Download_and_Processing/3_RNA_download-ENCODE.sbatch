#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/GROcap_download-ENCODE.log.out
#SBATCH -e logs/GROcap_download-ENCODE.log.err

# Download capped GRO-cap data from ENCODE

# data/BAM
#   |--GRO-cap_ENCODE_merge_hg38.bam

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
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