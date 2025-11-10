#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/1_run_strainid.log.out
#SBATCH -e logs/1_run_strainid.log.err

# Execute GenoPipe's StrainID module on all samples in the INPUTBAM directory (2024-Krebs_Science/data/BAM)

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/01_Run_GenoPipe
GENOPIPE=/path/to/GenoPipe
###

set -exo
module load samtools
module load anaconda3
source activate genopipe # load conda environment with GenoPipe dependencies

# Inputs and outputs
STRAINID=$GENOPIPE/StrainID
GENOME=$WRK/../data/hg38_files/hg38.fa
INPUTBAM=$WRK/../data/BAM
OUTPUT=$WRK/strainid

cd $WRK
[ -d logs ] || mkdir logs
[ -d $OUTPUT ] || mkdir $OUTPUT

cd $STRAINID
bash identify-Strain.sh -i $INPUTBAM -v $STRAINID/hg38_VCF/ -g $GENOME -o $OUTPUT
