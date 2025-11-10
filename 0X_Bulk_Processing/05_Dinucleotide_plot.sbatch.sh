#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/Dinucleotide.log.out-%a
#SBATCH -e logs/Dinucleotide.log.err-%a
#SBATCH --array 1-94

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/X_Bulk_Processing/Dinucleotide_plot.txt
THREADS=4

# Dependencies
# - java
# - perl
# - python

###
set -exo
module load anaconda3
source activate bioinfo

# Fill in placeholder constants with your directories
GENOME="$WRK/data/hg38_files/hg38.fa"
# Script shortcuts
SCRIPTMANAGER="$WRK/bin/ScriptManager-v0.15.jar"
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
MOTIFSCAN=$WRK/bin/scan_FASTA_for_motif_as_binary_string.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library/
cd $WRK/X_Bulk_Processing

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }
# Determine the dinucleotide
DI=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`
IFS=';' read -ra PANELS <<< "$PANEL"

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} x ${DI}"
BASE=${DI}_${BED}

# ===============================================================================================================================

echo "Run Dinucleotode scan, make composite plot of dinucleotide"

java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $BEDFILE -o $WRK/Library/${BED}.fa
conda deactivate
python $MOTIFSCAN -i  $WRK/Library/${BED}.fa -m ${DI} -o $WRK/Library/${DI}_${BED}
conda activate bioinfo
rm $WRK/Library/${DI}_${BED}_anti.cdt
perl $COMPOSITE $WRK/Library/${DI}_${BED}_sense.cdt $WRK/Library/${DI}_${BED}.out
rm $WRK/Library/${BED}.fa

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/CDT"
    cp $WRK/Library/${DI}_${BED}_sense.cdt "$DIR/CDT/"
    mkdir -p "$DIR/Composites"
    cp $WRK/Library/${DI}_${BED}.out "$DIR/Composites"
done

rm $WRK/Library/${DI}_${BED}_sense.cdt
rm $WRK/Library/${DI}_${BED}.out





