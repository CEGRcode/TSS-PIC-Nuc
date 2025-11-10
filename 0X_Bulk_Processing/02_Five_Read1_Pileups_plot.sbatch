#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/02_Five_Read1_Pileups_plot.log.out-%a
#SBATCH -e logs/02_Five_Read1_Pileups_plot.log.err-%a
#SBATCH --array 1-334

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/X_Bulk_Processing/Five_Read1_Pileups_plot.txt
THREADS=4
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda3
source activate bioinfo

# Fill in placeholder constants with your directories
BAMDIR=$WRK/data/BAM

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
cd $WRK/X_Bulk_Processing
# Determine insert filters for the current job array index
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`

# Determine heatmap contrast threshold for the current job array index
CONTRAST_THRESH=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Determine scaling method
NORMALIZATION=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
BAM=`basename $BAMFILE ".bam"`
[[ -f $BAMFILE ]] || { echo "ERROR: BAM file not found: $BAMFILE"; exit 1; }
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $7}'`
IFS=';' read -ra PANELS <<< "$PANEL"

# ===============================================================================================================================

# Parse min/max params
TP_PARAMS="--cpu $THREADS -5 -1"
FILTER=""
[[ $MIN =~ "NaN" ]] || FILTER=${FILTER}-MIN${MIN}
[[ $MIN =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -n ${MIN}"
[[ $MAX =~ "NaN" ]] || FILTER=${FILTER}-MAX${MAX}
[[ $MAX =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -x ${MAX}"

# Parse normalization params and get scaling factor
FACTOR=1
if [ "$NORMALIZATION" == "NCIS" ]; then
    NFFILE=$BAMDIR/NormalizationFactors/${BAM}_NCISb_ScalingFactors.out
    if [[ ! -f $NFFILE ]]; then
        echo "ERROR: Missing Normalization factor file (NCIS): $NFFILE"
        exit 1
    fi
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
elif [ "$NORMALIZATION" == "TotalTag" ]; then
    NFFILE=$BAMDIR/NormalizationFactors/${BAM}_TotalTag_ScalingFactors.out
    if [[ ! -f $NFFILE ]]; then
        echo "ERROR: Missing Normalization factor file (TotalTag): $NFFILE"
        exit 1
    fi
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
else
    NORMALIZATION="Raw"
fi

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} x ${BAMFILE} x ${MIN} x ${MAX} x ${NORMALIZATION}"
BASE=${BAM}_${BED}_5read1${FILTER}
OUTFILE=$WRK/Library/${BASE}_Normalized.out

# ===============================================================================================================================

echo "Run Custom read1 pileup"

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE $TP_PARAMS -o $WRK/Library/${BASE}.out

# Scale matrix
java -jar $SCRIPTMANAGER read-analysis scale-matrix $WRK/Library/${BASE}.out -s $FACTOR -o "$OUTFILE" -l 1

# Remove unscaled matrix
rm $WRK/Library/${BASE}.out

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Composites"
    cp "$OUTFILE" "$DIR/Composites/"
    echo "Copied $OUTFILE to $DIR/Composites/"
done

# Remove scaled matrix from top-level Library folder
[ -f "$OUTFILE" ] && rm "$OUTFILE"
rm  $SCRIPTMANAGER
