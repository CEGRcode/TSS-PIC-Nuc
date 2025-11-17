#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/03_Five_Read1_Pileups_heatmap.log.out-%a
#SBATCH -e logs/03_Five_Read1_Pileups_heatmap.log.err-%a
#SBATCH --array 6-29

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/0X_Bulk_Processing/Five_Read1_Pileups_heatmap.txt
THREADS=4
###

# Dependencies
# - java
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
cd $WRK/0X_Bulk_Processing

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
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $7}'`
IFS=';' read -ra PANELS <<< "$PANEL"

# ===============================================================================================================================

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

# Set window labels
RLABEL="R"
LLABEL="L"
MLABLE="M"
WINDOW=`echo $BED | rev | cut -d"_" -f1 | rev`
[ $WINDOW == "1000bp" ] && LLABEL="-500" && RLABEL="+500" && MLABEL="0"
[ $WINDOW == "500bp" ] && LLABEL="-250" && RLABEL="+250" && MLABEL="0"
[ $WINDOW == "250bp" ] && LLABEL="-125" && RLABEL="+125" && MLABEL="0"
[ $WINDOW == "200bp" ] && LLABEL="-100" && RLABEL="+100" && MLABEL="0"
[ $WINDOW == "100bp" ] && LLABEL="-30" && RLABEL="+70" && MLABEL="+20"
# Parse min/max params
TP_PARAMS="--cpu $THREADS -5 -1"
FILTER=""
[[ $MIN =~ "NaN" ]] || FILTER=${FILTER}-MIN${MIN}
[[ $MIN =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -n ${MIN}"
[[ $MAX =~ "NaN" ]] || FILTER=${FILTER}-MAX${MAX}
[[ $MAX =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -x ${MAX}"

# Parse heatmap params
HM_PARAMS="-p 0.95"
[[ $CONTRAST_THRESH =~ "NaN" ]] || HM_PARAMS="-a $CONTRAST_THRESH"

# Parse normalization params and get scaling factor
FACTOR=1
if [ "$NORMALIZATION" == "NCIS" ]; then
    NFFILE=$WRK/data/NormalizationFactors/$BAM\_NCISb_ScalingFactors.out
    [[ -f $NFFILE ]] || exit "Missing Normalization factor file (NCIS): $NFFILE"
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
elif [ "$NORMALIZATION" == "TotalTag" ]; then
    NFFILE=$WRK/data/NormalizationFactors/$BAM\_TotalTag_ScalingFactors.out
    [[ -f $NFFILE ]] || exit "Missing Normalization factor file (TotalTag): $NFFILE"
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
else
    NORMALIZATION="Raw"
fi

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} x ${BAMFILE} x ${MIN} x ${MAX} x ${NORMALIZATION}"
BASE=${BAM}_${BED}_5read1${FILTER}
OUTFILE=$WRK/Library/${BASE}_Normalized_merge_label.svg

# ===============================================================================================================================

echo "Run Custom read1 pileup"

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE $TP_PARAMS -M $WRK/Library/$BASE

# Scale matrix
java -jar $SCRIPTMANAGER read-analysis scale-matrix $WRK/Library/${BASE}_anti.cdt  -s $FACTOR -o $WRK/Library/${BASE}_Normalized_anti.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $WRK/Library/${BASE}_sense.cdt -s $FACTOR -o $WRK/Library/${BASE}_Normalized_sense.cdt

rm $WRK/Library/${BASE}_anti.cdt 
rm $WRK/Library/${BASE}_sense.cdt
# Two-color heatmap & merge
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS --blue $WRK/Library/${BASE}_Normalized_sense.cdt -o $WRK/Library/${BASE}_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS --red  $WRK/Library/${BASE}_Normalized_anti.cdt  -o $WRK/Library/${BASE}_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $WRK/Library/${BASE}_sense_treeview.png $WRK/Library/${BASE}_anti_treeview.png -o $WRK/Library/${BASE}_merge.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $WRK/Library/${BASE}_merge.png \
	-l $LLABEL -m $MLABEL -r $RLABEL -w 1 -f 20 -o "$OUTFILE"

rm $WRK/Library/${BASE}_sense_treeview.png
rm $WRK/Library/${BASE}_anti_treeview.png
rm $WRK/Library/${BASE}_merge.png

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Heatmap"
    mkdir -p "$DIR/CDT"
    cp "$OUTFILE" "$DIR/Heatmap/"
    cp $WRK/Library/${BASE}_Normalized_anti.cdt "$DIR/CDT"
    cp $WRK/Library/${BASE}_Normalized_sense.cdt "$DIR/CDT"
    echo "Copied $OUTFILE to $DIR/Heatmap/ and $DIR/CDT" 
done

# Remove scaled matrix from top-level Library folder
[ -f "$OUTFILE" ] && rm "$OUTFILE"
rm $WRK/Library/${BASE}_Normalized_anti.cdt
rm $WRK/Library/${BASE}_Normalized_sense.cdt

rm $SCRIPTMANAGER
