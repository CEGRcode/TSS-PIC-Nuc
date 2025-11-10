#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/04_ProT.log.out-%a
#SBATCH -e logs/04_ProT.log.err-%a
#SBATCH --array 1-7

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/X_Bulk_Processing/ProT-plot-heatmap.txt
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
COMPOSITE=$WRK/bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
cd $WRK/X_Bulk_Processing/

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
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
[ $WINDOW == "200bp" ] && LLABEL="-100" && RLABEL="+100" && MLABEL="0"

# ===============================================================================================================================

echo "Run ProT "
java -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed -p -o $WRK/Library/${BED} $GENOME  $BEDFILE
java -jar "$SCRIPTMANAGER" figure-generation three-color -an -9 -ad -6 -ax -5 $WRK/Library/${BED}_PropT.cdt -o $WRK/Library/${BED}_PropT.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $WRK/Library/${BED}_PropT.png -f 20 -l $LLABEL -m $MLABEL -r $RLABEL -o $WRK/Library/${BED}_PropT.svg
# Make composite
perl $COMPOSITE $WRK/Library/${BED}_PropT.cdt $WRK/Library/${BED}_PropT.out

rm $WRK/Library/${BED}_PropT.png
# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Heatmap"
    mkdir -p "$DIR/CDT"
    mkdir -p "$DIR/Composites"
    cp $WRK/Library/${BED}_PropT.svg "$DIR/Heatmap/"
    cp $WRK/Library/${BED}_PropT.cdt "$DIR/CDT"
    cp $WRK/Library/${BED}_PropT.out "$DIR/Composites"
done

# Remove scaled matrix from top-level Library folder
rm $WRK/Library/${BED}_PropT.svg
rm $WRK/Library/${BED}_PropT.cdt
rm $WRK/Library/${BED}_PropT.out
rm $SCRIPTMANAGER