#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/04_4-color-heatmap.log.out-%a
#SBATCH -e logs/04_4-color-heatmap.log.err-%a
#SBATCH --array 1-8

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/X_Bulk_Processing/4-color-plot.txt
THREADS=4
###
set -exo
module load anaconda3
source activate bioinfo

# Fill in placeholder constants with your directories
GENOME="$WRK/data/hg38_files/hg38.fa"
# Script shortcuts
SCRIPTMANAGER="$WRK/bin/ScriptManager-v0.15.jar"
RESIZE=$WRK/bin/resize_png.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library/

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
[ $WINDOW == "6bp" ] && LLABEL="-3" && RLABEL="+2" && MLABEL="0"
[ $WINDOW == "60bp" ] && LLABEL="-10" && RLABEL="+50" && MLABEL="+20"
[ $WINDOW == "32bp" ] && LLABEL="-15" && RLABEL="+15" && MLABEL="0"
[ $WINDOW == "200bp" ] && LLABEL="-100" && RLABEL="+100" && MLABEL="0"

# ===============================================================================================================================

# MAIN CONDITIONAL BLOCK
if    [[ "$WINDOW" == "200bp" ]]; then
    # --- 200bp 2 color case---- #
    java -jar "$SCRIPTMANAGER" sequence-analysis fasta-extract "$GENOME" "$BEDFILE" -o "$WRK/Library/${BED}.fa"

    # Generate four-color plot from sequence
    java -jar "$SCRIPTMANAGER" figure-generation four-color -c 2EFE99 2EFE99 0000FF 0000FF "$WRK/Library/${BED}.fa" -o "$WRK/Library/${BED}.png"

    # Resize PNG to standard dimensions
    python "$RESIZE" -r 600 -c 200 -i "$WRK/Library/${BED}.png" -o "$WRK/Library/${BED}_RESIZE.png"

    # Generate labeled heatmap
    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$WRK/Library/${BED}_RESIZE.png" \
        -l "$LLABEL" -m "$MLABEL" -r "$RLABEL" -w 1 -f 20 -o "$WRK/Library/${BED}_RESIZE.svg"

    # Cleanup
    rm -f "$WRK/Library/${BED}_RESIZE.png" "$WRK/Library/${BED}.fa"

if [[ "$WINDOW" != "32bp" ]]; then
    # --- standard case ---
    java -jar "$SCRIPTMANAGER" sequence-analysis fasta-extract "$GENOME" "$BEDFILE" -o "$WRK/Library/${BED}.fa"

    # Generate four-color plot from sequence
    java -jar "$SCRIPTMANAGER" figure-generation four-color "$WRK/Library/${BED}.fa" -o "$WRK/Library/${BED}.png"

    # Resize PNG to standard dimensions
    python "$RESIZE" -r 600 -c 200 -i "$WRK/Library/${BED}.png" -o "$WRK/Library/${BED}_RESIZE.png"

    # Generate labeled heatmap
    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$WRK/Library/${BED}_RESIZE.png" \
        -l "$LLABEL" -m "$MLABEL" -r "$RLABEL" -w 1 -f 20 -o "$WRK/Library/${BED}_RESIZE.svg"

    # Cleanup
    rm -f "$WRK/Library/${BED}_RESIZE.png" "$WRK/Library/${BED}.fa"

else
    # --- 32bp special case ---
    awk 'BEGIN{OFS="\t"} {if($6=="-") print $1,$2+1,$3,$4,$5,$6; else print $1,$2,$3-1,$4,$5,$6}' \
        "$WRK/Library/${BED}.bed" > "$WRK/Library/${BED}_31bp.bed"

    java -jar "$SCRIPTMANAGER" sequence-analysis fasta-extract "$GENOME" "$WRK/Library/${BED}_31bp.bed" -o "$WRK/Library/${BED}.fa"

    java -jar "$SCRIPTMANAGER" figure-generation four-color "$WRK/Library/${BED}.fa" -o "$WRK/Library/${BED}.png"

    python "$RESIZE" -r 600 -c 200 -i "$WRK/Library/${BED}.png" -o "$WRK/Library/${BED}_RESIZE.png"

    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$WRK/Library/${BED}_RESIZE.png" \
        -l "$LLABEL" -m "$MLABEL" -r "$RLABEL" -w 1 -f 20 -o "$WRK/Library/${BED}_RESIZE.svg"

    rm -f "$WRK/Library/${BED}_RESIZE.png" "$WRK/Library/${BED}.fa" "$WRK/Library/${BED}_31bp.bed"
fi

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Heatmap"
    cp $WRK/Library/${BED}_RESIZE.svg "$DIR/Heatmap/"
done

rm $WRK/Library/${BED}_RESIZE.svg



