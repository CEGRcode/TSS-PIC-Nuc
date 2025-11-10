#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/05_Dinucleotide_heatmap.log.out-%a
#SBATCH -e logs/05_Dinucleotide_heatmap.log.err-%a
#SBATCH --array 1-4

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/X_Bulk_Processing/Dinucleotide_heatmap.txt
THREADS=4
###
# Dependencies
# - java
# - python

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
cd $WRK/X_Bulk_Processing/

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }
# Determine the dinucleotide
DI=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`

# Determine heatmap contrast threshold for the current job array index
CONTRAST_THRESH=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Determine color
COLOR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
IFS=';' read -ra PANELS <<< "$PANEL"

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} x ${DI}"
BASE=${DI}_${BED}

# ===============================================================================================================================

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

# Set window labels
RLABEL="R"
LLABEL="L"
MLABLE="M"
WINDOW=`echo $BED | rev | cut -d"_" -f1 | rev`
[ $WINDOW == "500bp" ] && LLABEL="-250" && RLABEL="+250" && MLABEL="0"

# Parse heatmap params
HM_PARAMS="-p 0.95"
[[ $CONTRAST_THRESH =~ "NaN" ]] || HM_PARAMS="-a $CONTRAST_THRESH"

# ===============================================================================================================================


echo "Run Dinucleotode scan, make heatmap"

java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $BEDFILE -o $WRK/Library/${BED}.fa
conda deactivate
python $MOTIFSCAN -i  $WRK/Library/${BED}.fa -m ${DI} -o $WRK/Library/${DI}_${BED}
source activate bioinfo
rm $WRK/Library/${DI}_${BED}_anti.cdt

java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS -c $COLOR  $WRK/Library/${DI}_${BED}_sense.cdt -o $WRK/Library/${DI}_${BED}_sense.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $WRK/Library/${DI}_${BED}_sense.png \
	-l $LLABEL -m $MLABEL -r $RLABEL -w 1 -f 20 \
	#-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $WRK/Library/${DI}_${BED}_label.svg
rm $WRK/Library/${DI}_${BED}_sense.cdt
rm $WRK/Library/${DI}_${BED}_sense.png
rm $WRK/Library/${BED}.fa

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Heatmap"
    cp $WRK/Library/${DI}_${BED}_label.svg "$DIR/Heatmap/"
done

rm $WRK/Library/${DI}_${BED}_label.svg 





