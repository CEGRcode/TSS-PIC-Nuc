#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/06_Conservation-phylo_pileups.log.out-%a
#SBATCH -e logs/06_Conservation-phylo_pileups.log.err-%a
#SBATCH --array 1-13

### CHANGE ME
WRK=/Path/to/Title/
METADATA=$WRK/0X_Bulk_Processing/Conservation-SNP.txt
THREADS=4
###

# Dependencies
# - java
# - perl
# - python
# - pyBigWig

set -exo
module load anaconda
source activate bioinfo

# Fill in placeholder constants with your directories
CDIR=$WRK/data/Conservation-SNP

# Script shortcuts
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
PILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT.py
# Set up output directories
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
cd $WRK/0X_Bulk_Processing/

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE ".bed"`
[[ -f $BEDFILE ]] || { echo "ERROR: BED file not found: $BEDFILE"; exit 1; }

# Determine output figure panel(s)
PANEL=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
IFS=';' read -ra PANELS <<< "$PANEL"

# Pileup conservation scores
CONSERVATION=$CDIR/hg38.phyloP30way.bw
CONS=phyloP30way

echo $BED x $CONS

# Pileup BigWig
python $PILEUPBW -i $CONSERVATION -r $BEDFILE -o $WRK/Library/${CONS}_$BED.cdt
# Make composite
perl $COMPOSITE $WRK/Library/${CONS}_$BED.cdt $WRK/Library/${CONS}_$BED.out

# Pileup SNPs
DBSNP=$CDIR/dbSnp153_snv.bw
SNP=SNP

echo $BED x $SNP

# Pileup BigWig
python $PILEUPBW -i $DBSNP -r $BEDFILE -o $WRK/Library/${SNP}_$BED.cdt
# Make composite
perl $COMPOSITE $WRK/Library/${SNP}_$BED.cdt $WRK/Library/${SNP}_$BED.out

# Loop over each panel name
for P in "${PANELS[@]}"; do
    DIR=$WRK/Library/$P
    mkdir -p "$DIR/Composites"
    cp $WRK/Library/${CONS}_$BED.out "$DIR/Composites/"
    cp $WRK/Library/${SNP}_$BED.out "$DIR/Composites/"
    mkdir -p "$DIR/CDT"
    cp $WRK/Library/${CONS}_$BED.cdt "$DIR/CDT/"
    cp $WRK/Library/${SNP}_$BED.cdt "$DIR/CDT/"
done
# Remove scaled matrix from top-level Library folder
rm $WRK/Library/${CONS}_$BED.cdt
rm $WRK/Library/${SNP}_$BED.cdt
rm $WRK/Library/${CONS}_$BED.out
rm $WRK/Library/${SNP}_$BED.out