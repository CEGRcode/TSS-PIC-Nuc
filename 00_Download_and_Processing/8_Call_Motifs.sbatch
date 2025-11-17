#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/2_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/2_FIMO_Motifs_from_Genome.log.err-%a
#SBATCH --array 1-11

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME

WRK=/Path/to/Title/
CALL_RefPT=$WRK/05_Call_RefPT/
TEMP=$CALL_RefPT/temp_Filter_and_Sort_by_occupancy
WebLogos=$CALL_RefPT/WebLogos
BAMDIR=$WRK/data/BAM

[ -d $WebLogos ] || mkdir $WebLogos
[ -d $TEMP ] || mkdir $TEMP

cd $CALL_RefPT

# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load anaconda3
source activate meme


# Inputs and outputs
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.txt
blacklist=$WRK/data/hg38_files/hg38-blacklist.bed
background=$WRK/data/hg38_files/human_background_model.txt
EXCLUSION=100

ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER
ORIGINAL_CHEXMIX=$WRK/bin/chexmix.v0.52.public.jar
CHEXMIX=$WRK/bin/chexmix.v0.52.public-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_CHEXMIX $CHEXMIX

# Define PWM file path based on SLURM_ARRAY_TASK_ID index
PWM=$(ls PWM/*.meme.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# Parse TF from PWM filename
TF=$(basename "$PWM" "_M1.meme.txt")
[ -d logs ] || mkdir logs
[ -d FIMO ] || mkdir FIMO
[ -d FIMO/$TF ] || mkdir FIMO/$TF
[ -d Chexmix ] || mkdir Chexmix
[ -d Chexmix/$TF ] || mkdir Chexmix/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences
fimo --verbosity 1 --thresh 1.0E-4 --max-strand --oc FIMO/$TF $PWM $GENOME

##get weblogo
ceqlogo -i $PWM -m 1  -o $WebLogos/$TF\_logo.eps -f EPS
ceqlogo -i $PWM -m 1 -r -o $WebLogos/$TF\_logoRC.eps -f EPS

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed FIMO/$TF/fimo.gff -o FIMO/$TF/${TF}_unsorted_noproximity_unfiltered.bed

source deactivate
source activate bioinfo
# Filter out motifs that overlap with blacklist regions

bedtools intersect -v -a FIMO/$TF/${TF}_unsorted_noproximity_unfiltered.bed -b $BLACKLIST | \
awk -v TF="${TF}" '{if ($4 ~ /MEME-1/) print $0 > "FIMO/" TF "/motif1_unsorted_noproximity_temp.bed"}' 

# Filter out motifs that overlap with blacklist regions
awk -v TF="${TF}"  '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6,TF "_M1"}' FIMO/$TF/motif1_unsorted_noproximity_temp.bed | \
awk -v TF="${TF}" '{if ($6 == "+") print $0 > "FIMO/" TF "/motif1_unsorted_noproximity+.bed"; else print $0 > "FIMO/" TF "/motif1_unsorted_noproximity-.bed"}' 
rm FIMO/$TF/motif1_unsorted_noproximity_temp.bed
bedtools intersect -v -a FIMO/$TF/motif1_unsorted_noproximity-.bed -b FIMO/$TF/motif1_unsorted_noproximity+.bed | cat FIMO/$TF/motif1_unsorted_noproximity+.bed - > FIMO/$TF/motif1_unsorted_noproximity.bed
rm  FIMO/$TF/motif1_unsorted_noproximity-.bed FIMO/$TF/motif1_unsorted_noproximity+.bed 

# Apply proximity filter
java -jar $SCRIPTMANAGER peak-analysis filter-bed  FIMO/$TF/motif1_unsorted_noproximity.bed -o FIMO/$TF/motif1_unsorted

# Rename motif-1 file and cleanup
awk -v TF="${TF}"  '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,TF "_M1"}' FIMO/$TF/motif1_unsorted-FILTER.bed > FIMO/$TF/$TF\_motif1_unsorted.bed
rm FIMO/$TF/motif1_unsorted-FILTER.bed
rm FIMO/$TF/motif1_unsorted_noproximity.bed  FIMO/$TF/motif1_unsorted-CLUSTER.bed


# =====Sort by occupancy=====

# Expand WINDOW
for file in FIMO/${TF}/${TF}_*_unsorted.bed ; do
    filename=$(basename "$file" ".bed")
    motif=`basename $file ".bed" | cut -d "_" -f 2`
    BAMFILE=$BAMDIR/K562_${TF}_BX_rep1_hg38.bam
    BAM=`basename $BAMFILE ".bam"`
    WINDOW=100
    awk 'BEGIN{OFS="\t"}{if ($1 !~ /alt|random|chrUn/) print}' $file > $TEMP/${filename}_temp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c ${WINDOW} $TEMP/${filename}_temp.bed -o $TEMP/${filename}_${WINDOW}.bed
    rm $TEMP/${filename}_temp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -s 6 --combined $TEMP/${filename}_${WINDOW}.bed $BAMFILE -M $TEMP/${BAM}_${filename}_${WINDOW}
    java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -c $WINDOW $TEMP/${filename}_${WINDOW}.bed $TEMP/${BAM}_${filename}_${WINDOW}_combined.cdt -o FIMO/$TF/${TF}_${motif}_sorted
    java -jar $SCRIPTMANAGER figure-generation heatmap --black -r 1 -l 2 -p .95 FIMO/$TF/${TF}_${motif}_sorted.cdt -o FIMO/$TF/${TF}_${motif}_sorted.png
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 FIMO/$TF/${TF}_${motif}_sorted.bed -o FIMO/$TF/${TF}_${motif}_sorted_1bp.bed
    rm FIMO/$TF/${TF}_${motif}_sorted.bed
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum FIMO/$TF/${TF}_${motif}_sorted.cdt -o FIMO/$TF/${TF}_${motif}_sorted.out
    tail -n +2 FIMO/$TF/${TF}_${motif}_sorted.out | cut -f 2 | paste  FIMO/$TF/${TF}_${motif}_sorted_1bp.bed - > FIMO/$TF/${TF}_${motif}_Occupancy.bed
    rm FIMO/$TF/${TF}_${motif}_sorted_1bp.bed
    rm FIMO/$TF/${TF}_${motif}_sorted.cdt
done

conda deactivate

rm $SCRIPTMANAGER
## run checkmix of TF 
java -Xmx16G -jar $CHEXMIX --expt $BAMDIR/K562_${TF}_BX_rep1_hg38.bam --format BAM --ctrl "$BAMDIR/K562_IgG_BX_merge_hg38.bam" --threads 4 --geninfo $Genome --back $background --meme1proc --noread2 --seq $GENOME --scalewin 10000 --round 3 --minmodelupdateevents 50 --prlogconf -4 --alphascale 1.0 --betascale 0.05 --epsilonscale 0.2 --exclude $BLACKLIST --nomotifs --out Chexmix/$TF > Chexmix/${TF}.txt



