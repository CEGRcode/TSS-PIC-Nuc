#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory
module load anaconda3
source activate meme

### CHANGE ME
WRK=/Path/to/Title/
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.text
###

# outputs
Reference=$WRK/0X_Bulk_Processing/Reference
# Create output directories if they don't exist
[ -d $Reference ] || mkdir $Reference

bedtools shift -i $WRK/02_TSS_NFR/TSS_4color.bed -g $Genome -p -1 -m 1 > $Reference/TSS_4color_up1.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 6  $Reference/TSS_4color_up1.bed -o $Reference/TSS_4color_up1_6bp.bed
rm $Reference/TSS_4color_up1.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/02_TSS_NFR/nonInr.bed -o $Reference/nonInr_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/02_TSS_NFR/Inr.bed -o $Reference/Inr_500bp.bed

for BEDFILE in $WRK/02_TSS_NFR/Inr_group/TSS_*_all.bed ; do
	BED=$(basename "$BEDFILE" ".bed")
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $BEDFILE -o $Reference/${BED}_500bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/03_core-promoter/FixedTATA_TSS_same.bed -o $Reference/FixedTATA_TSS_same_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/03_core-promoter/FixedTATA_TSS_oppo.bed -o $Reference/FixedTATA_TSS_oppo_1000bp.bed


bedtools shift -i $WRK/03_core-promoter/TATA_TSS-strand_tsssort.bed -g $Genome -p 20 -m -20 > $Reference/TATA_TSS-strand_tsssort_down20bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 60 $Reference/TATA_TSS-strand_tsssort_down20bp.bed -o $Reference/TATA_TSS-strand_tsssort_down20bp_60bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $Reference/TATA_TSS-strand_tsssort_down20bp.bed -o $Reference/TATA_TSS-strand_tsssort_down20bp_100bp.bed
rm $Reference/TATA_TSS-strand_tsssort_down20bp.bed

for file in $WRK/03_core-promoter/TATA_TSS_*_same_*_5prime.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 ${filename}.bed -o $Reference/${filename}_100bp.bed
  rm ${filename}_down30.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $WRK/03_core-promoter/TSS_TATAsame_oppo_noTATA.bed -o $Reference/TSS_TATAsame_oppo_noTATA_200bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/03_core-promoter/TSS_TATAsame_oppo_noTATA.bed -o $Reference/TSS_TATAsame_oppo_noTATA_500bp.bed


for file in   $WRK/03_core-promoter/TSS_noTATA_*.bed  $WRK/03_core-promoter/TSS_TATA_oppo.bed $WRK/03_core-promoter/TSS_TATA_same.bed ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 "$file" -o $Reference/${filename}_1000bp.bed
done

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed $WRK/04_plusoneNucleosome/nophase_ajusted_+1Nuc_Di.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

for file in $WRK/04_plusoneNucleosome/AdjrNuc.bed $WRK/04_plusoneNucleosome/Adj+2Nuc_+1Nuc.bed $WRK/04_plusoneNucleosome/Adj+1Nuc_TSS_all.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 250 $WRK/04_plusoneNucleosome/Adj+1Nuc_TSS_all.bed -o $Reference/Adj+1Nuc_TSS_all_250bp.bed

for file in  $WRK/04_plusoneNucleosome/adj+1Nuc_allphase_TSS.bed ; do
    filename=$(basename "$file" ".bed")
    bedtools shift -i $file -g $Genome -p -100 -m 100 > $Reference/${filename}_up100.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 120 $Reference/${filename}_up100.bed -o $Reference/${filename}_up100_120bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 120 $file -o $Reference/${filename}_1000bp.bed
    rm  $Reference/${filename}_up100.bed
done


for file in  $WRK/04_plusoneNucleosome/adj+1Nuc_phase48.bed $WRK/04_plusoneNucleosome/adj+1Nuc_phase93.bed   ; do
    filename=$(basename "$file" ".bed")
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done


for file in $WRK/05_Call_Motifs/*_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20 $WRK/05_Call_Motifs/${TF}/${TF}_M1.bed -o $Reference/${TF}_M1_20bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/05_Call_Motifs/${TF}/nearestTSS_${TF}_M1_same-oppo.bed -o $Reference/nearestTSS_${TF}_M1_same-oppo_500bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/${TF}/nearestTSS_${TF}_M1_same-oppo.bed -o $Reference/nearestTSS_${TF}_M1_same-oppo_1000bp.bed
    cp $WRK/05_Call_Motifs/${TF}/${TF}_M1.bed  $Reference/${TF}_M1.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${TF}_Occupancy_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $file -o $Reference/${TF}_Occupancy_32bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  $WRK/05_Call_Motifs/WDR5/TSS_same_WDR5_M1.bed -o $Reference/TSS_same_WDR5_M1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/WDR5/WDR5_M1_same_TSS.bed -o $Reference/WDR5_M1_same_TSS_1000bp.bed

for file in $WRK/05_Call_Motifs/SP1_Occupancy_1bp.bed $WRK/05_Call_Motifs/NFYC_Occupancy_1bp.bed $WRK/05_Call_Motifs/GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/${TF}/${TF}_M1_TSS_oppo_up100.bed -o $Reference/${TF}_M1_TSS_oppo_up100_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/${TF}/${TF}_M1_TSS_same_up100.bed -o $Reference/${TF}_M1_TSS_same_up100_1000bp.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/${TF}/${TF}_M1_TSS_same_noTATA.bed -o $Reference/${TF}_M1_TSS_same_noTATA_1000bp.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/${TF}/${TF}_M1_TSS_same_TATA.bed -o $Reference/${TF}_M1_TSS_same_TATA_1000bp.bed 
done

for file in $WRK/05_Call_Motifs/3MOTIF/*with*.bed ; do
    filename=$(basename "$file" ".bed ")
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/YY1/YY1_M1_TSS_oppo_overlap.bed -o $Reference/YY1_M1_TSS_oppo_overlap_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/YY1/YY1_M1_TSS_same_overlap.bed -o $Reference/YY1_M1_TSS_same_overlap_1000bp.bed 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/YY1/YY1_M1_TSS_same_1.bed -o $Reference/YY1_M1_TSS_same_1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/YY1/YY1_M1_TSS_same_2.bed -o $Reference/YY1_M1_TSS_same_2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_Motifs/YY1/YY1_M1_TSS_same_3.bed -o $Reference/YY1_M1_TSS_same_3_1000bp.bed















