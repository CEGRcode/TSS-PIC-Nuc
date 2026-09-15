#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory
module load anaconda3
source activate bioinfo

### CHANGE ME
WRK=/Path/to/Title/
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.txt
###

# outputs
Reference=$WRK/0X_Bulk_Processing/Reference
# Create output directories if they don't exist
[ -d $Reference ] || mkdir $Reference
cat $WRK/02_TSS_NFR/Inr_group/TSS_CA_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_TA_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_GA_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_AA_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_CG_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_TG_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_GG_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_AG_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_CC_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_TC_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_GC_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_AC_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_CT_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_TT_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_GT_all.bed $WRK/02_TSS_NFR/Inr_group/TSS_AT_all.bed > $Reference/TSS_4color.bed
bedtools shift -i $Reference/TSS_4color.bed -g $Genome -p -1 -m 1 > $Reference/TSS_4color_up1.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 6  $Reference/TSS_4color_up1.bed -o $Reference/TSS_4color_up1_6bp.bed
rm $Reference/TSS_4color_up1.bed $Reference/TSS_4color.bed
echo "TSS_4color_up1_6bp.bed : $(wc -l < $Reference/TSS_4color_up1_6bp.bed ) sites"
#TSS_4color_up1_6bp.bed : 54917 sites

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/02_TSS_NFR/Inr_group/nonInr.bed -o $Reference/nonInr_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/02_TSS_NFR/Inr_group/Inr.bed -o $Reference/Inr_500bp.bed

for BEDFILE in $WRK/02_TSS_NFR/Inr_group/TSS_*_all.bed ; do
	BED=$(basename "$BEDFILE" ".bed")
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $BEDFILE -o $Reference/${BED}_500bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/03_core-promoter/FixedTATA_TSS_same.bed -o $Reference/FixedTATA_TSS_same_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/03_core-promoter/FixedTATA_TSS_oppo.bed -o $Reference/FixedTATA_TSS_oppo_1000bp.bed
echo "FixedTATA_TSS_oppo_1000bp.bed : $(wc -l < $Reference/FixedTATA_TSS_oppo_1000bp.bed ) sites"
#FixedTATA_TSS_oppo_1000bp.bed : 2071 sites
echo "FixedTATA_TSS_same_1000bp.bed : $(wc -l < $Reference/FixedTATA_TSS_same_1000bp.bed ) sites"
#FixedTATA_TSS_same_1000bp.bed : 5507 sites

bedtools shift -i $WRK/03_core-promoter/TATA_TSS-strand_tsssort.bed -g $Genome -p 20 -m -20 > $Reference/TATA_TSS-strand_tsssort_down20bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 60 $Reference/TATA_TSS-strand_tsssort_down20bp.bed -o $Reference/TATA_TSS-strand_tsssort_down20bp_60bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $Reference/TATA_TSS-strand_tsssort_down20bp.bed -o $Reference/TATA_TSS-strand_tsssort_down20bp_100bp.bed
rm $Reference/TATA_TSS-strand_tsssort_down20bp.bed

for file in $WRK/03_core-promoter/TATA_TSS_*_same_*_5prime.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $file -o $Reference/${filename}_100bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $WRK/03_core-promoter/TSS_TATAsame_oppo_noTATA.bed -o $Reference/TSS_TATAsame_oppo_noTATA_200bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/03_core-promoter/TSS_TATAsame_oppo_noTATA.bed -o $Reference/TSS_TATAsame_oppo_noTATA_500bp.bed

for file in   $WRK/03_core-promoter/TSS_noTATA_*.bed  $WRK/03_core-promoter/TSS_TATA_oppo.bed $WRK/03_core-promoter/TSS_TATA_same.bed ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 "$file" -o $Reference/${filename}_1000bp.bed
done

cat $Reference/TSS_noTATA_1_1000bp.bed  $Reference/TSS_noTATA_2_1000bp.bed  $Reference/TSS_noTATA_3_1000bp.bed $Reference/TSS_noTATA_4_1000bp.bed > $Reference/TSS_noTATA_all_1000bp.bed

echo "TSS_TATA_same_1000bp.bed : $(wc -l < $Reference/TSS_TATA_same_1000bp.bed ) sites"
#TSS_TATA_same_1000bp.bed : 2557 sites
echo "TSS_TATA_oppo_1000bp.bed : $(wc -l < $Reference/TSS_TATA_oppo_1000bp.bed ) sites"
#TSS_TATA_oppo_1000bp.bed : 833 sites

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed $WRK/04_plusoneNucleosome/nophase_ajusted_+1Nuc_Di.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

for file in $WRK/04_plusoneNucleosome/AdjrNuc.bed $WRK/04_plusoneNucleosome/Adj+2Nuc_+1Nuc.bed $WRK/04_plusoneNucleosome/Adj+1Nuc_TSS_all.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

echo "Adj+1Nuc_TSS_all_1000bp.bed : $(wc -l <  $Reference/Adj+1Nuc_TSS_all_1000bp.bed ) sites"
#Adj+1Nuc_TSS_all_1000bp.bed : 45548 sites
echo "TSS_TATA_oppo_1000bp.bed : $(wc -l < $Reference/TSS_TATA_oppo_1000bp.bed ) sites"
#TSS_TATA_oppo_1000bp.bed : 833 sites

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 250 $WRK/04_plusoneNucleosome/Adj+1Nuc_TSS_all.bed -o $Reference/Adj+1Nuc_TSS_all_250bp.bed
### TSS of 7 group of Nuc for function analysis
awk -v file="+1Nuc_" '{ if ( $8 ~ /YRWS/  && $13> 0 && $12> 0 && $11> 0 && $10> 0 && $21 !~ /Divergent/  && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160    )  { print $0 > ( file "YRWS.bed")
     } else if ( $8 ~ /sameYR_lowWS/  && $11> 0 && $10> 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160  ) { print $0 > ( file "sameYR_lowWS.bed")
     } else if ( $8 ~ /sameYR_antiWS/ && $13< 0 && $12< 0 && $11> 0 && $10> 0 && $21 !~ /Divergent/  && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160   ) { print $0 > ( file "sameYR_antiWS.bed")
     } else if ( $8 ~ /antiYR_antiWS/ && $13< 0 && $12< 0 && $11< 0 && $10< 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160  ) { print $0 > ( file "antiYR_antiWS.bed")
     } else if ( $8 ~ /antiYR_sameWS/ && $13> 0 && $12> 0 && $11< 0 && $10< 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160   )  { print $0 > ( file "antiYR_sameWS.bed")
     } else if ( $8 ~ /lessDNAencode/ && $13 == "0" && $12 == "0"  && $11 == "0" && $10 == "0" && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/  && $7 > 0 && $5 >=50 && $5 <=160  )  { print $0 > ( file "lessDNAencode.bed")
     } else if ($8 ~ /lowYR_sameWS/ && $12> 0 && $13> 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/  && $7 > 0 && $5 >=50 && $5 <=160  )  { print $0 > ( file "lowYR_sameWS.bed")
     }
        }' $WRK/04_plusoneNucleosome/Adj+1Nuc_Di_TSS_all_phase.bed

awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_YRWS.bed   >  $Reference/TSS_+1Nuc_YRWS.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_sameYR_lowWS.bed  > $Reference/TSS_+1Nuc_sameYR_lowWS.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_sameYR_antiWS.bed   > $Reference/TSS_+1Nuc_sameYR_antiWS.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_antiYR_antiWS.bed   > $Reference/TSS_+1Nuc_antiYR_antiWS.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_antiYR_sameWS.bed   > $Reference/TSS_+1Nuc_antiYR_sameWS.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_lessDNAencode.bed  > $Reference/TSS_+1Nuc_lessDNAencode.bed
awk '{OFS="\t"; print $14,$15,$16,$17,$18,$19}' +1Nuc_lowYR_sameWS.bed   > $Reference/TSS_+1Nuc_lowYR_sameWS.bed
cd $Reference/
cat TSS_+1Nuc_YRWS.bed TSS_+1Nuc_sameYR_lowWS.bed TSS_+1Nuc_sameYR_antiWS.bed TSS_+1Nuc_antiYR_antiWS.bed TSS_+1Nuc_antiYR_sameWS.bed TSS_+1Nuc_lessDNAencode.bed TSS_+1Nuc_lowYR_sameWS.bed > $Reference/TSS_+1Nuc_all.bed
## 
rm +1Nuc_*.bed
cp $WRK/04_plusoneNucleosome/TSS_all_phase_adj+1Nuc_Di.bed $Reference/
cd ..
## take each group Nuc for all TSS, no filter
awk -v file="+1Nuc_" '{ if ($8 ~ /YRWS/ && $7 > 0 ) { print $0 > ( file "YRWS.bed")
     } else if ($8 ~ /sameYR_lowWS/ && $7 > 0 ) { print $0 > ( file "sameYR_lowWS.bed")
     } else if ($8 ~ /sameYR_antiWS/ && $7 > 0 ) { print $0 > ( file "sameYR_antiWS.bed")
     } else if ($8 ~ /antiYR_antiWS/ && $7 > 0 ) { print $0 > ( file "antiYR_antiWS.bed")
     } else if ($8 ~ /antiYR_sameWS/ && $7 > 0 ) { print $0 > ( file "antiYR_sameWS.bed")
     } else if ($8 ~ /lessDNAencode/ && $7 > 0 ) { print $0 > ( file "lessDNAencode.bed")
     } else if ($8 ~ /lowYR_sameWS/ && $7 > 0 ) { print $0 > ( file "lowYR_sameWS.bed")
     }
        }' $WRK/04_plusoneNucleosome/Adj+1Nuc_Di_TSS_all_phase.bed
mv +1Nuc_*.bed $Reference/

## function analysis for WDR5 motif same strand with TSS
mv +1Nuc_*.bed $Reference/
cp $WRK/05_Call_RefPT/WDR5/TSS_same_WDR5_M1.bed $Reference/
echo "TSS_same_WDR5_M1.bed : $(wc -l < $Reference/TSS_same_WDR5_M1.bed ) sites"
#TSS_same_WDR5_M1.bed : 54 sitess

## for TSS distance filtered Nuc
for file in  $WRK/04_plusoneNucleosome/adj+1Nuc_allphase_TSS.bed ; do
    filename=$(basename "$file" ".bed")
    bedtools shift -i $file -g $Genome -p -100 -m 100 > $Reference/${filename}_up100.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 120 $Reference/${filename}_up100.bed -o $Reference/${filename}_up100_120bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
    rm  $Reference/${filename}_up100.bed
done
echo "adj+1Nuc_allphase_TSS_up100_120bp.bed : $(wc -l < $Reference/adj+1Nuc_allphase_TSS_up100_120bp.bed ) sites"
#adj+1Nuc_allphase_TSS_up100_120bp.bed : 27838 sites

for file in  $WRK/04_plusoneNucleosome/adj+1Nuc_phase48.bed $WRK/04_plusoneNucleosome/adj+1Nuc_phase93.bed   ; do
    filename=$(basename "$file" ".bed")
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

## TFs
for file in $WRK/05_Call_RefPT/*_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20 $WRK/05_Call_RefPT/${TF}/${TF}_M1.bed -o $Reference/${TF}_M1_20bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $WRK/05_Call_RefPT/${TF}/nearestTSS_${TF}_M1_same-oppo.bed -o $Reference/nearestTSS_${TF}_M1_same-oppo_500bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/${TF}/nearestTSS_${TF}_M1_same-oppo.bed -o $Reference/nearestTSS_${TF}_M1_same-oppo_1000bp.bed
    cp $WRK/05_Call_RefPT/${TF}/${TF}_M1.bed  $Reference/${TF}_M1.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${TF}_Occupancy_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $file -o $Reference/${TF}_Occupancy_32bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  $WRK/05_Call_RefPT/WDR5/TSS_same_WDR5_M1.bed -o $Reference/TSS_same_WDR5_M1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/WDR5/WDR5_M1_same_TSS.bed -o $Reference/WDR5_M1_same_TSS_1000bp.bed

for file in $WRK/05_Call_RefPT/SP1_Occupancy_1bp.bed $WRK/05_Call_RefPT/NFYC_Occupancy_1bp.bed $WRK/05_Call_RefPT/GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/${TF}/${TF}_M1_TSS_oppo_up100.bed -o $Reference/${TF}_M1_TSS_oppo_up100_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/${TF}/${TF}_M1_TSS_same_up100.bed -o $Reference/${TF}_M1_TSS_same_up100_1000bp.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/${TF}/${TF}_M1_TSS_same_noTATA.bed -o $Reference/${TF}_M1_TSS_same_noTATA_1000bp.bed 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/${TF}/${TF}_M1_TSS_same_TATA.bed -o $Reference/${TF}_M1_TSS_same_TATA_1000bp.bed 
done

for file in $WRK/05_Call_RefPT/3MOTIF/*with*.bed ; do
    filename=$(basename "$file" ".bed")
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $Reference/${filename}_1000bp.bed
done

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/YY1/YY1_M1_TSS_oppo_overlap.bed -o $Reference/YY1_M1_TSS_oppo_overlap_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/YY1/YY1_M1_TSS_same_overlap.bed -o $Reference/YY1_M1_TSS_same_overlap_1000bp.bed 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/YY1/YY1_M1_TSS_same_1.bed -o $Reference/YY1_M1_TSS_same_1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/YY1/YY1_M1_TSS_same_2.bed -o $Reference/YY1_M1_TSS_same_2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/05_Call_RefPT/YY1/YY1_M1_TSS_same_3.bed -o $Reference/YY1_M1_TSS_same_3_1000bp.bed



