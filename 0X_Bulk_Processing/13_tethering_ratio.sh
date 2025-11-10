module load anaconda3
source activate bioinfo
# Script for F2a and E5

### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/F5h ] || mkdir -p $WRK/Library/F5h
[ -d $WRK/Library/F5i ] || mkdir -p $WRK/Library/F5i


cd $WRK/Library/F5h




for file in NFYC_Occupancy_1bp.bed  SP1_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    filename=$(basename "$file" "_1bp.bed")
    TFBAM=$BAMDIR/K562_${TF}_BX_rep1_hg38.bam
    BAM=K562_${TF}_BX_rep1_hg38.bam
    FACTOR=`grep 'Scaling factor' $NormDIR/K562_${TF}_BX_rep1_hg38_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`
    bedtools sort -i ${TF}/TSS_same_${TF}_M1_TATA.bed > TSS_same_${TF}_M1_TATA_sort.bed
    bedtools intersect -u -a ../Nuc1/TSS_all_phase_adj+1Nuc_Di.bed -b  TSS_same_${TF}_M1_TATA_sort.bed | \
    awk '{OFS="\t"} { print $17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ${TF}/TATA_TSS_same_${TF}_M1_sort.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_TATA_sort.bed > ${TF}/${TF}_M1_TSS_same_TATA_sort.bed
    rm TSS_same_${TF}_M1_TATA_sort.bed
    bedtools sort -i ${TF}/TSS_same_${TF}_M1_noTATA.bed > TSS_same_${TF}_M1_noTATA_sort.bed
    bedtools intersect -u -a ../Nuc1/TSS_all_phase_adj+1Nuc_Di.bed -b  TSS_same_${TF}_M1_noTATA_sort.bed | \
    awk '{OFS="\t"} { print $17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ${TF}/noTATA_TSS_same_${TF}_M1_sort.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_noTATA_sort.bed > ${TF}/${TF}_M1_TSS_same_noTATA_sort.bed
    rm TSS_same_${TF}_M1_noTATA_sort.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 40 ${TF}/TATA_TSS_same_${TF}_M1_sort.bed -o TATA_TSS_same_${TF}_M1_sort_40bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 40 ${TF}/noTATA_TSS_same_${TF}_M1_sort.bed -o noTATA_TSS_same_${TF}_M1_sort_40bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 40 ${TF}/${TF}_M1_TSS_same_TATA_sort.bed -o ${TF}_M1_TSS_same_TATA_sort_40bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 40 ${TF}/${TF}_M1_TSS_same_noTATA_sort.bed -o ${TF}_M1_TSS_same_noTATA_sort_40bp.bed
done

for file in NFYC_Occupancy_1bp.bed  SP1_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed   ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)

    for file in  $BAMDIR/*_TBP_*.bam  $BAMDIR/*_GTF2A1_*.bam  $BAMDIR/*_GTF2B_*.bam $BAMDIR/*_Taf4_*.bam ; do
     BAM=$(basename "$file" ".bam")
     FACTOR=`grep 'Scaling factor' $NormDIR/${BAM}_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`
     java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  TATA_TSS_same_${TF}_M1_sort_40bp.bed "$file" --cpu 4 -s 6 --combined -M ${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1
     java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  noTATA_TSS_same_${TF}_M1_sort_40bp.bed  "$file" --cpu 4 -s 6 --combined  -M ${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1
     java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  ${TF}_M1_TSS_same_TATA_sort_40bp.bed "$file" --cpu 4 -s 6 --combined -M ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1
     java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${TF}_M1_TSS_same_noTATA_sort_40bp.bed  "$file" --cpu 4 -s 6 --combined  -M ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1

     java -jar "$SCRIPTMANAGER" read-analysis scale-matrix ${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1_combined.cdt -s "$FACTOR" -o  ${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1_combined_Normalized.cdt
     java -jar "$SCRIPTMANAGER" read-analysis scale-matrix ${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1_combined.cdt -s "$FACTOR" -o  ${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_Normalized.cdt
     java -jar "$SCRIPTMANAGER" read-analysis scale-matrix ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined.cdt -s "$FACTOR" -o  ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined_Normalized.cdt
     java -jar "$SCRIPTMANAGER" read-analysis scale-matrix ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined.cdt -s "$FACTOR" -o  ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_Normalized.cdt

     rm ${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1_combined.cdt
     rm ${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1_combined.cdt
     rm ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined.cdt
     rm ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined.cdt

     java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1_combined_Normalized.cdt -o ${TF}/${BAM}_TATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out 
     java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_Normalized.cdt -o ${TF}/${BAM}_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out 
     java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined_Normalized.cdt -o ${TF}/${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined_SCORES.out 
     java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_Normalized.cdt -o ${TF}/${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_SCORES.out 
     rm ${BAM}_TATA_TSS_same_${TF}_M1_20bp_read1_combined_Normalized.cdt
     rm ${BAM}_noTATA_TSS_same_${TF}_M1_20bp_read1_combined_Normalized.cdt
     rm ${BAM}_${TF}_M1_TSS_same_TATA_40bp_read1_combined_Normalized.cdt
     rm ${BAM}_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_Normalized.cdt
    done
done




for file in NFYC_Occupancy_1bp.bed  SP1_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    paste ${TF}/*_TBP_*_TATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out  ${TF}/*_TBP_*_${TF}_M1_TSS_same_TATA_40bp_read1_combined_SCORES.out ${TF}/*_GTF2A1_*_TATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_GTF2A1_*_${TF}_M1_TSS_same_TATA_40bp_read1_combined_SCORES.out ${TF}/*_GTF2B_*_TATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_GTF2B_*_${TF}_M1_TSS_same_TATA_40bp_read1_combined_SCORES.out ${TF}/*_Taf4_*_TATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_Taf4_*_${TF}_M1_TSS_same_TATA_40bp_read1_combined_SCORES.out | \
    tail -n +2  | cut -f 2,4,6,8,10,12,14,16 | paste ${TF}/${TF}_M1_TSS_same_TATA_sort.bed - > ${TF}/${TF}_M1_TSS_same_TATA_sort_score.bed

    paste ${TF}/*_TBP_*_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out  ${TF}/*_TBP_*_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_SCORES.out ${TF}/*_GTF2A1_*_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_GTF2A1_*_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_SCORES.out ${TF}/*_GTF2B_*_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_GTF2B_*_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_SCORES.out  ${TF}/*_Taf4_*_noTATA_TSS_same_${TF}_M1_40bp_read1_combined_SCORES.out ${TF}/*_Taf4_*_${TF}_M1_TSS_same_noTATA_40bp_read1_combined_SCORES.out | \
    tail -n +2  | cut -f  2,4,6,8,10,12,14,16  | paste ${TF}/${TF}_M1_TSS_same_noTATA_sort.bed - > ${TF}/${TF}_M1_TSS_same_noTATA_sort_score.bed

    rm *_40bp.bed
    rm *_20bp.bed
done
