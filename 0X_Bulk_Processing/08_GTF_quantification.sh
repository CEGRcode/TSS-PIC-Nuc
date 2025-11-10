module load anaconda3
source activate bioinfo
# Script for F2a and E5

### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
[ -d $WRK/Library/F2a ] || mkdir -p $WRK/Library/F2a
[ -d $WRK/Library/E5 ] || mkdir -p $WRK/Library/E5

BAMDIR=$WRK/data/BAM
cd $WRK/Library/F2a

### calculate the occupancy upstream 100bp vs downstream 100bp


for file in $BAMDIR/*_TAF4B_*.bam $BAMDIR/*_TAF4_*.bam   $BAMDIR/*_MED12_*.bam $BAMDIR/*_GTF2B_*.bam  $BAMDIR/*_GTF2A1_*.bam $BAMDIR/*_TBP_*.bam  $BAMDIR/*_TAF3_*.bam $BAMDIR/*_TAF1_*.bam $BAMDIR/*_TAF9_*.bam $BAMDIR/*_TAF6_*.bam $BAMDIR/*_TAF11_*.bam  $BAMDIR/*_TAF13_*.bam $BAMDIR/*_TAF10_*.bam  ; do
  BAM=$(basename "$file" ".bam")
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 401-500  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_U.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 401-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_U.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_D.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_D.out

  paste  ${BAM}_TSS_TATA_same_U.out ${BAM}_TSS_TATA_same_D.out  ${BAM}_TSS_noTATA_4_U.out  ${BAM}_TSS_noTATA_4_D.out | \
  awk '{OFS="\t"} {print ($2-$1)/($1+$2),($4-$3)/($3+$4)}' > ${BAM}_TSS_same_noTATA_UD.out
  rm  ${BAM}_TSS_TATA_same_U.out ${BAM}_TSS_TATA_same_D.out  ${BAM}_TSS_noTATA_4_U.out  ${BAM}_TSS_noTATA_4_D.out
done


for file in $BAMDIR/*_GTF2A1_*.bam $BAMDIR/*_TBP_*.bam ; do
    BAM=$(basename "$file" ".bam")
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done

for file in  $BAMDIR/*_Pol2_*.bam   ; do
  BAM=$(basename "$file" ".bam")
  
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done

for file in $BAMDIR/*_TAF4B_*.bam  ; do
  BAM=$(basename "$file" ".bam")
  
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done

cd $WRK/Library/E5
    
for file in $BAMDIR/*_GTF2A2_*.bam  $BAMDIR/*_DR1_*.bam  ; do
    BAM=$(basename "$file" ".bam")
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done


for file in  $BAMDIR/*_NELFA_*.bam  $BAMDIR/*_NELFE_*.bam   ; do
  BAM=$(basename "$file" ".bam")
  
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done

for file in $BAMDIR/*_TAF4_*.bam   $BAMDIR/*_MED12_*.bam ; do
  BAM=$(basename "$file" ".bam")
  
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done


for file in $BAMDIR/*_GTF2B_*.bam ; do
  BAM=$(basename "$file" ".bam")
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done


for file in  $BAMDIR/*_TAF1_*.bam  $BAMDIR/*_TAF3_*.bam $BAMDIR/*_TAF11_*.bam $BAMDIR/*_TAF13_*.bam  $BAMDIR/*_TAF10_*.bam  $BAMDIR/*_TAF9_*.bam $BAMDIR/*_TAF6_*.bam  $BAMDIR/*_TAF9B_*.bam  $BAMDIR/*_ERCC3_*.bam   $BAMDIR/*_GTF2H1_*.bam ; do
 BAM=$(basename "$file" ".bam")
  
  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_F.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_F.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/${BAM}_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_1_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_2_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_3_B.out
  tail -n +2 Composites/${BAM}_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_noTATA_4_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_same_B.out
  tail -n +2 Composites/${BAM}_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${BAM}_TSS_TATA_oppo_B.out

  paste  ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${BAM}_TSS_same_oppo_noTATA.out
  rm ${BAM}_TSS_TATA_same_F.out ${BAM}_TSS_TATA_same_B.out ${BAM}_TSS_TATA_oppo_F.out ${BAM}_TSS_TATA_oppo_B.out ${BAM}_TSS_noTATA_1_F.out ${BAM}_TSS_noTATA_1_B.out ${BAM}_TSS_noTATA_2_F.out ${BAM}_TSS_noTATA_2_B.out ${BAM}_TSS_noTATA_3_F.out ${BAM}_TSS_noTATA_3_B.out  ${BAM}_TSS_noTATA_4_F.out  ${BAM}_TSS_noTATA_4_B.out
done













