module load anaconda3
source activate bioinfo


### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
[ -d $WRK/Library/F2a ] || mkdir -p $WRK/Library/F2a
[ -d $WRK/Library/F2b ] || mkdir -p $WRK/Library/F2b
[ -d $WRK/Library/E5a ] || mkdir -p $WRK/Library/E5a
[ -d $WRK/Library/F5d ] || mkdir -p $WRK/Library/F5d

BAMDIR=$WRK/data/BAM/
cd $WRK/Library/F2b

### calculate the occupancy upstream 100bp vs downstream 100bp

for file in $BAMDIR/K562_TAF4B_*.bam $BAMDIR/K562_TAF4_*.bam   $BAMDIR/K562_MED12_*.bam $BAMDIR/K562_GTF2B_*.bam  $BAMDIR/K562_GTF2A1_*.bam $BAMDIR/K562_TBP_*.bam  $BAMDIR/K562_TAF3_*.bam $BAMDIR/K562_TAF1_*.bam $BAMDIR/K562_TAF9_*.bam $BAMDIR/K562_TAF6_*.bam $BAMDIR/K562_TAF11_*.bam  $BAMDIR/K562_TAF13_*.bam $BAMDIR/K562_TAF10_*.bam  ; do
  TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 401-500  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_U.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 401-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_U.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_D.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_D.out

  paste ${TF}_TSS_TATA_same_U.out ${TF}_TSS_TATA_same_D.out  ${TF}_TSS_noTATA_4_U.out  ${TF}_TSS_noTATA_4_D.out | \
  awk '{OFS="\t"} {print ($2-$1)/($1+$2),($4-$3)/($3+$4)}' > $BAM\_TSS_same_noTATA_UD.out
  rm  ${TF}_TSS_TATA_same_U.out ${TF}_TSS_TATA_same_D.out  ${TF}_TSS_noTATA_4_U.out  ${TF}_TSS_noTATA_4_D.out
done

cd $WRK/Library/F2a

for file in $BAMDIR/K562_GTF2A1_*.bam $BAMDIR/K562_TBP_*.bam ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done

for file in $BAMDIR/K562_PolII_*.bam  ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done

for file in  $BAMDIR/K562_TAF4B_*.bam   ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done

cd $WRK/Library/E5a
    
for file in $BAMDIR/K562_GTF2A2_*.bam  $BAMDIR/K562_DR1_*.bam  ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 441-500 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-560 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done


for file in  $BAMDIR/K562_NELFA_*.bam  $BAMDIR/K562_NELFE_*.bam   ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 501-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-700 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done

for file in $BAMDIR/K562_TAF4_*.bam   $BAMDIR/K562_MED12_*.bam ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 351-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-850 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done


for file in $BAMDIR/K562_GTF2B_*.bam ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 451-550 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 551-650 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done


for file in  $BAMDIR/K562_TAF1_*.bam  $BAMDIR/K562_TAF3_*.bam $BAMDIR/K562_TAF11_*.bam $BAMDIR/K562_TAF13_*.bam  $BAMDIR/K562_TAF10_*.bam  $BAMDIR/K562_TAF9_*.bam $BAMDIR/K562_TAF6_*.bam  $BAMDIR/K562_TAF9B_*.bam  $BAMDIR/K562_ERCC3_*.bam   $BAMDIR/K562_GTF2H1_*.bam ; do
    TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_F.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 401-600 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_F.out

  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_1_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_1_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_2_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_2_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_3_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_3_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_noTATA_4_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_noTATA_4_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_same_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_same_B.out
  tail -n +2 Composites/K562_${TF}_BX_rep1_hg38_TSS_TATA_oppo_1000bp_5read1_Normalized.out | cut -f 601-800 | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_TSS_TATA_oppo_B.out

  paste  ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out | \
  awk '{OFS="\t"} {print $1-$2,$3-$4,$5-$6,$7-$8,$9-$10,$11-$12}' > ${TF}_TSS_same_oppo_noTATA.out
  rm ${TF}_TSS_TATA_same_F.out ${TF}_TSS_TATA_same_B.out ${TF}_TSS_TATA_oppo_F.out ${TF}_TSS_TATA_oppo_B.out ${TF}_TSS_noTATA_1_F.out ${TF}_TSS_noTATA_1_B.out ${TF}_TSS_noTATA_2_F.out ${TF}_TSS_noTATA_2_B.out ${TF}_TSS_noTATA_3_F.out ${TF}_TSS_noTATA_3_B.out  ${TF}_TSS_noTATA_4_F.out  ${TF}_TSS_noTATA_4_B.out
done

cd $WRK/Library/F2d
### calculate nucleosome engagement Triptolide treatment vs DMSO treatment

for file in  $BAMDIR/TriptolideK562_TFIIB_*.bam $BAMDIR/TriptolideK562_PolII_*.bam $BAMDIR/TriptolideK562_TBP_*.bam   $BAMDIR/TriptolideK562_GTF2A1_*.bam  ; do
   TF=`basename $file ".bam" | cut -d "_" -f 2`
  tail -n +2 Composites/TriptolideK562_${TF}_BX_rep1_hg38_Adj+1Nuc_TSS_all_1000bp_5read1_Normalized.out | cut -f 426-575  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_Adj+1Nuc_Trip.out
  tail -n +2 Composites/TriptolideK562_${TF}_BX_rep1_hg38_Adj+1Nuc_TSS_all_1000bp_5read1_Normalized.out | cut -f 251-750  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_all_Trip.out
  tail -n +2 Composites/DMSOK562_${TF}_BX_rep1_hg38_Adj+1Nuc_TSS_all_1000bp_5read1_Normalized.out | cut -f 426-575  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_Adj+1Nuc_DMSO.out
  tail -n +2 Composites/DMSOK562_${TF}_BX_rep1_hg38_Adj+1Nuc_TSS_all_1000bp_5read1_Normalized.out | cut -f 251-750  | awk '{
    for (i = 1; i <= NF; i++) sum += $i } END { print sum }' > ${TF}_all_DMSO.out
  paste ${TF}_Adj+1Nuc_Trip.out ${TF}_all_Trip.out ${TF}_Adj+1Nuc_DMSO.out ${TF}_all_DMSO.out | \
  awk '{OFS="\t"} {print ($1)/($2),($3)/($4)}' > ${TF}_Tris-DMSO.out
  rm  ${TF}_Adj+1Nuc_Trip.out ${TF}_all_Trip.out ${TF}_Adj+1Nuc_DMSO.out ${TF}_all_DMSO.out 
done














