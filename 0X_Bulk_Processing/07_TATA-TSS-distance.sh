
module load anaconda3
source activate bioinfo


### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
##DATA
BAMDIR=$WRK/data/BAM/
NormDIR="$WRK/data/NormalizationFactors"
TBP="$BAMDIR/K562_TBP_BX_rep1_hg38.bam"
GTF2A="$BAMDIR/K562_GTF2A1_BX_rep1_hg38.bam"
TBPFACTOR=`grep 'Scaling factor' $NormDIR/K562_TBP_BX_rep1_hg38_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`
GTF2AFACTOR=`grep 'Scaling factor' $NormDIR/K562_GTF2A1_BX_rep1_hg38_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`

## Determine output
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
[ -d $WRK/Library/F1C ] || mkdir -p $WRK/Library/F1C

cd $WRK/Library/F1C


## GET TATA box distance and TBP/TFIIA occupancy info

for file in $WRK/03_core-promoter/TATA_TSS_same_5prime.bed; do
    filename=$(basename "$file" .bed)
    awk -v filename="$filename" '{
        distance = $22 + 0  # Ensure distance is treated as a number
        if (distance >= -31 && distance <= -21) {
            print $0 >> (filename"_"distance".bed")
        }
    }' "$file"
done

wc -l TATA_TSS_same_5prime_*.bed

     310 TATA_TSS_same_5prime_-21.bed
     306 TATA_TSS_same_5prime_-22.bed
     411 TATA_TSS_same_5prime_-23.bed
     484 TATA_TSS_same_5prime_-24.bed
     711 TATA_TSS_same_5prime_-25.bed
     904 TATA_TSS_same_5prime_-26.bed
     746 TATA_TSS_same_5prime_-27.bed
     569 TATA_TSS_same_5prime_-28.bed
     342 TATA_TSS_same_5prime_-29.bed
     215 TATA_TSS_same_5prime_-30.bed
     161 TATA_TSS_same_5prime_-31.bed

for file in TATA_TSS_same_5prime_*.bed ; do
  BED=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 60 $file -o ${BED}_60bp.bed
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_60bp.bed $TBP --cpu 4 -1 --combined -o TBP_${BED}_read1_combined_60bp.out
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix -r 1 -l 1 TBP_${BED}_read1_combined_60bp.out -s "$TBPFACTOR" -o TBP_${BED}_read1_combined_60bp_Normalized.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_60bp.bed $GTF2A --cpu 4 -1 --combined -o GTF2A_${BED}_read1_combined_60bp.out
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix -r 1 -l 1 GTF2A_${BED}_read1_combined_60bp.out -s "$GTF2AFACTOR" -o GTF2A_${BED}_read1_combined_60bp_Normalized.out
  rm TBP_${BED}_read1_combined_60bp.out  
  rm GTF2A_${BED}_read1_combined_60bp.out
  rm ${BED}.bed ${BED}_60bp.bed
done

#!/bin/bash

# Output file
out_file="TBP_occupancy_TSS-TATAdistance.out"
> "$out_file"  # clear file if it exists

# Loop over the files you want to process
for i in {-21..-31}
do
    file="TBP_TATA_TSS_same_5prime_${i}_read1_combined_60bp_Normalized.out"
    if [[ -f "$file" ]]; then
        # Run awk to sum the 2nd line
        sum=$(awk 'NR == 2 {s=0; for(i=1; i<=NF; i++) s += $i} END {print s}' "$file")
        echo "$file $sum" >> "$out_file"
    else
        echo "$file not found" >> "$out_file"
    fi
done

rm TBP_TATA_TSS_same_5prime_*_read1_combined_60bp_Normalized.out

#!/bin/bash

# Output file
out_file="GTF2A_occupancy_TSS-TATAdistance.out"
> "$out_file"  # clear file if it exists

# Loop over the files you want to process
for i in {-21..-31}
do
    file="GTF2A_TATA_TSS_same_5prime_${i}_read1_combined_60bp_Normalized.out"
    if [[ -f "$file" ]]; then
        # Run awk to sum the 2nd line
        sum=$(awk 'NR == 2 {s=0; for(i=1; i<=NF; i++) s += $i} END {print s}' "$file")
        echo "$file $sum" >> "$out_file"
    else
        echo "$file not found" >> "$out_file"
    fi
done

rm GTF2A_TATA_TSS_same_5prime_*_read1_combined_60bp_Normalized.out


cat  $WRK/03_core-promoter/TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed $WRK/03_core-promoter/TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed | awk '{ if (($22 >= -56) && ($22 <= -4) && ($8 ~ /TATAsame/)) print $0 > "TSS_first_second_sameTATA.bed"  }'


#!/bin/bash

# Input BED file
input="TSS_first_second_sameTATA.bed"
filename=$(basename "$input" .bed)

# 1️⃣ Split into separate distance files (-46 to -6)
awk -v filename="$filename" '{
    distance = $22 + 0  # Ensure numeric
    if (distance >= -46 && distance <= -6) {
        print $0 >> (filename"_"distance".bed")
    }
}' "$input"

# 2️⃣ Create output summary file
out="TSS-TATA-distance-count.out"
> "$out"

# 3️⃣ Count lines in each generated BED file and save to output
for f in ${filename}_-*.bed; do
    if [[ -f "$f" ]]; then
        count=$(wc -l < "$f")
        echo -e "$f\t$count" >> "$out"
    fi
done

# 4️⃣ Sort by distance (optional, for tidy output)
sort -t_ -k4,4n "$out" -o "$out"

echo "✅ Report written to $out"

rm TSS_first_second_sameTATA_*.bed 

sort -k22,22n TSS_first_second_sameTATA.bed | awk '{if ( $22 >= -31 && $22 <= -21 ) print $0 > "TSS_fixTATA.bed" ; else if ( $22 >= -46 && $22 <= -31 ) print $0 > "TSS_nonfixTATA_1.bed" ; else if ( $22 >= -21 && $22 <= -6 ) print $0 > "TSS_nonfixTATA_2.bed"   }'
sort -k22,22n TSS_first_second_sameTATA.bed | awk '{OFS="\t"} {print $16,$17,$18,$19,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$22 }'  >  TATA_TSS_first_second_same_sort.bed

awk -v filename="TATA_TSS_first_second_same" '{
    distance = $22 + 0  # Ensure numeric
    if (distance >= -46 && distance <= -6) {
        print $0 >> (filename"_"distance".bed")
    }
}' TATA_TSS_first_second_same_sort.bed



for file in TATA_TSS_first_second_same_*.bed ; do
  BED=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $file -o ${BED}_200bp.bed
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_200bp.bed $TBP --cpu 4 -1 -o TBP_${BED}_read1_200bp.out
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix -r 1 -l 1 TBP_${BED}_read1_200bp.out -s "$TBPFACTOR" -o TBP_${BED}_read1_200bp_Normalized.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_200bp.bed $GTF2A --cpu 4 -1 -o GTF2A_${BED}_read1_200bp.out
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix -r 1 -l 1 GTF2A_${BED}_read1_200bp.out -s "$GTF2AFACTOR" -o GTF2A_${BED}_read1_200bp_Normalized.out
  rm GTF2A_${BED}_read1_200bp.out
  rm TBP_${BED}_read1_200bp.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_200bp.bed $CoPRO --cpu 4 -2 -o CoPRO_${BED}_read2_200bp.out
  rm ${BED}_200bp.bed
done


# Pileup (read 1)

for file in TATA_TSS_first_second_same_sort.bed ; do
  BED=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $file -o ${BED}_100bp.bed
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_100bp.bed $TBP --cpu 4 -1 -M TBP_${BED}_read1_100bp
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix  TBP_${BED}_read1_100bp_sense.cdt -s "$TBPFACTOR" -o TBP_${BED}_read1_100bp_sense_Normalized.cdt
  java -jar "$SCRIPTMANAGER" read-analysis scale-matrix TBP_${BED}_read1_100bp_anti.cdt -s "$TBPFACTOR" -o TBP_${BED}_read1_100bp_anti_Normalized.cdt
  java -jar  $SCRIPTMANAGER figure-generation heatmap  --blue TBP_${BED}_read1_100bp_sense_Normalized.cdt -o  TBP_${BED}_read1_100bp_sense_treeview.png
  java -jar  $SCRIPTMANAGER figure-generation heatmap  --red TBP_${BED}_read1_100bp_anti_Normalized.cdt -o  TBP_${BED}_read1_100bp_anti_treeview.png
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${BED}_100bp.bed $CoPRO --cpu 4 -2 -M CoPRO_${BED}_read2_100bp
  java -jar  $SCRIPTMANAGER figure-generation heatmap  --blue CoPRO_${BED}_read2_100bp_sense.cdt -o  CoPRO_${BED}_read2_100bp_sense.png
  java -jar  $SCRIPTMANAGER figure-generation heatmap  --red CoPRO_${BED}_read2_100bp_anti.cdt -o  CoPRO_${BED}_read2_100bp_anti.png
  java -jar $SCRIPTMANAGER figure-generation merge-heatmap TBP_${BED}_read1_100bp_sense_treeview.png TBP_${BED}_read1_100bp_anti_treeview.png -o TBP_${BED}_read1_100bp_merge.png
  java -jar  $SCRIPTMANAGER figure-generation merge-heatmap CoPRO_${BED}_read2_100bp_sense.png CoPRO_${BED}_read2_100bp_anti.png -o CoPRO_${BED}_read2_100bp_merge.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap  CoPRO_${BED}_read2_100bp_merge.png -l -50 -m 0 -r 50 -w 1 -f 20 -o CoPRO_${BED}_read2_100bp.svg
  java -jar $SCRIPTMANAGER figure-generation label-heatmap  TBP_${BED}_read1_100bp_merge.png -l -50 -m 0 -r 50 -w 1 -f 20 -o TBP_${BED}_read1_100bp.svg
  #rm TBP_${BED}_read1_100bp_sense_treeview.png TBP_${BED}_read1_100bp_anti_treeview.png
  #rm CoPRO_${BED}_read2_100bp_sense.png CoPRO_${BED}_read2_100bp_anti.png
  #rm  TBP_${BED}_read1_100bp_*.cdt CoPRO_${BED}_read2_100bp_*.cdt
  rm ${BED}_100bp.bed
done






