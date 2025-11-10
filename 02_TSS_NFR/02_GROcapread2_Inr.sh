module load anaconda3
source activate bioinfo
# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system
WRK=/Path/to/Title/
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
#Genome
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.text
Nuc="$WRK/Annotation/BNase-Nucleosomes.bed"
## data reference
CpGBED="$WRK/Annotation/CpGIslands.bed"
TSSBED="$WRK/Annotation/hg38_TSS-ALL_GENCODEV47_SORT.bed"
BLACKLIST=$WRK/hg38_files/hg38-blacklist.bed
tRNA="$WRK/Annotation/gencode.v47.tRNAs.bed"
NoncodingRNABED="$WRK/Annotation/gencode_v47.long_noncoding_RNAs_TSS.bed"
Promoter="$WRK/Annotation/EncodeBroadHmmK562_Promoter.bed"
Enhancer="$WRK/Annotation/EncodeBroadHmmK562_Enhancer.bed"
Insulator="$WRK/Annotation/EncodeBroadHmmK562_Insulator.bed"
Repressive="$WRK/Annotation/EncodeBroadHmmK562_Repressive.bed"
Transcription="$WRK/Annotation/EncodeBroadHmmK562_Transcription.bed"
Other="$WRK/Annotation/EncodeBroadHmmK562_Other.bed"

CONSERCATION=$WRK/CONSERVATION/hg38.phyloP30way.bw
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
PILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT.py
PILEUPBG=$WRK/bin/pileup_BedGraph_on_RefPT.py
DBSNP=$WRK/CONSERVATION/dbSnp153_snv.bw
#RNA
RNABAM=$WRK/BAM/*_Grocap_hg38.bam
CoPROBAMFILE=$WRK/BAM/*_CoPRO_hg38.bam

### CHANGE ME
TSS="$WRK/2025_Inr/02_TSS_NFR"

[ -d "$TSS" ] || mkdir -p "$TSS"

cd $TSS

java -jar $SCRIPTMANAGER bam-format-converter bam-to-bed -2 -o CoPRO_Read2_TSS.bed -p $CoPROBAMFILE
tail -n +2 CoPRO_Read2_TSS.bed | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' | awk '{if ($6 == "+") print $0 > "CoPRO_Read2_TSS_+.bed"; else print $0 > "CoPRO_Read2_TSS_-.bed" }' 

awk '{OFS="\t"} {print $1,$2+1,$2+1,$4,0,$6}' CoPRO_Read2_TSS_+.bed | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6}'  | uniq > CoPRO_Read2_TSS_+_1bp.bed 
awk '{OFS="\t"} {print $1,$3-1,$3-1,$4,0,$6}' CoPRO_Read2_TSS_-.bed  | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' | uniq > CoPRO_Read2_TSS_-_1bp.bed 


cat CoPRO_Read2_TSS_+_1bp.bed  CoPRO_Read2_TSS_-_1bp.bed  | bedtools sort -i | uniq > CoPRO_Read2_TSS_1bp.bed 
rm  CoPRO_Read2_TSS_+_1bp.bed  CoPRO_Read2_TSS_-_1bp.bed
rm CoPRO_Read2_TSS_+.bed CoPRO_Read2_TSS_-.bed 
bedtools intersect -v -a  CoPRO_Read2_TSS_1bp.bed  -b $BLACKLIST | awk '{if (($6 == "+") && ($1 != "ChrM") && ($2 != "-1" ) && ($3 != "1" ) && ($2 != "0" ) && ($3 != "0" )) print $0 > "CoPRO_Read2_TSS_1bp+.bed"; else if (($6 == "-") && ($1 != "ChrM") && ($2 != "-1" ) && ($3 != "1" ) && ($2 != "0" ) && ($3 != "0" )) print $0 > "CoPRO_Read2_TSS_1bp-.bed" }' 

### get CoPRO signal on each TSS
#!/bin/bash

set -e  # Exit on error

# Constants
CHUNK_SIZE=500000
INPUT_FILES=("CoPRO_Read2_TSS_1bp+.bed" "CoPRO_Read2_TSS_1bp-.bed")

mkdir -p process SCORES
mv CoPRO_Read2_TSS_1bp.bed process/

# Function to split file into chunks
split_bed_file() {
    local file="$1"
    local prefix="${file%.bed}"  # Remove .bed extension
    local counter=1
    local start=1

    total_lines=$(wc -l < "$file")

    while [ "$start" -le "$total_lines" ]; do
        end=$((start + CHUNK_SIZE - 1))
        output_file="${prefix}_${counter}.bed"
        sed -n "${start},${end}p" "$file" > "$output_file"
        ((counter++))
        start=$((end + 1))
    done
}

# Split input files
for file in "${INPUT_FILES[@]}"; do
    split_bed_file "$file"
done

# Process all split BED files
for file in CoPRO_Read2_TSS_1bp+_*.bed CoPRO_Read2_TSS_1bp-_*.bed; do
    filename="${file%.bed}"

    # Example of where analysis would go:
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2 "$file" -o "${filename}_2bp.bed"
    java -jar $SCRIPTMANAGER read-analysis tag-pileup --cpu 4 "${filename}_2bp.bed" -2 $CoPROBAMFILE -M SCORES/CoPRO_${filename}_2bp_read2

    tail -n +2 SCORES/${filename}_2bp_ENCFF663UAN_CoPRO_hg38_5read2_sense.cdt \
        | cut -f 3-4 \
        | paste "$file" - \
        | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
        > "${filename}_CoPRO.bed"

    rm "$file"
    # Optionally clean up intermediate files if uncommenting:
    rm "${filename}_2bp.bed"
    rm SCORES/CoPRO_${filename}_2bp_read2_anti.cdt
done

### take sites with TSS signal higher than 3, confirm strand
cat   CoPRO_Read2_TSS_1bp-_*_CoPRO.bed | bedtools sort -i | uniq | awk '{if (($7 >= 3 ) && ($7 >= $8 ) && ($1 !~ "U") && ($1 !~ "ChrM") && ($1 !~ "Y")) print $0 > "CoPRO_Read2_TSS-_F2_1.bed" }' 
cat   CoPRO_Read2_TSS_1bp+_*_CoPRO.bed | bedtools sort -i | uniq | awk '{if (($7 >= 3 ) && ($7 >= $8 ) && ($1 !~ "U") && ($1 !~ "ChrM") && ($1 !~ "Y")) print $0 > "CoPRO_Read2_TSS+_F2_1.bed" }' 
cat CoPRO_Read2_TSS-_F2_1.bed CoPRO_Read2_TSS+_F2_1.bed | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$7,$6}' | bedtools sort -i | uniq | bedtools intersect -v -a - -b $BLACKLIST > TSS.bed

cat CoPRO_Read2_TSS-_F2_1.bed CoPRO_Read2_TSS+_F2_1.bed | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$7,$6}' | bedtools intersect -v -a - -b $BLACKLIST | bedtools sort -i | uniq |  awk '{if (($1 !~ "U") && ($1 !~ "M") && ($1 !~ "Y") && ($1 !~ "random")) print $0 > "TSS.bed" }' 
awk '{if ( $6 == "+" ) print $0 > "TSS+.bed" ; else if ($6 == "-" ) print $0 > "TSS-.bed"  }' TSS.bed

### determine TSS flank sequence
bedtools shift -i TSS.bed -g $Genome -p -1 -m 1 > TSS_up1.bed

rm CoPRO_Read2_TSS-_F2_1.bed CoPRO_Read2_TSS+_F2_1.bed

java -jar $SCRIPTMANAGER expand-bed -c 2 TSS_up1.bed -o TSS_up1_2bp.bed
java -jar $SCRIPTMANAGER expand-bed -c 2 TSS_down1.bed -o TSS_down1_2bp.bed
java -jar $SCRIPTMANAGER expand-bed -c 2 TSS_up3.bed -o TSS_up3_2bp.bed
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME TSS_up1_2bp.bed -o TSS_up1_2bp.fa
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME TSS_down1_2bp.bed -o TSS_down1_2bp.fa
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME TSS_up3_2bp.bed -o TSS_up3_2bp.fa

rm TSS_up1.bed TSS_up1_2bp.bed TSS_down1.bed TSS_down1_2bp.bed  TSS_up3.bed TSS_up3_2bp.bed

#!/bin/bash
input_file="TSS_up1_2bp.fa"
output_file="TSS_up1_2bp.txt"
# Convert all sequences to uppercase and remove lines starting with '>'
awk '!/^>/ { print toupper($0) }' "$input_file" > "$output_file"
echo "Conversion complete. Output written to $output_file"

#!/bin/bash
input_file="TSS_down1_2bp.fa"
output_file="TSS_down1_2bp.txt"
# Convert all sequences to uppercase and remove lines starting with '>'
awk '!/^>/ { print toupper($0) }' "$input_file" > "$output_file"
echo "Conversion complete. Output written to $output_file"

#!/bin/bash
input_file="TSS_up3_2bp.fa"
output_file="TSS_up3_2bp.txt"
# Convert all sequences to uppercase and remove lines starting with '>'
awk '!/^>/ { print toupper($0) }' "$input_file" > "$output_file"
echo "Conversion complete. Output written to $output_file"

mkdir -p Inr_group
cd Inr_group
cut -f 1 ../TSS_up1_2bp.txt | paste ../TSS.bed - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"-1_"$7"_0"}' > TSS_temp.bed
cut -f 1 ../TSS_down1_2bp.txt | paste TSS_temp.bed - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_1_"$8"_2"}' > TSS_temp2.bed
cut -f 1 ../TSS_up3_2bp.txt | paste TSS_temp2.bed - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"-3_"$8"_-2_"$7}' | bedtools sort -i | uniq > TSS_test.bed

rm TSS_temp.bed TSS_temp2.bed 

java -jar $SCRIPTMANAGER peak-analysis filter-bed -e 4 -o TSS_test TSS_test.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"c"}' TSS_test_4bp-FILTER.bed > a.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"f"}' TSS_test_4bp-CLUSTER.bed > b.bed

cat a.bed b.bed | bedtools sort -i | uniq > c.bed

cut -f 7 c.bed | paste TSS_test.bed - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_"$8}' > TSS.bed

rm TSS_test.bed a.bed b.bed c.bed 
mv TSS_1bp_center.bed TSS_1bp_flanking.bed  TSS_4bp-CLUSTER.bed   TSS_4bp-FILTER.bed TSS.bed ../process/

cat TSS.bed | \
awk '{
    if ($7 ~ /_c/)  print $0 > "TSS_1bp_center.bed";
    else if ($7 ~ /_f/) print $0 > "TSS_1bp_flanking.bed";
}' 

cat TSS_1bp_center.bed | bedtools sort -i | uniq > TSS_Inr.bed
awk '{if ($6 == "+" ) print $0 > "TSS_Inr_+.bed" ; else print $0 > "TSS_Inr_-.bed" }' TSS_Inr.bed

java -jar $SCRIPTMANAGER peak-analysis filter-bed -e 40 -o TSS_Inr_+ TSS_Inr_+.bed
java -jar $SCRIPTMANAGER peak-analysis filter-bed -e 40 -o TSS_Inr_- TSS_Inr_-.bed

cat TSS_Inr_+-FILTER.bed TSS_Inr_--FILTER.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"main"}'  > a.bed
cat TSS_Inr_+-CLUSTER.bed TSS_Inr_--CLUSTER.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"alternative"}'  > b.bed

cat a.bed b.bed | bedtools sort -i | uniq > c.bed

cut -f 7 c.bed | paste TSS_Inr.bed - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_"$8}' | \
awk '{
    if ($7 ~ /_main/)  print $0 > "Main_Inr_0_1bp.bed";
    else if ($7 ~ /_alternative/) print $0 > "Alternative_Inr_0_1bp.bed";
}' 

rm a.bed b.bed c.bed TSS_Inr_+-FILTER.bed TSS_Inr_--FILTER.bed TSS_Inr_+-CLUSTER.bed TSS_Inr_--CLUSTER.bed TSS_Inr_+.bed TSS_Inr_-.bed

## Assoicate nucleosome to TSS

bedtools closest -a Main_Inr_0_1bp.bed -b "$Nuc" -d -D a -t first | bedtools sort -i | uniq | sort -k14,14n > K562_TSS_Near_Nuc_1.bed 

awk '{
  if (($14 > -90) && ($14 < 40)) 
    print $0 > "K562_TSS_Nuc_temp.bed";
  else if ($14 >= 40) 
    print $0 > "K562_TSS_NFR_1.bed";
  else 
    print $0 > "K562_TSS_NFR_temp.bed";
}' K562_TSS_Near_Nuc_1.bed 

awk '{print $0"\tNDR"}' K562_TSS_Nuc_temp.bed > K562_TSS_NDR.bed
rm K562_TSS_Nuc_temp.bed
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' K562_TSS_NFR_temp.bed | bedtools sort -i | bedtools closest -a - -b "$Nuc" -d -iu -D a -t first | cat - K562_TSS_NFR_1.bed > K562_TSS_NFR_2.bed
rm K562_TSS_NFR_temp.bed K562_TSS_NFR_1.bed
awk '{print $0"\tNFR"}' K562_TSS_NFR_2.bed > K562_TSS_NFR.bed
rm K562_TSS_NFR_2.bed 

cat K562_TSS_NFR.bed K562_TSS_NDR.bed | bedtools sort -i - | uniq | sort -k14,14n | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_"$15,$8,$9,$10,$11,$12,$13,$14 }' > K562_TSS_Near_Nuc.bed
rm K562_TSS_NFR.bed K562_TSS_NDR.bed

bedtools sort -i K562_TSS_Near_Nuc.bed | uniq | awk '{if ($7 ~/NDR/) print $0 > "K562_TSS_NDR.bed" ; else  print $0 > "K562_TSS_NFR_+1Nuc.bed"}' 

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' K562_TSS_NFR_+1Nuc.bed |  bedtools sort -i - | uniq | bedtools closest -a - -b "$Nuc" -d -io -id -D a | cut -f 8-14 | paste K562_TSS_NFR_+1Nuc.bed - | \
awk '{
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$21}' > K562_TSS_NFR_+1Nuc_-1Nuc.bed
rm K562_TSS_NFR_+1Nuc.bed 
awk '{OFS="\t"} {print $8,$9,$10,$11,$12,$6,$1,$2,$3,$4,$5,$6,$7}' K562_TSS_NDR.bed | bedtools sort -i - | uniq | bedtools closest -a - -b "$Nuc" -d -iu -io -D a -t first | awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}'  > K562_TSS_Dynamic_+1Nuc.bed 
awk '{OFS="\t"} {print $8,$9,$10,$11,$12,$6,$1,$2,$3,$4,$5,$6,$7}' K562_TSS_NDR.bed | bedtools sort -i - | uniq | bedtools closest -a - -b "$Nuc" -d -id -io -D a -t first | awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}'  > K562_TSS_Dynamic_-1Nuc.bed

 cat K562_TSS_Dynamic_-1Nuc.bed | cut -f 8-14 | paste K562_TSS_Dynamic_+1Nuc.bed - | awk '{
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$21}' > K562_TSS_NDR_+1Nuc_-1Nuc.bed

rm K562_TSS_Dynamic_-1Nuc.bed K562_TSS_Dynamic_+1Nuc.bed K562_TSS_NDR.bed

cat K562_TSS_NFR_+1Nuc_-1Nuc.bed K562_TSS_NDR_+1Nuc_-1Nuc.bed | bedtools sort -i - | uniq > K562_TSS_+1Nuc_-1Nuc.bed
rm K562_TSS_NFR_+1Nuc_-1Nuc.bed K562_TSS_NDR_+1Nuc_-1Nuc.bed

# Take NFR of each TSS

awk '{if (($6 == "+") && ($10 != "-1") &&  ($15 != "-1"))  print $0 > "TSS+_+1Nuc_-1Nuc.bed" ; else if (($6 == "-") && ($10 != "-1") &&  ($15 != "-1")) print $0 > "TSS-_+1Nuc_-1Nuc.bed" }' K562_TSS_+1Nuc_-1Nuc.bed
awk '{OFS="\t"} {print $1,$14,$9,$1"_"$14"_"$9,$9-$14,"+",$2,$3,$4,$5,$6,$7}' TSS+_+1Nuc_-1Nuc.bed | awk '{if ($5 >= 0) print $0 > "TSS+_NFR_temp.bed" }' 
bedtools intersect -v -a TSS+_NFR_temp.bed -b $tRNA | awk '{OFS="\t"} {print $1,$7,$8,$9,$10,$11,$12,$2,$3,$4,$5,$6}' | uniq | awk '{if ( $11 > 10 ) print $0 > "TSS+_NFR.bed" }' 
rm TSS+_NFR_temp.bed
awk '{OFS="\t"} {print $1,$9,$14,$1"_"$9"_"$14,$14-$9,"-",$2,$3,$4,$5,$6,$7}' TSS-_+1Nuc_-1Nuc.bed | awk '{if ($5 >= 0) print $0 > "TSS-_NFR_temp.bed" }' 
bedtools intersect -v -a TSS-_NFR_temp.bed -b $tRNA | awk '{OFS="\t"} {print $1,$7,$8,$9,$10,$11,$12,$2,$3,$4,$5,$6}' | uniq |  awk '{if ( $11 > 10 ) print $0 > "TSS-_NFR.bed" }'
rm TSS-_NFR_temp.bed

rm TSS+_+1Nuc_-1Nuc.bed TSS-_+1Nuc_-1Nuc.bed

#!/bin/bash
input_file="TSS+_NFR.bed"  
output_file="TSS+_NFR_1.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $8
    value = $5
    if (name2 in rows && value > rows[name2]) {
        rows[name2] = value
        lines[name2] = $0
    } else if (!(name2 in rows)) {
        rows[name2] = value
        lines[name2] = $0
    }
} END {
    for (name2 in lines) {
        print lines[name2]
    }
}' "$input_file" | uniq > "$output_file"

#!/bin/bash
input_file="TSS+_NFR_1.bed"  
output_file="TSS+_NFR_2.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $9
    value = $5
    if (name2 in rows && value > rows[name2]) {
        rows[name2] = value
        lines[name2] = $0
    } else if (!(name2 in rows)) {
        rows[name2] = value
        lines[name2] = $0
    }
} END {
    for (name2 in lines) {
        print lines[name2]
    }
}' "$input_file" |  uniq > "$output_file"

#!/bin/bash
input_file="TSS-_NFR.bed"  
output_file="TSS-_NFR_1.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $9
    value = $5
    if (name2 in rows && value > rows[name2]) {
        rows[name2] = value
        lines[name2] = $0
    } else if (!(name2 in rows)) {
        rows[name2] = value
        lines[name2] = $0
    }
} END {
    for (name2 in lines) {
        print lines[name2]
    }
}' "$input_file" | uniq > "$output_file"

#!/bin/bash
input_file="TSS-_NFR_1.bed" 
output_file="TSS-_NFR_2.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $8
    value = $5
    if (name2 in rows && value > rows[name2]) {
        rows[name2] = value
        lines[name2] = $0
    } else if (!(name2 in rows)) {
        rows[name2] = value
        lines[name2] = $0
    }
} END {
    for (name2 in lines) {
        print lines[name2]
    }
}' "$input_file" | uniq > "$output_file"

rm TSS-_NFR_1.bed TSS+_NFR_1.bed

cat TSS-_NFR.bed TSS+_NFR.bed > TSS_NFR.bed

cat TSS-_NFR_2.bed TSS+_NFR_2.bed | bedtools intersect -u -a - -b  TSS_NFR.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_First",$8,$9,$10,$11}' | bedtools sort -i | uniq | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$6}' | awk '{if ($11 >= 100) print $0 > "First_TSS_NFR.bed" }' 
rm TSS-_NFR_2.bed TSS+_NFR_2.bed TSS-_NFR.bed TSS+_NFR.bed
bedtools intersect -v -a TSS_NFR.bed -b First_TSS_NFR.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_Second",$8,$9,$10,$11}' | bedtools sort -i | uniq | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$6}' > Second_TSS_NFR.bed 
rm TSS_NFR.bed

## Take the distance to closest nucleosome
bedtools sort -i K562_TSS_Near_Nuc.bed | bedtools closest -a First_TSS_NFR.bed -b - -s  | awk '{if ($4 == $16) print $0 > "Temp.bed"}'
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$26,$7,$8,$9,$10,$11}' Temp.bed > First_TSS_NFR_Nuc.bed 
rm Temp.bed
bedtools sort -i K562_TSS_Near_Nuc.bed | bedtools closest -a Second_TSS_NFR.bed -b - -s  | awk '{if ($4 == $16) print $0 > "Temp.bed"}'
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$26,$7,$8,$9,$10,$11}' Temp.bed > Second_TSS_NFR_Nuc.bed
rm Temp.bed

## find TSSs share same NFR, divergent or convergent
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' First_TSS_NFR_Nuc.bed  | awk '{if ($6 == "+") print $0 > "First_TSS+.bed" ; else print $0 > "First_TSS-.bed"}'  

bedtools closest -a First_TSS+.bed -b First_TSS-.bed -d -id -io -D a -t first | \
awk '{
    if (($5 >= $15) && ($9 == $20))
        print $0 > "TSS+_TSS-_1shareNuc_NFR+.bed"; 
    else if (($5 < $15) && ($9 == $20))
        print $0 > "TSS+_TSS-_1shareNuc_NFR-.bed"; 
    else if (($5 >= $15) && ($9 >= $19) && ($9 < $20) && ($20 <= $10)) 
        print $0 > "TSS+_TSS-_shareNFR+.bed";
    else if (($5 < $15) && ($9 >= $19) && ($9 < $20) && ($20 <= $10)) 
        print $0 > "TSS+_TSS-_shareNFR-.bed";
}'
rm First_TSS+.bed First_TSS-.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Reference",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_reference_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Reference",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_reference_TSS+_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_divergent_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_divergent_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Reference",$19,$10,$1"_"$19"_"$10,$10-$19,$16,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR-.bed > First_reference_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Reference",$19,$10,$1"_"$19"_"$10,$10-$19,$16,"0shareNuc"}' TSS+_TSS-_shareNFR-.bed > First_reference_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$16,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR-.bed > First_divergent_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$16,"0shareNuc"}' TSS+_TSS-_shareNFR-.bed > First_divergent_TSS+_0shareNuc.bed

rm TSS+_TSS-_1shareNuc_NFR+.bed  TSS+_TSS-_shareNFR+.bed   TSS+_TSS-_1shareNuc_NFR-.bed TSS+_TSS-_shareNFR-.bed 
cat *_divergent_*.bed | bedtools sort -i | uniq  > First_Divergent_TSS_TU.bed
cat *_reference_*.bed | bedtools sort -i | uniq > First_Reference_TSS_TU.bed  

rm *_divergent_*.bed *_reference_*.bed 

# If TSS cannot find TSS pair, find it from second_TSS

cat First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8}' | bedtools intersect -v -a  First_TSS_NFR_Nuc.bed  -b - | bedtools sort -i | uniq | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{if ($6 == "+") print $0 > "First_TSS+.bed" ; else print $0 > "First_TSS-.bed"}' 

cat Second_TSS_NFR_Nuc.bed  | bedtools sort -i | uniq | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{if ($6 == "+") print $0 > "Second_TSS+.bed" ; else print $0 > "Second_TSS-.bed"}' 

bedtools closest -a First_TSS+.bed -b Second_TSS-.bed -d -id -io -D a | \
awk '{
    if ($9 == $20)
        print $0 > "TSS+_TSS-_1shareNuc_NFR+.bed"; 
    else if (($9 >= $19) && ($9 < $20) && ($20 <= $10)) 
        print $0 > "TSS+_TSS-_shareNFR+.bed";
}'

bedtools closest -a First_TSS-.bed -b Second_TSS+.bed -d -id -io -D a | \
awk '{
    if ($10 == $19)
        print $0 > "TSS-_TSS+_1shareNuc_NFR-.bed"; 
    else if (($10 > $19) && ($10 <= $20) && ($19 >= $9)) 
        print $0 > "TSS-_TSS+_shareNFR-.bed";
}'

rm First_TSS+.bed First_TSS-.bed Second_TSS-.bed Second_TSS+.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Secondpair",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_secondpair_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Secondpair",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_secondpair_TSS+_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > Second_divergent_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > Second_divergent_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Secondpair",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}'   TSS-_TSS+_1shareNuc_NFR-.bed > First_secondpair_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Secondpair",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}'   TSS-_TSS+_shareNFR-.bed > First_secondpair_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}' TSS-_TSS+_1shareNuc_NFR-.bed > Second_divergent_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Divergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}' TSS-_TSS+_shareNFR-.bed > Second_divergent_TSS+_0shareNuc.bed

cat Second_divergent_*.bed | bedtools sort -i | uniq  > Second_Divergent_TSS_TU.bed
cat First_secondpair_*.bed | bedtools sort -i | uniq > First_Secondpair_TSS_TU.bed  

rm Second_divergent_*.bed  First_secondpair_*.bed

rm TSS+_TSS-_1shareNuc_NFR+.bed TSS+_TSS-_shareNFR+.bed   TSS-_TSS+_1shareNuc_NFR-.bed  TSS-_TSS+_shareNFR-.bed

# If TSS cannot find divergent TSS pair first or second, find convergent first TSS pairs

cat First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed Second_Divergent_TSS_TU.bed First_Secondpair_TSS_TU.bed | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8}' | \
bedtools intersect -v -a First_TSS_NFR_Nuc.bed  -b - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' |  awk '{if ($6 == "+") print $0 > "First_TSS+.bed" ; else print $0 > "First_TSS-.bed"}' 

bedtools closest -a First_TSS+.bed -b First_TSS-.bed -d -iu -io -D a | \
awk '{
    if (($5 >= $15) && ($10 == $19))
        print $0 > "TSS+_TSS-_1shareNuc_NFR+.bed"; 
    else if (($5 < $15) && ($10 == $19))
        print $0 > "TSS+_TSS-_1shareNuc_NFR-.bed"; 
    else if (($5 >= $15) && ($10 >= $19) && ($9 <= $20) && ($19 >= $9)) 
        print $0 > "TSS+_TSS-_shareNFR+.bed";
    else if (($5 < $15) && ($10 >= $19) && ($9 <= $20) && ($19 >= $9)) 
        print $0 > "TSS+_TSS-_shareNFR-.bed";
}'

rm First_TSS+.bed  First_TSS-.bed 
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Convergentref",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_convergentpair_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Convergentref",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_convergentpair_TSS+_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_convergent_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_convergent_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergentref",$9,$20,$1"_"$9"_"$20,$20-$9,$16,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR-.bed > First_convergentpair_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergentref",$9,$20,$1"_"$9"_"$20,$20-$9,$16,"0shareNuc"}' TSS+_TSS-_shareNFR-.bed > First_convergentpair_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$16,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR-.bed > First_convergent_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$16,"0shareNuc"}' TSS+_TSS-_shareNFR-.bed > First_convergent_TSS+_0shareNuc.bed

cat First_convergent_*.bed | bedtools sort -i | uniq  > First_Convergent_TSS_TU.bed
cat First_convergentpair_*.bed | bedtools sort -i | uniq > First_Convergentpair_TSS_TU.bed  

rm First_convergent_*.bed  First_convergentpair_*.bed
rm TSS+_TSS-_1shareNuc_NFR+.bed TSS+_TSS-_shareNFR+.bed   TSS+_TSS-_1shareNuc_NFR-.bed TSS+_TSS-_shareNFR-.bed

# If TSS cannot find convergent first TSS pairs, find it from second TSS
cat First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed Second_Divergent_TSS_TU.bed First_Secondpair_TSS_TU.bed First_Convergent_TSS_TU.bed First_Convergentpair_TSS_TU.bed  | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8}' | \
bedtools intersect -v -a First_TSS_NFR_Nuc.bed  -b - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' |  awk '{if ($6 == "+") print $0 > "First_TSS+.bed" ; else print $0 > "First_TSS-.bed"}' 

bedtools intersect -v -a Second_TSS_NFR_Nuc.bed -b Second_Divergent_TSS_TU.bed | bedtools sort -i | uniq | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' |  awk '{if ($6 == "+") print $0 > "Second_TSS+.bed" ; else print $0 > "Second_TSS-.bed"}' 

bedtools closest -a First_TSS+.bed -b Second_TSS-.bed -d -iu -io -D a  | \
awk '{
    if ($10 == $19)
        print $0 > "TSS+_TSS-_1shareNuc_NFR+.bed"; 
    else if (($10 >= $19) && ($10 <= $20) && ($19 >= $9)) 
        print $0 > "TSS+_TSS-_shareNFR+.bed";
}'

bedtools closest -a First_TSS-.bed -b Second_TSS+.bed -d -iu -io -D a | \
awk '{
    if ($9 == $20)
        print $0 > "TSS-_TSS+_1shareNuc_NFR-.bed"; 
    else if (($9 >= $19) && ($9 <= $20) && ($20 <= $10)) 
        print $0 > "TSS-_TSS+_shareNFR-.bed";
}'

rm First_TSS+.bed Second_TSS-.bed First_TSS-.bed Second_TSS+.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_ConvergentrefSecondpair",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > First_convergentsecondpair_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_ConvergentrefSecondpair",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > First_convergentsecondpair_TSS+_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"1shareNuc"}' TSS+_TSS-_1shareNuc_NFR+.bed > Second_convergent_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$9,$20,$1"_"$9"_"$20,$20-$9,$6,"0shareNuc"}' TSS+_TSS-_shareNFR+.bed > Second_convergent_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_ConvergentrefSecondpair",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS-_TSS+_1shareNuc_NFR-.bed > First_convergentsecondpair_TSS-_1shareNuc.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_ConvergentrefSecondpair",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS-_TSS+_shareNFR-.bed > First_convergentsecondpair_TSS-_0shareNuc.bed

awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"1shareNuc"}' TSS-_TSS+_1shareNuc_NFR-.bed > Second_convergent_TSS+_1shareNuc.bed
awk '{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18"_Convergent",$19,$10,$1"_"$19"_"$10,$10-$19,$6,"0shareNuc"}' TSS-_TSS+_shareNFR-.bed > Second_convergent_TSS+_0shareNuc.bed

cat Second_convergent_*.bed | bedtools sort -i | uniq  > Second_Convergent_TSS_TU.bed
cat First_convergentsecondpair_*.bed | bedtools sort -i | uniq > First_ConvergentSecondpair_TSS_TU.bed  

rm Second_convergent_*.bed  First_convergentsecondpair_*.bed
rm TSS+_TSS-_1shareNuc_NFR+.bed TSS+_TSS-_shareNFR+.bed   TSS-_TSS+_1shareNuc_NFR-.bed  TSS-_TSS+_shareNFR-.bed
## Singletone TSS in TU 
cat First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed Second_Divergent_TSS_TU.bed First_Secondpair_TSS_TU.bed First_Convergent_TSS_TU.bed First_Convergentpair_TSS_TU.bed Second_Convergent_TSS_TU.bed First_ConvergentSecondpair_TSS_TU.bed | \
awk 'BEGIN {OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8}' | \
bedtools intersect -v -a First_TSS_NFR_Nuc.bed -b - | \
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Singleton",$9,$10,$1"_"$9"_"$10,$10-$9,$6,"0shareNuc"}' > First_Singleton_TSS_TU.bed

cat First_Singleton_TSS_TU.bed First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed Second_Divergent_TSS_TU.bed First_Secondpair_TSS_TU.bed | bedtools sort -i | uniq > TSS_bi_uni_sort.bed

bedtools closest -a TSS_bi_uni_sort.bed -b TSS_bi_uni_sort.bed -s -io -d -D a | \
awk '{
     if ((($10 == $23) && ($5 < $21)) || (($24 == $9) && ($5 < $21)))
        print $0 > "Flanking_TSS_to_center_TSS.bed";
}'

bedtools intersect -v -a TSS_bi_uni_sort.bed -b Flanking_TSS_to_center_TSS.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Center",$9,$10,$11,$12,$13,$14}' | \
awk '{
    if ($8 ~ /Singleton/) 
        print $0 > "First_Singleton_TSS_center_TU.bed";
    else if ($8 ~ /Reference/)
        print $0 > "First_Reference_TSS_center_TU.bed";
    else if ($8 ~ /Secondpair/)
        print $0 > "First_Secondpair_TSS_center_TU.bed";
}'

cat First_Singleton_TSS_center_TU.bed  First_Reference_TSS_center_TU.bed First_Secondpair_TSS_center_TU.bed | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13}' | \
bedtools intersect -u -a TSS_bi_uni_sort.bed -b - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Center",$9,$10,$11,$12,$13,$14}' | awk '{
    if ($8 ~ /Second_Divergent/)  
        print $0 > "Second_Divergent_TSS_center_TU.bed";
    else if ($8 ~ /Divergent/)  
        print $0 > "First_Divergent_TSS_center_TU.bed"
}'

bedtools intersect -v -a TSS_bi_uni_sort.bed -b *_center_TU.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Flanking",$9,$10,$11,$12,$13,$14}' | \
awk '{
    if ($8 ~ /Second_Divergent/)  
        print $0 > "Second_Divergent_TSS_Flanking_TU.bed";
    else if ($8 ~ /Singleton/) 
        print $0 > "First_Singleton_TSS_Flanking_TU.bed";
    else if ($8 ~ /Reference/) 
        print $0 > "First_Reference_TSS_Flanking_TU.bed";
    else if ($8 ~ /Secondpair/) 
        print $0 > "First_Secondpair_TSS_Flanking_TU.bed";
    else if ($8 ~ /Divergent/)
        print $0 > "First_Divergent_TSS_Flanking_TU.bed"
}'

rm Flanking_TSS_to_center_TSS.bed TSS_bi_uni_sort.bed

mkdir -p  First_TSS_NFR

mv *_Flanking_TU.bed    First_TSS_NFR/
mv First_Singleton_TSS_TU.bed First_Divergent_TSS_TU.bed First_Reference_TSS_TU.bed Second_Divergent_TSS_TU.bed First_Secondpair_TSS_TU.bed     First_TSS_NFR/
mv First_Convergent_TSS_TU.bed First_Convergentpair_TSS_TU.bed  Second_Convergent_TSS_TU.bed First_ConvergentSecondpair_TSS_TU.bed  First_TSS_NFR/

## Take TU with orientated TSS above threthold and TU in cononical distance

cat First_Secondpair_TSS_center_TU.bed First_Reference_TSS_center_TU.bed First_Singleton_TSS_center_TU.bed | awk '{
    if (($8 ~ /Singleton/) && ($5 >= 5) && ($12 <= 2000))
        print $0 > "Singleton_oriabove_center_TU_temp.bed";
    else if (($5 >= 5) && ($12 <= 2000))
        print $0 > "Reference_oriabove_center_TU_temp.bed";
    else
        print $0 > "Orientated_oribelow_center_TU_temp.bed"
}'

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_oriabove",$9,$10,$11,$12,$13,$14}' Reference_oriabove_center_TU_temp.bed | bedtools sort -i | uniq > Reference_oriabove_center_TU.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_oribelow",$9,$10,$11,$12,$13,$14}' Orientated_oribelow_center_TU_temp.bed | bedtools sort -i | uniq > Orientated_oribelow_center_TU.bed

awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13}' Reference_oriabove_center_TU.bed > Oriabove_center_TU.bed
awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13}' Orientated_oribelow_center_TU.bed > Oribelow_center_TU.bed

rm Reference_oriabove_center_TU_temp.bed Orientated_oribelow_center_TU_temp.bed

cat Second_Divergent_TSS_center_TU.bed First_Divergent_TSS_center_TU.bed | bedtools intersect -u -a - -b Oriabove_center_TU.bed | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_oriabove",$9,$10,$11,$12,$13,$14}' | bedtools sort -i | uniq > Divergent_oriabove_center_TU.bed

cat Second_Divergent_TSS_center_TU.bed First_Divergent_TSS_center_TU.bed | bedtools intersect -u -a - -b Oribelow_center_TU.bed | bedtools intersect -v -a - -b Divergent_oriabove_center_TU.bed | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_oribelow",$9,$10,$11,$12,$13,$14}' | bedtools sort -i | uniq > Divergent_oribelow_center_TU.bed

rm Oriabove_center_TU.bed Oribelow_center_TU.bed

mv First_Secondpair_TSS_center_TU.bed First_Reference_TSS_center_TU.bed First_Singleton_TSS_center_TU.bed Second_Divergent_TSS_center_TU.bed First_Divergent_TSS_center_TU.bed First_TSS_NFR/
mv Divergent_oribelow_center_TU.bed Orientated_oribelow_center_TU.bed First_TSS_NFR/

## calculate the anti in TU of Singleton_oriabove_center
awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13}' Singleton_oriabove_center_TU_temp.bed > TU_Singleton_Temp.bed

java -jar $SCRIPTMANAGER read-analysis tag-pileup --cpu 4 -2 TU_Singleton_Temp.bed $CoPROBAMFILE -M SCORES/TU_Singleton_Temp_CoPRO_read2
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/TU_Singleton_Temp_CoPRO_read2_anti.cdt -o SCORES/
rm SCORES/TU_Singleton_Temp_CoPRO_read2_sense.cdt
tail -n +2 SCORES/TU_Singleton_Temp_CoPRO_read2_anti_SCORES.out | cut -f 2 | paste Singleton_oriabove_center_TU_temp.bed - | awk '{
    if (($5 > $15) && ($15 < 3) )
        print $0 > "Real_Singleton_oriabove_center_TU_temp.bed";
    else 
        print $0 > "Unreal_Singleton_oriabove_center_TU_temp.bed"
}'
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_oriabove",$9,$10,$11,$12,$13,$14}' Real_Singleton_oriabove_center_TU_temp.bed > Singleton_oriabove_center_TU.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_unreal",$9,$10,$11,$12,$13,$14}' Unreal_Singleton_oriabove_center_TU_temp.bed > Singleton_unreal_center_TU.bed

mv Singleton_unreal_center_TU.bed First_TSS_NFR/
rm  Real_Singleton_oriabove_center_TU_temp.bed Unreal_Singleton_oriabove_center_TU_temp.bed   Singleton_oriabove_center_TU_temp.bed   TU_Singleton_Temp.bed

cat Singleton_oriabove_center_TU.bed Divergent_oriabove_center_TU.bed Reference_oriabove_center_TU.bed | bedtools sort -i | uniq > TSS_oriabove_center_TU.bed

#!/bin/bash
input_file="TSS_oriabove_center_TU.bed" 
output_file="TSS_oriabove_center_TU_Uniq.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $4
    value = $12
    if (name2 in rows && value < rows[name2]) {
        rows[name2] = value
        lines[name2] = $0
    } else if (!(name2 in rows)) {
        rows[name2] = value
        lines[name2] = $0
    }
} END {
    for (name2 in lines) {
        print lines[name2]
    }
}' "$input_file" | bedtools sort -i | uniq > "$output_file"  

rm TSS_oriabove_center_TU.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 80 TSS_oriabove_center_TU_Uniq.bed -o TSS_oriabove_center_TU_Uniq_80bp.bed
bedtools intersect -c -a TSS_oriabove_center_TU_Uniq_80bp.bed -b Alternative_Inr_0_1bp.bed -s > TSS_oriabove_center_TU_Uniq_80bp_Alternative.bed
cat TSS_oriabove_center_TU_Uniq_80bp_Alternative.bed | cut -f 15  | paste TSS_oriabove_center_TU_Uniq.bed -  > TSS_oriabove_center_TU_Uniq_Alter.bed
rm TSS_oriabove_center_TU_Uniq_80bp_Alternative.bed TSS_oriabove_center_TU_Uniq_80bp.bed


## plot
sort -k5,5nr TSS_oriabove_center_TU_Uniq_Alter.bed |  awk '{
    if ($8 ~ /-1_AA/) print $0 > "TSS_AA_all.bed";
    else if ($8 ~ /-1_AT/) print $0 > "TSS_AT_all.bed";
    else if ($8 ~ /-1_AG/) print $0 > "TSS_AG_all.bed";
    else if ($8 ~ /-1_AC/) print $0 > "TSS_AC_all.bed";
    else if ($8 ~ /-1_TA/) print $0 > "TSS_TA_all.bed";
    else if ($8 ~ /-1_TT/) print $0 > "TSS_TT_all.bed";
    else if ($8 ~ /-1_TG/) print $0 > "TSS_TG_all.bed";
    else if ($8 ~ /-1_TC/) print $0 > "TSS_TC_all.bed";
    else if ($8 ~ /-1_GA/) print $0 > "TSS_GA_all.bed";
    else if ($8 ~ /-1_GT/) print $0 > "TSS_GT_all.bed";
    else if ($8 ~ /-1_GG/) print $0 > "TSS_GG_all.bed";
    else if ($8 ~ /-1_GC/) print $0 > "TSS_GC_all.bed";
    else if ($8 ~ /-1_CA/) print $0 > "TSS_CA_all.bed";
    else if ($8 ~ /-1_CT/) print $0 > "TSS_CT_all.bed";
    else if ($8 ~ /-1_CG/) print $0 > "TSS_CG_all.bed";
    else if ($8 ~ /-1_CC/) print $0 > "TSS_CC_all.bed"
}' 


wc -l TSS_*_all.bed

     603 TSS_AA_all.bed
     107 TSS_AC_all.bed
     516 TSS_AG_all.bed
     176 TSS_AT_all.bed
   26268 TSS_CA_all.bed
    1469 TSS_CC_all.bed
    3699 TSS_CG_all.bed
    1129 TSS_CT_all.bed
    2132 TSS_GA_all.bed
     219 TSS_GC_all.bed
    1130 TSS_GG_all.bed
      76 TSS_GT_all.bed
    6445 TSS_TA_all.bed
     924 TSS_TC_all.bed
    4936 TSS_TG_all.bed
     350 TSS_TT_all.bed
   50179 total

cat TSS_CA_all.bed TSS_TA_all.bed TSS_GA_all.bed TSS_AA_all.bed TSS_CG_all.bed TSS_TG_all.bed TSS_GG_all.bed TSS_AG_all.bed TSS_CC_all.bed TSS_TC_all.bed TSS_GC_all.bed TSS_AC_all.bed TSS_CT_all.bed TSS_TT_all.bed TSS_GT_all.bed TSS_AT_all.bed > TSS_4color.bed


cd Inr_group/
cat TSS_CA_all.bed TSS_TA_all.bed TSS_GA_all.bed TSS_AA_all.bed TSS_CG_all.bed TSS_TG_all.bed TSS_GG_all.bed TSS_AG_all.bed TSS_CC_all.bed TSS_TC_all.bed TSS_GC_all.bed  > Inr.bed
cat TSS_AC_all.bed TSS_CT_all.bed TSS_TT_all.bed TSS_GT_all.bed TSS_AT_all.bed > nonInr.bed

cd .. 

## lable with ChromHMM annotation, tss annotation, CpG annotation

cat TSS_oriabove_center_TU_Uniq_Alter.bed | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8,$15}' | bedtools sort -i | uniq >  center_TU_TSS-PIC.bed 

bedtools intersect -u -a center_TU_TSS-PIC.bed  -b $TSSBED | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_codingTSS",$8,$9,$10,$11,$12,$13,$14,$15}' | bedtools sort -i | uniq > Center_TU_proteincoding_TSS-PIC.bed
bedtools intersect -v -a center_TU_TSS-PIC.bed  -b $TSSBED | bedtools intersect -u -a - -b $NoncodingRNABED | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_noncodingTSS",$8,$9,$10,$11,$12,$13,$14,$15}' | bedtools sort -i | uniq > Center_TU_noncoding_TSS-PIC.bed 
bedtools intersect -v -a center_TU_TSS-PIC.bed  -b $TSSBED $NoncodingRNABED  | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_other",$8,$9,$10,$11,$12,$13,$14,$15}' | bedtools sort -i | uniq > Center_TU_other_TSS-PIC.bed 

cat Center_TU_proteincoding_TSS-PIC.bed Center_TU_noncoding_TSS-PIC.bed  Center_TU_other_TSS-PIC.bed | bedtools intersect -u -a - -b $CpGBED | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_CpG",$8,$9,$10,$11,$12,$13,$14,$15}' > Center_TU_CpG_TSS-PIC.bed 
cat Center_TU_proteincoding_TSS-PIC.bed Center_TU_noncoding_TSS-PIC.bed  Center_TU_other_TSS-PIC.bed | bedtools intersect -v -a - -b Center_TU_CpG_TSS-PIC.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7"_noCpG",$8,$9,$10,$11,$12,$13,$14,$15}' > Center_TU_noCpG_TSS-PIC.bed 

rm Center_TU_proteincoding_TSS-PIC.bed Center_TU_noncoding_TSS-PIC.bed  Center_TU_other_TSS-PIC.bed


## check TU overlap ChromHMM labeling
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -u -a - -b $Enhancer | awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_Enhancer",$15}' > TSS_center_TU_Enhancer.bed
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -v -a - -b $Enhancer | bedtools intersect -u -a - -b  $Promoter    | awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_Promoter",$15}' > TSS_center_TU_Promoter.bed
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -v -a - -b $Enhancer  $Promoter  | bedtools intersect -u -a - -b  $Transcription | awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_Transcription",$15}' > TSS_center_TU_Transcription.bed
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -v -a - -b $Enhancer  $Promoter $Transcription | bedtools intersect -u -a - -b $Insulator | awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_Insulator",$15}' > TSS_center_TU_Insulator.bed
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -v -a - -b $Enhancer  $Promoter $Transcription $Insulator | bedtools intersect -u -a - -b $Repressive |  awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_Repressive",$15}'  > TSS_center_TU_Repressive.bed
cat Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed | bedtools intersect -v -a - -b $Enhancer  $Promoter $Transcription $Insulator $Repressive | awk '{OFS="\t"} {print $1,$8,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7"_NonHMM",$15}'> TSS_center_TU_NonHMM.bed

rm  Center_TU_CpG_TSS-PIC.bed  Center_TU_noCpG_TSS-PIC.bed

cat TSS_center_TU_Promoter.bed TSS_center_TU_Enhancer.bed TSS_center_TU_Transcription.bed TSS_center_TU_Insulator.bed TSS_center_TU_Repressive.bed TSS_center_TU_NonHMM.bed | \
awk '{if (($8 !~ /Divergent/) && ($14 ~ /Enhancer/)) print $0 > "Orientated_TSS_center_TU_Enhancer_PIC.bed" ; \
else if (($8 !~ /Divergent/) && ($14 ~ /Promoter/)) print $0 > "Orientated_TSS_center_TU_Promoter_PIC.bed" ; \
else if (($8 !~ /Divergent/) && ($14 ~ /Transcription/)) print $0 > "Orientated_TSS_center_TU_Transcription_PIC.bed" ; \
else if (($8 !~ /Divergent/) && ($14 ~ /Insulator/)) print $0 > "Orientated_TSS_center_TU_Insulator_PIC.bed" ; \
else if (($8 !~ /Divergent/) && ($14 ~ /Repressive/)) print $0 > "Orientated_TSS_center_TU_Repressive_PIC.bed"  ; \
else if (($8 !~ /Divergent/) && ($14 ~ /NonHMM/)) print $0 > "Orientated_TSS_center_TU_NonHMM_PIC.bed" }' 

wc -l Orientated_TSS_center_TU_Enhancer_PIC.bed
13554 Orientated_TSS_center_TU_Enhancer_PIC.bed

wc -l Orientated_TSS_center_TU_Promoter_PIC.bed
11549 Orientated_TSS_center_TU_Promoter_PIC.bed

wc -l Orientated_TSS_center_TU_NonHMM_PIC.bed
124 Orientated_TSS_center_TU_NonHMM_PIC.bed

wc -l Orientated_TSS_center_TU_Transcription_PIC.bed
1367 Orientated_TSS_center_TU_Transcription_PIC.bed

wc -l Orientated_TSS_center_TU_Insulator_PIC.bed
 66 Orientated_TSS_center_TU_Insulator_PIC.bed

wc -l Orientated_TSS_center_TU_Repressive_PIC.bed
937 Orientated_TSS_center_TU_Repressive_PIC.bed

## Calculate Bidirectional/unidirectional numbers of each  Chrom HMM labeling
cat Orientated_TSS_center_TU_Enhancer_PIC.bed Orientated_TSS_center_TU_Promoter_PIC.bed Orientated_TSS_center_TU_NonHMM_PIC.bed Orientated_TSS_center_TU_Transcription_PIC.bed Orientated_TSS_center_TU_Insulator_PIC.bed Orientated_TSS_center_TU_Repressive_PIC.bed | \
awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$2,$3,$4,$5,$6,$7,$8,$15}' | sort -k5,5n | \
awk '{if (($7 ~ /Promoter/) && ($14 ~ /_Singleton/)) print $0 > "Uni_Promoter.bed" ; \
else if ($7 ~ /Promoter/) print $0 > "Bi_Promoter.bed"; \
else if (($7 ~ /Enhancer/) && ($14 ~ /Singleton/)) print $0 > "Uni_Enhancer.bed" ; \
else if ($7 ~ /Enhancer/) print $0 > "Bi_Enhancer.bed"; \
else if (($7 ~ /Transcription/) && ($14 ~ /Singleton/)) print $0 > "Uni_Transcription.bed" ; \
else if ($7 ~ /Transcription/) print $0 > "Bi_Transcription.bed"; \
else if (($7 ~ /Insulator/) && ($14 ~ /Singleton/)) print $0 > "Uni_Insulator.bed" ; \
else if ($7 ~ /Insulator/) print $0 > "Bi_Insulator.bed"; \
else if (($7 ~ /Repressive/) && ($14 ~ /Singleton/)) print $0 > "Uni_Repressive.bed" ; \
else if ($7 ~ /Repressive/) print $0 > "Bi_Repressive.bed"; \
else if (($7 ~ /NonHMM/) && ($14 ~ /_Singleton_/)) print $0 > "Uni_NonHMM.bed" ; \
else if ($7 ~ /NonHMM/) print $0 > "Bi_NonHMM.bed" }'

wc -l Uni_Promoter.bed
wc -l Bi_Promoter.bed
wc -l Uni_Enhancer.bed
wc -l Bi_Enhancer.bed
wc -l Uni_Transcription.bed
wc -l Bi_Transcription.bed
wc -l Uni_Insulator.bed
wc -l Bi_Insulator.bed
wc -l Uni_Repressive.bed
wc -l Bi_Repressive.bed
wc -l Uni_NonHMM.bed
wc -l Bi_NonHMM.bed

    1394 Uni_Promoter.bed
   10155 Bi_Promoter.bed
    3826 Uni_Enhancer.bed
    9728 Bi_Enhancer.bed
    1332 Uni_Transcription.bed
      35 Bi_Transcription.bed
      65 Uni_Insulator.bed
       1 Bi_Insulator.bed
     800 Uni_Repressive.bed
     137 Bi_Repressive.bed
      60 Uni_NonHMM.bed
      64 Bi_NonHMM.bed

rm Uni_*.bed
rm Bi_*.bed

mkdir -p annotation_region

cat TSS_center_TU_Promoter.bed TSS_center_TU_Enhancer.bed TSS_center_TU_Transcription.bed TSS_center_TU_Insulator.bed TSS_center_TU_Repressive.bed TSS_center_TU_NonHMM.bed | bedtools sort -i | uniq > TSS_center_TU_PIC.bed

rm Orientated_TSS_center_TU_Enhancer_PIC.bed Orientated_TSS_center_TU_Promoter_PIC.bed Orientated_TSS_center_TU_NonHMM_PIC.bed Orientated_TSS_center_TU_Transcription_PIC.bed Orientated_TSS_center_TU_Insulator_PIC.bed Orientated_TSS_center_TU_Repressive_PIC.bed
mv TSS_center_TU_Promoter.bed TSS_center_TU_Enhancer.bed TSS_center_TU_Transcription.bed TSS_center_TU_Insulator.bed TSS_center_TU_Repressive.bed TSS_center_TU_NonHMM.bed annotation_region/

