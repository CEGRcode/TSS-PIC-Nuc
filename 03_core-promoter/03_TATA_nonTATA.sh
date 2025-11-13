module load gcc
#module load samtools
module load anaconda3
source activate bioinfo
# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/Path/to/Title/

###
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
MOTIFSCAN=$WRK/bin/scan_FASTA_for_motif_as_binary_string.py
PILEUPBW=$WRK//bin/pileup_BigWig_on_RefPT.py
COMPOSITE=$WRK/bin/sum_Col_CDT.pl 

GENOME=$WRK/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.txt
BLACKLIST=$WRK/hg38_files/hg38-blacklist.bed

#Conservation and SNPs
CONSERVATION=$WRK/data/Conservation-SNP/hg38.phyloP30way.bw
DBSNP=$WRK/data/Conservation-SNP/dbSnp153_snv.bw

## Determine RNA-seq and BAM file for the current job array index
BAMDIR=$WRK/data/BAM
NormDir=$WRK/data/NormalizationFactors

RNABAM=$WRK/data/BAM/*_Grocap_*.bam
CoPROBAMFILE=$WRK/data/BAM/ENCFF663UAN_CoPRO_hg38.bam
TBP=$WRK/data/BAM/*_TBP_hg38.bam
TFIIA=$WRK/data/BAM/*_GTFIIA1_hg38.bam
TBPFACTOR=`grep 'Scaling factor' $NormDir/*_TBP_*_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`
GTF2AFACTOR=`grep 'Scaling factor' $NormDir/*_GTFIIA1_*_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`

## Determine core-promoter (TATA and TATA-less) location to each TSS
corepromoter=$WRK/03_core-promoter

[ -d "$corepromoter" ] || mkdir -p "$corepromoter"

cd $corepromoter
#check if the TSS have TATA box 0 ,1,2 mismatch at upstream 
## search motif TATAWA with 2 mismatch and lable it
java -jar $SCRIPTMANAGER sequence-analysis search-motif -m TATAWAWR -n 2 -o TATAWAWR_2Mismatch_hg38.bed $GENOME
# lable it with 0 or 1 2 Mismatch
cat ../02_TSS_NFR/center_TU_TSS-PIC.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' | bedtools intersect -u -a TATAWAWR_2Mismatch_hg38.bed -b - | \
bedtools intersect -v -a - -b $BLACKLIST | bedtools sort -i | uniq > TATAWAWR.bed 
wc -l TATAWAWR.bed
250406 TATAWAWR.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 TATAWAWR.bed -o TATAWAWR_1bp.bed

bedtools sort -i $TSSBED | bedtools closest -a - -b TATAWAWR_1bp.bed -id -s -d -D a -t first  | awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13}' | \
awk '{ if (($13 > -100 ) && ($2 != "-1" )) print $0 > "TATAWAWR_codingTSS.bed" }' 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 TATAWAWR_codingTSS.bed -o TATAWAWR_codingTSS_100bp.bed
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup TATAWAWR_codingTSS_100bp.bed "$CoPROBAMFILE" -2 --cpu 4 -o CoPRO_TATAWAWR_codingTSS_100bp_read2.out
rm TATAWAWR_codingTSS_100bp.bed TATAWAWR_codingTSS.bed

for file in TATAWAWR_1bp.bed ; do
    filename=$(basename "$file" ".bed")
    # shift to down 27
    bedtools shift -i $file -g $Genome -p 27 -m -27 > ${filename}_down27.bed
    # shift to up 27
    bedtools shift -i $file -g $Genome -p -26 -m 26 > ${filename}_up26.bed
    #eampand
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 12 ${filename}_down27.bed -o ${filename}_down27_12.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 12 ${filename}_up26.bed -o ${filename}_up26_12.bed
    # Perform tag pileup analysis on the expanded BED file
    java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${filename}_down27_12.bed "$CoPROBAMFILE" -2 --cpu 4 -M CoPRO_${filename}_down27_read2
    java -jar "$SCRIPTMANAGER" read-analysis tag-pileup ${filename}_up26_12.bed "$CoPROBAMFILE" -2 --cpu 4 -M CoPRO_${filename}_up26_read2
    # Aggregate the tag pileup results
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum CoPRO_${filename}_down27_read2_sense.cdt -o SCORES/
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum CoPRO_${filename}_up26_read2_anti.cdt -o SCORES/
    # Process the scores to extract relevant information and save it
    tail -n +2 SCORES/CoPRO_${filename}_down27_read2_sense_SCORES.out | \
        cut -f 2 | \
        paste $file - > SCORES/${filename}_sense.bed
    tail -n +2 SCORES/CoPRO_${filename}_up26_read2_anti_SCORES.out | \
        cut -f 2 | \
        paste SCORES/${filename}_sense.bed - > SCORES/${filename}_sense_anti.bed
    rm ${filename}_sense.bed CoPRO_${filename}_down27_read2_anti.cdt CoPRO_${filename}_up27_read2_sense.cdt CoPRO_${filename}_down27_read2_sense.cdt CoPRO_${filename}_up27_read2_anti.cdt
    rm ${filename}_down27_12.bed ${filename}_up27_12.bed ${filename}_down27.bed  ${filename}_up27.bed
done


awk '{ if ($5 == "2") print $0 > "TATAWAWR_2mis_sense_anti.bed" ; else if ($5 == "1") print $0 > "TATAWAWR_1mis_sense_anti.bed" ; else print $0 > "TATAWAWR_0mis_sense_anti.bed" }' SCORES/TATAWAWR_1bp_sense_anti.bed

awk '{ if (($7 > $8 ) && ($7 > 2 )) print $0 > "TATAWAWR_0mis_same.bed" ; else if (($7 < $8 ) && ($8 > 2 )) print $0 > "TATAWAWR_0mis_oppo.bed" }' TATAWAWR_0mis_sense_anti.bed
awk '{ if (($7 > $8 ) && ($7 > 2 )) print $0 > "TATAWAWR_1mis_same.bed" ; else if (($7 < $8 ) && ($8 > 2 )) print $0 > "TATAWAWR_1mis_oppo.bed" }' TATAWAWR_1mis_sense_anti.bed
awk '{ if (($7 > $8 ) && ($7 > 2 )) print $0 > "TATAWAWR_2mis_same.bed" ; else if (($7 < $8 ) && ($8 > 2 )) print $0 > "TATAWAWR_2mis_oppo.bed" }' TATAWAWR_2mis_sense_anti.bed


cat ../02_TSS_NFR/TSS_center_TU_PIC.bed | awk '{ if (($14 ~ /Enhancer/ ) || ($14 ~ /Promoter/ )) print $0 > "TSS_oriabove_EP_TU_Uniq_Alter.bed" }' 

## take second TSS
cat ../02_TSS_NFR/Second_TSS_NFR_Nuc.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'| bedtools closest -a - -b TSS_oriabove_EP_TU_Uniq_Alter.bed -s -io -t first | awk '{if (($2 > $19) && ($2 < $20) && ($6 == $23)) print $0 > "temp.bed" ; else if (($2 > $19) && ($2 < $20)) print $0 > "temp2.bed" }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Refside",$9,$10,$1"_"$9"_"$10,$10-$9,$23,$24,$12}' temp.bed  > Second_TSS_NFR_Nuc_1.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_Divside",$9,$10,$1"_"$9"_"$10,$10-$9,$23,$24,$12}' temp2.bed  > Second_TSS_NFR_Nuc_2.bed
rm temp.bed temp2.bed

cat Second_TSS_NFR_Nuc_1.bed Second_TSS_NFR_Nuc_2.bed | bedtools sort -i | uniq > Second_TSS_oriabove_EP_TU_FirstTSS.bed
rm Second_TSS_NFR_Nuc_1.bed Second_TSS_NFR_Nuc_2.bed

## associat each TATA box to TSS
cat TATAWAWR_0mis_same.bed TATAWAWR_1mis_same.bed TATAWAWR_2mis_same.bed | bedtools shift -i - -g $Genome -p 27 -m -27 | bedtools sort -i | uniq > TATAWAWR_same_downstream27.bed
cat TATAWAWR_0mis_same.bed TATAWAWR_1mis_same.bed TATAWAWR_2mis_same.bed | bedtools sort -i | uniq > TATAWAWR_same.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 12 TATAWAWR_same_downstream27.bed -o TATAWAWR_same_downstream27_12bp.bed

cat TSS_oriabove_EP_TU_Uniq_Alter.bed  Second_TSS_oriabove_EP_TU_FirstTSS.bed | bedtools intersect -u -a - -b TATAWAWR_same_downstream27_12bp.bed | \
bedtools sort -i | uniq | bedtools closest -a - -b TATAWAWR_same.bed -s -id -io -d -D a | awk '{ if (($24 <= -21 ) && ($24 >= -31 )) print $0 > "TSS_all_TATA_same.bed" }' 
rm TATAWAWR_same_downstream27.bed TATAWAWR_same_downstream27_12bp.bed

cat TATAWAWR_0mis_oppo.bed TATAWAWR_1mis_oppo.bed TATAWAWR_2mis_oppo.bed | bedtools shift -i - -g $Genome -p -26 -m 26 | bedtools sort -i | uniq > TATAWAWR_same_upstream26.bed
cat TATAWAWR_0mis_oppo.bed TATAWAWR_1mis_oppo.bed TATAWAWR_2mis_oppo.bed | bedtools sort -i | uniq > TATAWAWR_oppo.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 12 TATAWAWR_same_upstream26.bed -o TATAWAWR_same_upstream26_12bp.bed
cat TSS_oriabove_EP_TU_Uniq_Alter.bed  Second_TSS_oriabove_EP_TU_FirstTSS.bed | \
bedtools intersect -v -a - -b TSS_all_TATA_same.bed | bedtools intersect -u -a - -b TATAWAWR_same_upstream26_12bp.bed | \
bedtools sort -i | uniq | bedtools closest -a - -b TATAWAWR_oppo.bed -S -id -io -d -D a | awk '{ if (($24 <= -21 ) && ($24 >= -31 )) print $0 > "TSS_all_TATA_oppo.bed" }' 

rm TATAWAWR_same_upstream26.bed TATAWAWR_same_upstream26_12bp.bed
rm TATAWAWR_0mis_oppo.bed TATAWAWR_1mis_oppo.bed TATAWAWR_2mis_oppo.bed TATAWAWR_0mis_same.bed TATAWAWR_1mis_same.bed TATAWAWR_2mis_same.bed
 
## remove duplicated TATA if associated with different TSS. Keep the one with higher CoPRO data 
#!/bin/bash
input_file="TSS_all_TATA_same.bed"  
output_file="TSS_all_TATA_same_uniq.bed" 

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $19
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
}' "$input_file" | bedtools sort -i | uniq > "$output_file"

#!/bin/bash
input_file="TSS_all_TATA_oppo.bed"  
output_file="TSS_all_TATA_oppo_uniq_1.bed"  

awk -F'\t' 'BEGIN { OFS=FS } {
    name2 = $19
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
}' "$input_file" | bedtools sort -i | uniq > "$output_file"

rm TSS_all_TATA_same.bed  TSS_all_TATA_oppo.bed 

awk '{OFS="\t"} {print $16,$17,$18,$19,$20,$21,$24,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' TSS_all_TATA_oppo_uniq_1.bed | \
bedtools intersect -v -a - -b TSS_all_TATA_same_uniq.bed | awk '{OFS="\t"} {print $8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7}' > TSS_all_TATA_oppo_uniq.bed

rm TSS_all_TATA_oppo_uniq_1.bed

### sites with TATA box in NFR but not at core-promoter
cat TSS_oriabove_EP_TU_Uniq_Alter.bed  Second_TSS_oriabove_EP_TU_FirstTSS.bed | \
bedtools intersect -v -a - -b TSS_all_TATA_same_uniq.bed TSS_all_TATA_oppo_uniq.bed | bedtools sort -i | uniq > TSS_noTATAfixed.bed

bedtools closest -a TSS_noTATAfixed.bed -b TATAWAWR_same.bed -s -d -D a -t first | \
awk '{ if ((($17 >= $9) && ($17 <= $10) && ($24 >= -20)) || (($17 >= $9) && ($17 <= $10) && ($24 <= -32))) print $0 > "TSS_all_TATA_nonfixed_same.bed"}' 

bedtools intersect -v -a TSS_noTATAfixed.bed -b TSS_all_TATA_nonfixed_same.bed | bedtools closest -a - -b TATAWAWR_oppo.bed -S -d -D a -t first  | \
awk '{ if ((($17 >= $9) && ($17 <= $10) && ($24 >= -20)) || (($17 >= $9) && ($17 <= $10) && ($24 <= -32))) print $0 > "TSS_all_TATA_nonfixed_oppo.bed" }' 

rm  TATAWAWR_oppo.bed  TATAWAWR_same.bed  TSS_noTATAfixed.bed

## lable TSS
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_TATAsame"$20"mis",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$24}' TSS_all_TATA_same_uniq.bed | \
awk '{ if ($8 ~ /oriabove/ ) print $0 > "OriaboveTSS_TATAWAWR_same.bed" ; else if ($8 ~ /Second/ ) print $0 > "SecondTSS_TATAWAWR_same.bed"  }' 
rm TSS_all_TATA_same_uniq.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_TATAoppo"$20"mis",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' TSS_all_TATA_oppo_uniq.bed | \
awk '{ if ($8 ~ /oriabove/ ) print $0 > "OriaboveTSS_TATAWAWR_oppo.bed" ; else if ($8 ~ /Second/ ) print $0 > "SecondTSS_TATAWAWR_oppo.bed"  }' 
rm TSS_all_TATA_oppo_uniq.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_NonfixedTATAsame"$20"mis",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$24}' TSS_all_TATA_nonfixed_same.bed | \
awk '{ if ($8 ~ /oriabove/ ) print $0 > "OriaboveTSS_TATAWAWR_same_nonfixed.bed" ; else if ($8 ~ /Second/ ) print $0 > "SecondTSS_TATAWAWR_same_nonfixed.bed" }' 
rm TSS_all_TATA_nonfixed_same.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_NonfixedTATAoppo"$20"mis",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$24}' TSS_all_TATA_nonfixed_oppo.bed | \
awk '{ if ($8 ~ /oriabove/ ) print $0 > "OriaboveTSS_TATAWAWR_oppo_nonfixed.bed" ; else if ($8 ~ /Second/ ) print $0 > "SecondTSS_TATAWAWR_oppo_nonfixed.bed" }' 
rm TSS_all_TATA_nonfixed_oppo.bed

### back to tss lable based on TSS type

cat OriaboveTSS_TATAWAWR_*.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14"_TATA",$15,$16,$17,$18,$19,$20,$21,$22}' | bedtools sort -i | uniq > TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed
rm OriaboveTSS_TATAWAWR_*.bed
cat SecondTSS_TATAWAWR_*.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14"_TATA",$15,$16,$17,$18,$19,$20,$21,$22}' | bedtools sort -i | uniq >   Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed
rm SecondTSS_TATAWAWR_*.bed

cat TSS_oriabove_EP_TU_Uniq_Alter.bed  | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14"_noTATA",$15,$2,$3,$4,$5,$6,$7,$8"_noTATA"}' | \
bedtools intersect -v -a - -b TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed | \
awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8}' - > TSS_oriabove_EP_TU_Uniq_Alter_noTATA.bed

cat TSS_oriabove_EP_TU_Uniq_Alter.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_noTATA",$9,$10,$11,$12,$13,$14"_TATA",$15}' | \
bedtools intersect -v -a - -b TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed TSS_oriabove_EP_TU_Uniq_Alter_noTATA.bed > TSS_oriabove_EP_TU_Uniq_Alter_TATAcontaining.bed

cat Second_TSS_oriabove_EP_TU_FirstTSS.bed  | awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14"_noTATA",$15,$2,$3,$4,$5,$6,$7,$8"_noTATA"}' | \
bedtools intersect -v -a - -b TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed  | \
awk '{OFS="\t"} {print $1,$9,$10,$11,$12,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8}' - > Second_TSS_oriabove_EP_TU_FirstTSS_noTATA.bed

cat Second_TSS_oriabove_EP_TU_FirstTSS.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8"_noTATA",$9,$10,$11,$12,$13,$14"_TATA",$15}' | \
bedtools intersect -v -a - -b TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed Second_TSS_oriabove_EP_TU_FirstTSS_noTATA.bed > Second_TSS_oriabove_EP_TU_FirstTSS_TATAcontaining.bed

mkdir -p TSS_group
mv   Second_TSS_oriabove_EP_TU_FirstTSS*.bed  TSS_oriabove_EP_TU_Uniq_Alter*.bed  TSS_group/

## change TATA strand to TSS strand, make most 5' end as 0.


cat TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed | awk '{OFS="\t"} {print $16,$17,$18,$19,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$22}' | \
awk '{ if (($14 ~ /_TATAsame/ ) || ($14 ~ /_TATAoppo/ ))  print $0 > "FixedTATA_TSS.bed"  }'  
sort -k22,22nr FixedTATA_TSS.bed | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$12,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' | \
awk '{ if ($14 ~ /TATAsame/) print $0 > "FixedTATA_TSS_same.bed" ; else print $0 > "FixedTATA_TSS_oppo.bed"}'

bedtools shift -i FixedTATA_TSS_same.bed -g $Genome -p -3 -m 3 > TATA_TSS_same_5prime.bed
bedtools shift -i FixedTATA_TSS_oppo.bed -g $Genome -p -3 -m 3 > TATA_TSS_oppo_5prime.bed

wc -l TATA_TSS_same_5prime.bed
wc -l TATA_TSS_oppo_5prime.bed
5159 TATA_TSS_same_5prime.bed
1917 TATA_TSS_oppo_5prime.bed

cat TATA_TSS_same_5prime.bed TATA_TSS_oppo_5prime.bed > TATA_TSS-strand_tsssort.bed



## take sequence around  TATA-same TSS.

awk '{ if (($14 ~ /_-1_CA_0/ ) && ($22 == "-26"))  print $0 > "TATA_TSS_CA_same_mid_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_CA_0/ ) && ($22 == "-23"))  print $0 > "TATA_TSS_CA_same_left_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_CA_0/ ) && ($22 == "-29"))  print $0 > "TATA_TSS_CA_same_right_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TA_0/ ) && ($22 == "-26"))  print $0 > "TATA_TSS_TA_same_mid_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TA_0/ ) && ($22 == "-23"))  print $0 > "TATA_TSS_TA_same_left_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TA_0/ ) && ($22 == "-29"))  print $0 > "TATA_TSS_TA_same_right_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_CG_0/ ) && ($22 == "-26"))  print $0 > "TATA_TSS_CG_same_mid_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_CG_0/ ) && ($22 == "-23"))  print $0 > "TATA_TSS_CG_same_left_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_CG_0/ ) && ($22 == "-29"))  print $0 > "TATA_TSS_CG_same_right_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TG_0/ ) && ($22 == "-26"))  print $0 > "TATA_TSS_TG_same_mid_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TG_0/ ) && ($22 == "-23"))  print $0 > "TATA_TSS_TG_same_left_5prime.bed"  }'  TATA_TSS_same_5prime.bed
awk '{ if (($14 ~ /_-1_TG_0/ ) && ($22 == "-29"))  print $0 > "TATA_TSS_TG_same_right_5prime.bed"  }'  TATA_TSS_same_5prime.bed


## make meme PWM
conda activate meme
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20 FixedTATA_TSS_same.bed -o FixedTATA_TSS_same_20bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20 FixedTATA_TSS_oppo.bed -o FixedTATA_TSS_oppo_20bp.bed

java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME FixedTATA_TSS_same_20bp.bed -o FixedTATA_TSS_same_20bp.fa
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME FixedTATA_TSS_oppo_20bp.bed -o FixedTATA_TSS_oppo_20bp.fa

mkdir -p TATA_TSS_oppo/
mkdir -p TATA_TSS_same/
meme FixedTATA_TSS_same_20bp.fa -oc TATA_TSS_oppo/ -nmotifs 1 -maxw 20 -dna
meme FixedTATA_TSS_oppo_20bp.fa -oc TATA_TSS_same/ -nmotifs 1 -maxw 20 -dna

rm FixedTATA_TSS_same_20bp.bed  FixedTATA_TSS_oppo_20bp.bed FixedTATA_TSS_same_20bp.fa FixedTATA_TSS_oppo_20bp.fa


cat TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed | \
awk '{OFS="\t"} {print $16,$17,$18,$19,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$22}' >  TSS_TATA_fixednonfixed.bed

awk '{if ($8 ~ /_TATA/) print $0 > "TSS_fixedTATA.bed"}' TSS_TATA_fixednonfixed.bed

cat TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_TATA.bed TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_TATA.bed | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' |  awk '{if ($8 ~ /NonfixedTATA/) print $0 > "TSS_nonfixedTATA_temp.bed"}' 

bedtools shift -i TSS_nonfixedTATA_temp.bed -g $Genome -p -27 -m 27 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,"0",$6}'  | cut -f 1-6 | \
paste - TSS_nonfixedTATA_temp.bed  | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,"-26"}' > TSS_nonfixedTATA.bed

rm TSS_nonfixedTATA_temp.bed

cat TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_TATAcontaining.bed TSS_group/TSS_oriabove_EP_TU_Uniq_Alter_noTATA.bed TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_noTATA.bed TSS_group/Second_TSS_oriabove_EP_TU_FirstTSS_TATAcontaining.bed | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | bedtools sort -i | uniq > TSS_noTATA_TATAcontaining_temp.bed

bedtools shift -i TSS_noTATA_TATAcontaining_temp.bed -g $Genome -p -27 -m 27 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,"0",$6}'  | cut -f 1-6 | \
paste - TSS_noTATA_TATAcontaining_temp.bed  | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,"-26"}' > TSS_noTATA_TATAcontaining.bed

rm TSS_noTATA_TATAcontaining_temp.bed 

mkdir -p TSS_type_temp

mv TSS_*.bed TSS_type_temp/


for file in TSS_type_temp/TSS_noTATA_TATAcontaining.bed TSS_type_temp/TSS_TATA_fixednonfixed.bed TSS_type_temp/TSS_nonfixedTATA.bed ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 "$file" -o ${filename}_100bp.bed
  java -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed -p -o ${filename}_100bp $GENOME  ${filename}_100bp.bed
  rm ${filename}_100bp.bed
done

mv *.cdt SCORES/

tail -n +2 SCORES/TSS_TATA_fixednonfixed_100bp_PropT.cdt | cut -f 47-53 |  paste TSS_type_temp/TSS_TATA_fixednonfixed.bed - | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$23+$24+$25+$26+$27+$28+$29,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' | \
sort -k5,5n | awk '{if (($14 ~ /_Second_Refside_NonfixedTATA/) || ($14 ~ /_Second_Divside_NonfixedTATA/)) print $0 > "TATA_ProTsort_nonfixed_SecondTSS.bed"; else if (($14 ~ /_Second_Refside_TATA/) || ($14 ~ /_Second_Divside_TATA/)) print $0 > "TATA_ProTsort_fixed_SecondTSS.bed" ; else if (($14 ~ /_Divergent_/) && ($14 ~ /_NonfixedTATA/)) print $0 > "TATA_ProTsort_nonfixed_DivergentTSS.bed"; else if (($14 ~ /_Divergent_/) && ($14 ~ /_TATA/)) print $0 > "TATA_ProTsort_fixed_DivergentTSS.bed" ; else if ($14 ~ /_TATA/) print $0 > "TATA_ProTsort_fixed_OrientatedTSS.bed" ; else print $0 > "TATA_ProTsort_nonfixed_OrientatedTSS.bed" }' 

tail -n +2 SCORES/TSS_noTATA_TATAcontaining_100bp_PropT.cdt | cut -f 47-53 |  paste TSS_type_temp/TSS_noTATA_TATAcontaining.bed - | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$23+$24+$25+$26+$27+$28+$29,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' | \
sort -k5,5n | \
awk '{if (($14 ~ /_Second_Refside/) || ($14 ~ /_Second_Divside/)) print $0 > "PIC_ProTsort_TATAcontainnoTATA_SecondTSS.bed"; else if (($14 ~ /_Divergent/) && ($20 ~ /_TATA/))  print $0 > "PIC_ProTsort_TATAcontain_DivergentTSS.bed" ; else if ($14 ~ /_Divergent/)  print $0 > "PIC_ProTsort_noTATA_DivergentTSS.bed" ; else if ($20 ~ /_TATA/) print $0 > "PIC_ProTsort_TATAcontain_OrientatedTSS.bed"; else print $0 > "PIC_ProTsort_noTATA_OrientatedTSS.bed" }' 


tail -n +2 SCORES/TSS_nonfixedTATA_100bp_PropT.cdt | cut -f 47-53 |  paste TSS_type_temp/TSS_nonfixedTATA.bed - | \
awk '{OFS="\t"} {print $1,$2,$3,$4,$23+$24+$25+$26+$27+$28+$29,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' | \
sort -k5,5n | \
awk '{if (($14 ~ /_Second_Refside_NonfixedTATA/) || ($14 ~ /_Second_Divside_NonfixedTATA/)) print $0 > "PIC_ProTsort_nonfixed_SecondTSS.bed"; else if (($14 ~ /_Divergent_/) && ($14 ~ /_NonfixedTATA/)) print $0 > "PIC_ProTsort_nonfixed_DivergentTSS.bed"; else print $0 > "PIC_ProTsort_nonfixed_OrientatedTSS.bed" }' 


mkdir -p PIC
mv PIC_ProTsort*.bed TATA_ProTsort*.bed PIC/ 
wc -l PIC/PIC*.bed
wc -l PIC/TATA*.bed

cat PIC/TATA_ProTsort_fixed_OrientatedTSS.bed | sort -k5,5n |  awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$1,$2,$3,$4,$5,$6,$22}' | \
awk '{if  ($8 ~ /TATAsame2mis/)  print $0 > "TSS_TATA2_same.bed"; else if ($8 ~ /TATAsame1mis/)  print $0 > "TSS_TATA1_same.bed"; else if ($8 ~ /TATAsame0mis/)  print $0 > "TSS_TATA0_same.bed"; else if ($8 ~ /oppo/) print $0 > "TSS_TATA_oppo.bed" }'

cat PIC/PIC_ProTsort_noTATA_OrientatedTSS.bed | sort -k5,5n | awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$1,$2,$3,$4,$5,$6,$22}' | \
awk '{if  (($8 ~ /Reference/) || ($8 ~ /Singleton/))  print $0 > "TSS_noTATA.bed" }'

cat TSS_TATA0_same.bed TSS_TATA1_same.bed TSS_TATA2_same.bed > TSS_TATA_same.bed

cat TSS_TATA_same.bed TSS_TATA_oppo.bed TSS_noTATA.bed  > TSS_TATAsame_oppo_noTATA.bed

wc -l TSS_noTATA.bed
wc -l  TSS_TATA_same.bed
wc -l  TSS_TATA_oppo.bed

   12132 TSS_noTATA.bed
    2308 TSS_TATA_same.bed
     743 TSS_TATA_oppo.bed

head -n 3033 TSS_noTATA.bed > TSS_noTATA_1.bed
tail -n +3334 TSS_noTATA.bed | head -n 3033 > TSS_noTATA_2.bed
tail -n +6067 TSS_noTATA.bed | head -n 3033 > TSS_noTATA_3.bed
tail -n +9100 TSS_noTATA.bed  >  TSS_noTATA_4.bed

cat PIC/PIC_ProTsort_TATAcontain_DivergentTSS.bed PIC/PIC_ProTsort_TATAcontain_OrientatedTSS.bed PIC/PIC_ProTsort_noTATA_DivergentTSS.bed PIC/PIC_ProTsort_noTATA_OrientatedTSS.bed PIC/PIC_ProTsort_nonfixed_DivergentTSS.bed PIC/PIC_ProTsort_nonfixed_OrientatedTSS.bed PIC/TATA_ProTsort_fixed_OrientatedTSS.bed PIC/TATA_ProTsort_fixed_DivergentTSS.bed | \
awk '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6}'  | bedtools sort -i | uniq | \
awk '{if  ($6 == "+") print $0 > "TSS_all_+.bed"; else print $0 > "TSS_all_-.bed" }'

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$10-$2}' TSS_all_+.bed > TSS_all_+_+1Nuc.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$2-$9}' TSS_all_-.bed > TSS_all_-_+1Nuc.bed

cat TSS_all_+_+1Nuc.bed TSS_all_-_+1Nuc.bed | bedtools sort -i | uniq > TSS_all.bed
rm TSS_all_+.bed  TSS_all_+_+1Nuc.bed TSS_all_-.bed  TSS_all_-_+1Nuc.bed

