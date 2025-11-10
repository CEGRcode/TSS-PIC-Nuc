module load anaconda3
source activate bioinfo
# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system
WRK=/Path/to/Title/
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
MOTIFSCAN=$WRK/bin/scan_FASTA_for_motif_as_binary_string.py
PILEUPBW=$WRK//bin/pileup_BigWig_on_RefPT.py
COMPOSITE=$WRK/bin/sum_Col_CDT.pl 
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
#Genome
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.txt
Nuc="$WRK/Annotation/BNase-Nucleosomes.bed"

## Determine RNA-seq and BAM file for the current job array index
BAMDIR=$WRK/data/BAM
NormDir=$WRK/data/NormalizationFactors

Input=$WRK/data/BAM/BNase-seq_50U-10min_merge_hg38.bam
NakedDNA=$WRK/data/BAM/NakedDNA_BNase-seq_0.04U_1_hg38.bam
CoPROBAMFILE=$WRK/data/BAM/ENCFF663UAN_CoPRO_hg38.bam


for file in ../03_core-promoter/TSS_all.bed ; do
    filename=$(basename "$file" .bed)
    awk -v filename="$filename" '{
            if ($6 == "+") {
                print $0 >> (filename "_+.bed")
            } else {
                print $0 >> (filename "_-.bed")
            }
        }' "$file"
done

awk '{OFS="\t"} {print $1,$10,$10,$1"_"$10"_"$10,$10-$2,$6,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' TSS_all_+.bed > +1Nuc_TSS_+.bed
awk '{OFS="\t"} {print $1,$9,$9,$1"_"$9"_"$9,$2-$9,$6,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' TSS_all_-.bed > +1Nuc_TSS_-.bed

cat +1Nuc_TSS_+.bed +1Nuc_TSS_-.bed | bedtools sort -i | uniq > "+1Nuc_TSS_all.bed"
rm "+1Nuc_TSS_-.bed" "+1Nuc_TSS_+.bed"

for file in +1Nuc_TSS_all.bed ; do
    filename=$(basename "$file" .bed)
    awk -v filename="+1Nuc_TSS" '{
            if ($14  ~ /Divergent/) {
                print $0 >> (filename "_Divergent.bed")
            } else {
                print $0 >> (filename "_Reference.bed")
            }
        }' "$file"
done

for file in +1Nuc_TSS_Divergent.bed +1Nuc_TSS_Reference.bed  ; do
    filename=$(basename "$file" .bed)
    awk -v filename=$filename '{
            if ($5 >= 70) {
                print $0 >> (filename "_70.bed")
            } else {
                print $0 >> (filename "_half.bed")
            }
        }' "$file"
done

mkdir -p SCORES

for file in +1Nuc_TSS_*_70.bed +1Nuc_TSS_*_half.bed ; do
    filename=`basename $file ".bed"`
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $Input --cpu 4 -5 -1 -M SCORES/BZ_${filename}_300bp_read1_original
    rm ${filename}_300bp.bed
done

for file in +1Nuc_TSS_*_half.bed; do
    filename=$(basename "$file" ".bed")
    BAM="BZ"
    tail -n +2 SCORES/${BAM}_${filename}_300bp_read1_original_sense.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t";
        print  $56+$66+$76+$86+$96,
               $57+$67+$77+$87+$97,
               $58+$68+$78+$88+$98,
               $59+$69+$79+$89+$99,
               $60+$70+$80+$90+$100,
               $61+$71+$81+$91+$101,
               $52+$62+$72+$82+$92,
               $53+$63+$73+$83+$93,
               $54+$64+$74+$84+$94,
               +$55+$65+$75+$85+$95
    }' | \
    > ${filename}_read1_original_sense_peak.bed

    tail -n +2 SCORES/${BAM}_${filename}_300bp_read1_original_anti.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t"; print \
            $56+$66+$76+$86+$96+$107+$117+$127+$137+$147,
            $57+$67+$77+$87+$97+$108+$118+$128+$138+$148,
            $58+$68+$78+$88+$98+$109+$119+$129+$139+$149,
            $59+$69+$79+$89+$99+$110+$120+$130+$140+$150,
            $50+$60+$70+$80+$90+$101+$111+$121+$131+$141,
            $51+$61+$71+$81+$91+$102+$112+$122+$132+$142,
            $52+$62+$72+$82+$92+$103+$113+$123+$133+$143,
            $53+$63+$73+$83+$93+$104+$114+$124+$134+$144,
            $54+$64+$74+$84+$94+$105+$115+$125+$135+$145,
            $55+$65+$75+$85+$95+$106+$116+$126+$136+$146
    }'  | \
    paste ${filename}_read1_original_sense_peak.bed - | \
    awk '{OFS="\t"; print $1+$19,$2+$20,$3+$11,$4+$12,$5+$13,$6+$14,$7+$15,$8+$16,$9+$17,$10+$18 }' | \
    awk '{
    OFS="\t"; 
    total = 0; 
    for (i=1; i<=10; i++) total += $i; 
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,total 
}' | \
awk '{
    OFS="\t"; 
    if ($11 == 0) $11 = 1; 
    for (i=1; i<=10; i++) printf "%s\t", $i/($11); 
    print "" 
}' | \
paste "$file" - > ${filename}_read1_original_senseanti_peak.bed
    rm ${filename}_read1_original_sense_peak.bed


    awk -v filename="${filename}_read1_original" '
    {
        max_val = $29
        max_col = 29
        max_count = 1

        for (i = 30; i <= 38; i++) {
            if ($i > max_val) {
                max_val = $i
                max_col = i
                max_count = 1
            } else if ($i == max_val) {
                max_count++
            }
        }

        if (max_count == 1) {
            offset = max_col - 28
            if (offset == 10) offset = 0
            out_file = filename "_10x_+" offset ".bed"
            print $0 > out_file
        } else {
            print $0 > (filename "_nonunique_max.bed")
        }
    }'  ${filename}_read1_original_senseanti_peak.bed

    rm "${filename}_read1_original_senseanti_peak.bed"

done

for file in +1Nuc_TSS_*_70.bed; do
    filename=$(basename "$file" ".bed")
    BAM="BZ"
    tail -n +2 SCORES/${BAM}_${filename}_300bp_read1_original_sense.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t";
        print  $5+$15+$25+$35+$45+$56+$66+$76+$86+$96,
               $6+$16+$26+$36+$46+$57+$67+$77+$87+$97,
               $7+$17+$27+$37+$47+$58+$68+$78+$88+$98,
               $8+$18+$28+$38+$48+$59+$69+$79+$89+$99,
               $9+$19+$29+$39+$49+$60+$70+$80+$90+$100,
               $10+$20+$30+$40+$50+$61+$71+$81+$91+$101,
               $1+$11+$21+$31+$41+$52+$62+$72+$82+$92,
               $2+$12+$22+$32+$42+$53+$63+$73+$83+$93,
               $3+$13+$23+$33+$43+$54+$64+$74+$84+$94,
               $4+$14+$24+$34+$44+$55+$65+$75+$85+$95
    }'  | \
    > ${filename}_read1_original_sense_peak.bed

    tail -n +2 SCORES/${BAM}_${filename}_300bp_read1_original_anti.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t"; print \
            $56+$66+$76+$86+$96+$107+$117+$127+$137+$147,
            $57+$67+$77+$87+$97+$108+$118+$128+$138+$148,
            $58+$68+$78+$88+$98+$109+$119+$129+$139+$149,
            $59+$69+$79+$89+$99+$110+$120+$130+$140+$150,
            $50+$60+$70+$80+$90+$101+$111+$121+$131+$141,
            $51+$61+$71+$81+$91+$102+$112+$122+$132+$142,
            $52+$62+$72+$82+$92+$103+$113+$123+$133+$143,
            $53+$63+$73+$83+$93+$104+$114+$124+$134+$144,
            $54+$64+$74+$84+$94+$105+$115+$125+$135+$145,
            $55+$65+$75+$85+$95+$106+$116+$126+$136+$146
    }'  | \
    paste ${filename}_read1_original_sense_peak.bed - | \
    awk '{OFS="\t"; print $1+$19,$2+$20,$3+$11,$4+$12,$5+$13,$6+$14,$7+$15,$8+$16,$9+$17,$10+$18 }' | \
    awk '{
    OFS="\t"; 
    total = 0; 
    for (i=1; i<=10; i++) total += $i; 
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,total 
}' | \
awk '{
    OFS="\t"; 
    if ($11 == 0) $11 = 1; 
    for (i=1; i<=10; i++) printf "%s\t", $i/($11); 
    print "" 
}' | \
    paste "$file" - > ${filename}_read1_original_senseanti_peak.bed
    rm ${filename}_read1_original_sense_peak.bed


    awk -v filename="${filename}_read1_original" '
    {
        max_val = $29
        max_col = 29
        max_count = 1

        for (i = 30; i <= 38; i++) {
            if ($i > max_val) {
                max_val = $i
                max_col = i
                max_count = 1
            } else if ($i == max_val) {
                max_count++
            }
        }

        if (max_count == 1) {
            offset = max_col - 28
            if (offset == 10) offset = 0
            out_file = filename "_10x_+" offset ".bed"
            print $0 > out_file
        } else {
            print $0 > (filename "_nonunique_max.bed")
        }
    }'  ${filename}_read1_original_senseanti_peak.bed

    rm "${filename}_read1_original_senseanti_peak.bed"
done

mkdir -p temp
mv  +1Nuc_TSS_*_read1_original_10x_*.bed  +1Nuc_TSS_*_read1_original_nonunique_max.bed temp/
 cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29-$34}' > Adj+1Nuc_TSS_Divergent_0.bed
 cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+2.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +1 -m -1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$30-$35}' > Adj+1Nuc_TSS_Divergent_down1.bed
  cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+3.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +2 -m -2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$31-$36}'  > Adj+1Nuc_TSS_Divergent_down2.bed
 cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+4.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +3 -m -3 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$32-$37}' > Adj+1Nuc_TSS_Divergent_down3.bed
  cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+5.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +4 -m -4 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$33-$38}' > Adj+1Nuc_TSS_Divergent_down4.bed
 cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+6.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +5 -m -5 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$34-$29}' > Adj+1Nuc_TSS_Divergent_down5.bed
cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+0.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -1 -m +1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$38-$33}' > Adj+1Nuc_TSS_Divergent_up1.bed
  cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+9.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -2 -m +2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$37-$32}' > Adj+1Nuc_TSS_Divergent_up2.bed
cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+8.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -3 -m +3 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$36-$31}' > Adj+1Nuc_TSS_Divergent_up3.bed
  cat temp/+1Nuc_TSS_Divergent_*_read1_original_10x_+7.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -4 -m +4 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$35-$30}' > Adj+1Nuc_TSS_Divergent_up4.bed 
cat temp/+1Nuc_TSS_Divergent_*_read1_original_nonunique_max.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,0}' > Adj+1Nuc_Divergent_TSS_TSS_no10x.bed


 cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29-$34}' > Adj+1Nuc_TSS_Reference_0.bed
 cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+2.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +1 -m -1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$30-$35}' > Adj+1Nuc_TSS_Reference_down1.bed
  cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+3.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +2 -m -2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$31-$36}'  > Adj+1Nuc_TSS_Reference_down2.bed
 cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+4.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +3 -m -3 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$32-$37}' > Adj+1Nuc_TSS_Reference_down3.bed
  cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+5.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +4 -m -4 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$33-$38}' > Adj+1Nuc_TSS_Reference_down4.bed
 cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+6.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +5 -m -5 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$34-$29}' > Adj+1Nuc_TSS_Reference_down5.bed
cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+0.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -1 -m +1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$38-$33}' > Adj+1Nuc_TSS_Reference_up1.bed
  cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+9.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -2 -m +2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$37-$32}' > Adj+1Nuc_TSS_Reference_up2.bed
cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+8.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -3 -m +3 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$36-$31}' > Adj+1Nuc_TSS_Reference_up3.bed
  cat temp/+1Nuc_TSS_Reference_*_read1_original_10x_+7.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -4 -m +4 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$35-$30}' > Adj+1Nuc_TSS_Reference_up4.bed 
cat temp/+1Nuc_TSS_Reference_*_read1_original_nonunique_max.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,0}' > Adj+1Nuc_Reference_TSS_TSS_no10x.bed

cat Adj+1Nuc_TSS_Divergent_0.bed Adj+1Nuc_TSS_Divergent_down*.bed Adj+1Nuc_TSS_Divergent_up*.bed Adj+1Nuc_Divergent_TSS_TSS_no10x.bed | bedtools sort | uniq > Adj+1Nuc_TSS_Divergent.bed
cat Adj+1Nuc_TSS_Reference_0.bed Adj+1Nuc_TSS_Reference_down*.bed Adj+1Nuc_TSS_Reference_up*.bed Adj+1Nuc_Reference_TSS_TSS_no10x.bed | bedtools sort | uniq > Adj+1Nuc_TSS_Reference.bed

rm Adj+1Nuc_TSS_Divergent_0.bed Adj+1Nuc_TSS_Divergent_down*.bed Adj+1Nuc_TSS_Divergent_up*.bed Adj+1Nuc_Divergent_TSS_TSS_no10x.bed
rm Adj+1Nuc_TSS_Reference_0.bed Adj+1Nuc_TSS_Reference_down*.bed Adj+1Nuc_TSS_Reference_up*.bed Adj+1Nuc_Reference_TSS_TSS_no10x.bed

##check mucleosome dinucleotide peridocity

## remove nucleosome duplication
awk '{OFS="\t"; print $1,$2,$3,$4,"0",$6}' Adj+1Nuc_TSS_Divergent.bed | bedtools sort -i | uniq > UniqNuc_TSS_Div.bed
awk '{OFS="\t"; print $1,$2,$3,$4,"0",$6}' Adj+1Nuc_TSS_Reference.bed | bedtools sort -i | uniq > UniqNuc_TSS_Ref.bed

cat UniqNuc_TSS_Div.bed UniqNuc_TSS_Ref.bed | bedtools sort -i | uniq |  awk '{
            if ($6 == "+" ) {
                print $0 > "UniqNuc_+.bed"
            } else {
                print $0 > "UniqNuc_-.bed"
            }
        }'
rm UniqNuc_TSS_Div.bed UniqNuc_TSS_Ref.bed
for file in UniqNuc_+.bed UniqNuc_-.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_300bp.bed -o ${filename}_300bp.fa
        rm ${filename}_300bp.bed
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m RR -o SCORES/RR_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m YY -o SCORES/YY_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m SS -o SCORES/SS_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m WW -o SCORES/WW_${filename}
        rm SCORES/RR_${filename}_anti.cdt
        rm SCORES/YY_${filename}_anti.cdt
        rm SCORES/SS_${filename}_anti.cdt
        rm SCORES/WW_${filename}_anti.cdt
done
mv UniqNuc_+.bed UniqNuc_-.bed SCORES/
for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_150bp.bed -o ${filename}_150bp.fa
        rm ${filename}_150bp.bed
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m RR -o SCORES/RR_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m YY -o SCORES/YY_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m SS -o SCORES/SS_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m WW -o SCORES/WW_${filename}_150
        rm SCORES/RR_${filename}_150_anti.cdt
        rm SCORES/YY_${filename}_150_anti.cdt
        rm SCORES/SS_${filename}_150_anti.cdt
        rm SCORES/WW_${filename}_150_anti.cdt
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/RR_${filename}_150_sense.cdt -o RR_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/YY_${filename}_150_sense.cdt -o YY_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/SS_${filename}_150_sense.cdt -o SS_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/WW_${filename}_150_sense.cdt -o WW_${filename}_150_sense_SCORES.out
        rm SCORES/RR_${filename}_150_sense.cdt 
        rm SCORES/YY_${filename}_150_sense.cdt  
        rm SCORES/SS_${filename}_150_sense.cdt 
        rm SCORES/WW_${filename}_150_sense.cdt
done

for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/RR_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145)-($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150);
    }' | 
    > "${filename}_RR_score_temp.txt"

    tail -n +2 RR_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_RR_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_RR_score.txt"

    rm "${filename}_RR_score_temp.txt"
done

for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/YY_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150)-($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145);
    }' | \
    > "${filename}_YY_score_temp.txt"

    tail -n +2 YY_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_YY_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_YY_score.txt"

    rm "${filename}_YY_score_temp.txt"
done

for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/SS_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143);    
    }' | \
    > "${filename}_SS_score_temp.txt"
    tail -n +2 SS_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_SS_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_SS_score.txt"

    rm "${filename}_SS_score_temp.txt"
done

for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/WW_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143)-($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148);
    }' | \
     > "${filename}_WW_score_temp.txt"
     tail -n +2 WW_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_WW_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_WW_score.txt"

    rm "${filename}_WW_score_temp.txt"
done

rm temp.out
mv *_score.txt SCORES/
mv UniqNuc_+.bed UniqNuc_-.bed SCORES/
rm *.fa *.fa.fai
mv *_SCORES.out SCORES/


for file in SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed; do
    filename=$(basename "$file" .bed)
    cat SCORES/${filename}_WW_score.txt  > ${filename}_WW.txt
    cat SCORES/${filename}_SS_score.txt  > ${filename}_SS.txt
    cat SCORES/${filename}_YY_score.txt  > ${filename}_YY.txt
    cat SCORES/${filename}_RR_score.txt  > ${filename}_RR.txt
    paste $file ${filename}_YY.txt ${filename}_RR.txt ${filename}_WW.txt ${filename}_SS.txt | awk '{OFS="\t"; print $4,$7,$8,$9,$10}' > ${filename}_matrix.cdt
    rm ${filename}_YY.txt ${filename}_RR.txt ${filename}_WW.txt ${filename}_SS.txt
done


cat UniqNuc_+_matrix.cdt UniqNuc_-_matrix.cdt > UniqNuc_matrix.cdt
rm UniqNuc_+_matrix.cdt UniqNuc_-_matrix.cdt

python maxtrixcluster.py
cat SCORES/UniqNuc_+.bed SCORES/UniqNuc_-.bed  > temp.bed
tail -n +2 sites_kmeans_clustered.tsv | cut -f 6  | paste temp.bed - | awk -v filename="UniqNuc" '{
        if ($7 == "0" ) {
            print > (filename "_0.bed")  
        } else if ($7 == "1" ) {
            print > (filename "_1.bed")  
        } else if ($7 == "2" ) {
            print > (filename "_2.bed")  
        } else if ($7 == "3" ) {
            print > (filename "_3.bed")  
        } else if ($7 == "4" ) {
            print > (filename "_4.bed")  
        } else if ($7 == "5" ) {
            print > (filename "_5.bed")  
        } else if ($7 == "6" ) {
            print > (filename "_6.bed")  
        } else if ($7 == "7" ) {
            print > (filename "_7.bed")  
        } else if ($7 == "8" ) {
            print > (filename "_8.bed")  
        } else if ($7 == "9" ) {
            print > (filename "_9.bed")  
        } else if ($7 == "10" ) {
            print > (filename "_10.bed")  
        }  else if ($7 == "11" ) {
            print > (filename "_11.bed")  
        }  else if ($7 == "12" ) {
            print > (filename "_12.bed")  
        }  else if ($7 == "13" ) {
            print > (filename "_13.bed")  
        }  else if ($7 == "14" ) {
            print > (filename "_14.bed")  
        } }'


    
for file in UniqNuc_*.bed ; do 
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_300bp.bed -o ${filename}_300bp.fa
        rm ${filename}_300bp.bed
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m SS -o SS_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m WW -o WW_${filename}
        rm WW_${filename}_anti.cdt
        rm SS_${filename}_anti.cdt
        perl "$COMPOSITE" WW_${filename}_sense.cdt WW_${filename}.out
        perl "$COMPOSITE" SS_${filename}_sense.cdt SS_${filename}.out
        rm WW_${filename}_sense.cdt
        rm SS_${filename}_sense.cdt
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m YY -o YY_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m RR -o RR_${filename}
        rm YY_${filename}_anti.cdt
        rm RR_${filename}_anti.cdt
        perl "$COMPOSITE" YY_${filename}_sense.cdt YY_${filename}.out
        perl "$COMPOSITE" RR_${filename}_sense.cdt RR_${filename}.out
        rm YY_${filename}_sense.cdt
        rm RR_${filename}_sense.cdt
        rm ${filename}_300bp.fa ${filename}_300bp.fa.fai 
done

mv sites_kmeans_clustered.tsv UniqNuc_matrix.cdt SCORES/

cat  UniqNuc_*.bed > test.bed

bedtools sort -i test.bed | uniq | bedtools closest -a Adj+1Nuc_TSS_Reference.bed -b - -d -D a -s -t first | \
awk '{OFS="\t"; print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$1,$2,$3,$4,$5,$6,$29,$36}' > TSS_Reference_Adj+1Nuc_Di.bed

bedtools sort -i test.bed | uniq | bedtools closest -a Adj+1Nuc_TSS_Divergent.bed -b - -d -D a -s -t first  | \
awk '{OFS="\t"; print $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$1,$2,$3,$4,$5,$6,$29,$36}' > TSS_Divergent_Adj+1Nuc_Di.bed
mkdir -p check
mv UniqNuc_1.bed  UniqNuc_2.bed UniqNuc_10.bed UniqNuc_4.bed  UniqNuc_5.bed UniqNuc_6.bed UniqNuc_11.bed UniqNuc_7.bed UniqNuc_3.bed UniqNuc_8.bed UniqNuc_9.bed  UniqNuc_0.bed  UniqNuc_12.bed UniqNuc_13.bed UniqNuc_14.bed check/
mv  *.out check
## make matrix 
python Correlation.py
## check each group score
for file in check/UniqNuc_*.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_300bp.bed -o ${filename}_300bp.fa
        rm ${filename}_300bp.bed
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m RR -o SCORES/RR_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m YY -o SCORES/YY_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m SS -o SCORES/SS_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m WW -o SCORES/WW_${filename}
        rm SCORES/RR_${filename}_anti.cdt
        rm SCORES/YY_${filename}_anti.cdt
        rm SCORES/SS_${filename}_anti.cdt
        rm SCORES/WW_${filename}_anti.cdt
done
for file in check/UniqNuc_*.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_150bp.bed -o ${filename}_150bp.fa
        rm ${filename}_150bp.bed
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m RR -o SCORES/RR_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m YY -o SCORES/YY_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m SS -o SCORES/SS_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m WW -o SCORES/WW_${filename}_150
        rm SCORES/RR_${filename}_150_anti.cdt
        rm SCORES/YY_${filename}_150_anti.cdt
        rm SCORES/SS_${filename}_150_anti.cdt
        rm SCORES/WW_${filename}_150_anti.cdt
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/RR_${filename}_150_sense.cdt -o RR_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/YY_${filename}_150_sense.cdt -o YY_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/SS_${filename}_150_sense.cdt -o SS_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/WW_${filename}_150_sense.cdt -o WW_${filename}_150_sense_SCORES.out
        rm SCORES/RR_${filename}_150_sense.cdt
        rm SCORES/YY_${filename}_150_sense.cdt
        rm SCORES/SS_${filename}_150_sense.cdt 
        rm SCORES/WW_${filename}_150_sense.cdt
done

for file in check/UniqNuc_*.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/RR_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145)-($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150);
    }' | 
    > "${filename}_RR_score_temp.txt"

    tail -n +2 RR_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_RR_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_RR_score.txt"

    rm "${filename}_RR_score_temp.txt"
done

for file in check/UniqNuc_*.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/YY_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150)-($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145);
    }' | \
    > "${filename}_YY_score_temp.txt"

    tail -n +2 YY_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_YY_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_YY_score.txt"

    rm "${filename}_YY_score_temp.txt"
done

for file in check/UniqNuc_*.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/SS_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143);
    }' | \
    > "${filename}_SS_score_temp.txt"
    tail -n +2 SS_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_SS_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_SS_score.txt"

    rm "${filename}_SS_score_temp.txt"
done

for file in check/UniqNuc_*.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/WW_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143)-($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148);
    }' | \
     > "${filename}_WW_score_temp.txt"
     tail -n +2 WW_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_WW_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_WW_score.txt"

    rm "${filename}_WW_score_temp.txt"
done

rm temp.out
rm *.fa *.fa.fai
rm *_SCORES.out 

for file in check/UniqNuc_*.bed ; do
    filename=$(basename "$file" .bed)

    paste ${filename}_YY_score.txt \
          ${filename}_RR_score.txt \
          ${filename}_WW_score.txt \
          ${filename}_SS_score.txt | \
    awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
        count++;
    } END {
        avg1 = sum1 / count;
        avg2 = sum2 / count;
        avg3 = sum3 / count;
        avg4 = sum4 / count;
        print avg1, avg2, avg3, avg4;
    }' > ${filename}_YY_RR_WW_SS_means.out
    rm ${filename}_YY_score.txt
    rm ${filename}_RR_score.txt
    rm ${filename}_WW_score.txt
    rm ${filename}_SS_score.txt

done
mv *_means.out check/



cat TSS_Divergent_Adj+1Nuc_Di.bed TSS_Reference_Adj+1Nuc_Di.bed | sort -k29,29nr | awk '{OFS="\t"; print $23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' | \
awk '{
            if ( $7 != "0" ) {
                print $0 > "Onlyphase_ajusted_+1Nuc_Di.bed"
            } else if ( $7 == "0" ) {
                print $0 > "nophase_ajusted_+1Nuc_Di.bed"
            } 
        }' 


awk '{
            if ($8 == "9" || $8 == "12" || $8 == "14" ) {
                print $0 > "YRWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "1" || $8 == "6" ) {
                print $0 > "sameYRlowWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "8" || $8 == "11" ) {
                print $0 > "lowYRsameWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "0" || $8 == "3" || $8 == "7" ) {
                print $0 > "sameYRantiWS_Adj+1Nuc_TSS_1.bed"
        }  else if ($8 == "2" ) {
            print > "antiYRsameWS_Adj+1Nuc_TSS_1.bed"
        } else if ($8 == "13") {
            print > "antiYRantiWS_Adj+1Nuc_TSS_1.bed"
        } else if ($8 == "5" ||$8 == "4" || $8 == "10" ) {
            print > "lessDNAencode_Adj+1Nuc_TSS_1.bed"
        }
        }' Onlyphase_ajusted_+1Nuc_Di.bed

awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"YRWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' YRWS_Adj+1Nuc_TSS_1.bed > YRWS_Adj+1Nuc_TSS.bed
rm YRWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"sameYR_lowWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' sameYRlowWS_Adj+1Nuc_TSS_1.bed > sameYRlowWS_Adj+1Nuc_TSS.bed
rm sameYRlowWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"lowYR_sameWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' lowYRsameWS_Adj+1Nuc_TSS_1.bed > lowYRsameWS_Adj+1Nuc_TSS.bed
rm lowYRsameWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"sameYR_antiWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' sameYRantiWS_Adj+1Nuc_TSS_1.bed > sameYRantiWS_Adj+1Nuc_TSS.bed
rm sameYRantiWS_Adj+1Nuc_TSS_1.bed 
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"antiYR_sameWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' antiYRsameWS_Adj+1Nuc_TSS_1.bed > antiYRsameWS_Adj+1Nuc_TSS.bed
rm antiYRsameWS_Adj+1Nuc_TSS_1.bed 
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"antiYR_antiWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' antiYRantiWS_Adj+1Nuc_TSS_1.bed > antiYRantiWS_Adj+1Nuc_TSS.bed
rm antiYRantiWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"lessDNAencode",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' lessDNAencode_Adj+1Nuc_TSS_1.bed > lessDNAencode_Adj+1Nuc_TSS.bed
rm lessDNAencode_Adj+1Nuc_TSS_1.bed

awk '{
            if ($8 == "9" || $8 == "12" || $8 == "14" ) {
                print $0 > "YRWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "1" || $8 == "6" ) {
                print $0 > "sameYRlowWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "8" || $8 == "11" ) {
                print $0 > "lowYRsameWS_Adj+1Nuc_TSS_1.bed"
            } else if ($8 == "0" || $8 == "3" || $8 == "7" ) {
                print $0 > "sameYRantiWS_Adj+1Nuc_TSS_1.bed"
        }  else if ($8 == "2" ) {
            print > "antiYRsameWS_Adj+1Nuc_TSS_1.bed"
        } else if ($8 == "13") {
            print > "antiYRantiWS_Adj+1Nuc_TSS_1.bed"
        } else if ($8 == "5" ||$8 == "4" || $8 == "10" ) {
            print > "lessDNAencode_Adj+1Nuc_TSS_1.bed"
        }
        }' nophase_ajusted_+1Nuc_Di.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"YRWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' YRWS_Adj+1Nuc_TSS_1.bed > YRWS_noAdj+1Nuc_TSS.bed
rm YRWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"sameYR_lowWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' sameYRlowWS_Adj+1Nuc_TSS_1.bed > sameYRlowWS_noAdj+1Nuc_TSS.bed
rm sameYRlowWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"lowYR_sameWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' lowYRsameWS_Adj+1Nuc_TSS_1.bed > lowYRsameWS_noAdj+1Nuc_TSS.bed
rm lowYRsameWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"sameYR_antiWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' sameYRantiWS_Adj+1Nuc_TSS_1.bed > sameYRantiWS_noAdj+1Nuc_TSS.bed
rm sameYRantiWS_Adj+1Nuc_TSS_1.bed 
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"antiYR_sameWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' antiYRsameWS_Adj+1Nuc_TSS_1.bed > antiYRsameWS_noAdj+1Nuc_TSS.bed
rm antiYRsameWS_Adj+1Nuc_TSS_1.bed 
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"antiYR_antiWS",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' antiYRantiWS_Adj+1Nuc_TSS_1.bed > antiYRantiWS_noAdj+1Nuc_TSS.bed
rm antiYRantiWS_Adj+1Nuc_TSS_1.bed
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,"lessDNAencode",$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' lessDNAencode_Adj+1Nuc_TSS_1.bed > lessDNAencode_noAdj+1Nuc_TSS.bed
rm lessDNAencode_Adj+1Nuc_TSS_1.bed


mv *_noAdj+1Nuc_TSS.bed *_Adj+1Nuc_TSS.bed  SCORES/

wc -l SCORES/*_Adj+1Nuc_TSS.bed
    9379 SCORES/YRWS_Adj+1Nuc_TSS.bed
    1745 SCORES/antiYRantiWS_Adj+1Nuc_TSS.bed
    2425 SCORES/antiYRsameWS_Adj+1Nuc_TSS.bed
   10725 SCORES/lessDNAencode_Adj+1Nuc_TSS.bed
    4336 SCORES/lowYRsameWS_Adj+1Nuc_TSS.bed
    7537 SCORES/sameYRantiWS_Adj+1Nuc_TSS.bed
    5815 SCORES/sameYRlowWS_Adj+1Nuc_TSS.bed
   41962 total
wc -l SCORES/*_noAdj+1Nuc_TSS.bed
     568 SCORES/YRWS_noAdj+1Nuc_TSS.bed
     614 SCORES/antiYRantiWS_noAdj+1Nuc_TSS.bed
     420 SCORES/antiYRsameWS_noAdj+1Nuc_TSS.bed
    1539 SCORES/lessDNAencode_noAdj+1Nuc_TSS.bed
     435 SCORES/lowYRsameWS_noAdj+1Nuc_TSS.bed
    1030 SCORES/sameYRantiWS_noAdj+1Nuc_TSS.bed
     457 SCORES/sameYRlowWS_noAdj+1Nuc_TSS.bed
    5063 total

## Conservation score

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
        filename=$(basename "$file" ".bed")
        awk '{
    key = $1"\t"$2"\t"$3"\t"$4;
    if (!(key in max) || $7+0 > max[key]) {
        max[key] = $7;
        line[key] = $0;
    }
}
END {
    for (k in line) print line[k];
}'  $file > ${filename}_dedup.bed
done 

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 ${filename}_dedup.bed -o ${filename}_300bp.bed
        python $PILEUPBW -i $CONSERVATION -r ${filename}_300bp.bed -o SCORES/phyloP30way_${filename}.cdt
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 ${filename}_dedup.bed -o ${filename}_150bp.bed
        python $PILEUPBW -i $CONSERVATION -r ${filename}_150bp.bed -o phyloP30way_${filename}_150.cdt
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum phyloP30way_${filename}_150.cdt -o SCORES/phyloP30way_${filename}_150_SCORES.out 
        rm phyloP30way_${filename}_150.cdt ${filename}_150bp.bed
        rm  ${filename}_300bp.bed
done


for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    tail -n +2 SCORES/"phyloP30way_${filename}.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147)-($8+$18+$28+$38+$48+$59+$69+$79+$89+$100+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142);
    }'  > "${filename}_phyloP30way_score_temp.txt"
    paste  ${filename}_dedup.bed "${filename}_phyloP30way_score_temp.txt" > ${filename}.bed
    mv ${filename}.bed SCORES/
    rm "${filename}_phyloP30way_score_temp.txt"
    rm ${filename}_dedup.bed
done

## YRWS pattern
for file in SCORES/*_Adj+1Nuc_TSS.bed  SCORES/*_noAdj+1Nuc_TSS.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_300bp.bed -o ${filename}_300bp.fa
        rm ${filename}_300bp.bed
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m RR -o SCORES/RR_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m YY -o SCORES/YY_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m SS -o SCORES/SS_${filename}
        python $MOTIFSCAN -i  ${filename}_300bp.fa -m WW -o SCORES/WW_${filename}
        rm SCORES/RR_${filename}_anti.cdt
        rm SCORES/YY_${filename}_anti.cdt
        rm SCORES/SS_${filename}_anti.cdt
        rm SCORES/WW_${filename}_anti.cdt
        rm ${filename}_300bp.fa ${filename}_300bp.fa.fai 
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   ${filename}_150bp.bed -o ${filename}_150bp.fa
        rm ${filename}_150bp.bed
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m RR -o SCORES/RR_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m YY -o SCORES/YY_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m SS -o SCORES/SS_${filename}_150
        python $MOTIFSCAN -i  ${filename}_150bp.fa -m WW -o SCORES/WW_${filename}_150
        rm SCORES/RR_${filename}_150_anti.cdt
        rm SCORES/YY_${filename}_150_anti.cdt
        rm SCORES/SS_${filename}_150_anti.cdt
        rm SCORES/WW_${filename}_150_anti.cdt
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/RR_${filename}_150_sense.cdt -o SCORES/RR_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/YY_${filename}_150_sense.cdt -o SCORES/YY_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/SS_${filename}_150_sense.cdt -o SCORES/SS_${filename}_150_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/WW_${filename}_150_sense.cdt -o SCORES/WW_${filename}_150_sense_SCORES.out
        rm SCORES/RR_${filename}_150_sense.cdt
        rm SCORES/YY_${filename}_150_sense.cdt
        rm SCORES/SS_${filename}_150_sense.cdt 
        rm SCORES/WW_${filename}_150_sense.cdt
done


for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/RR_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145)-($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150);
    }' | 
    > "${filename}_RR_score_temp.txt"

    tail -n +2 SCORES/RR_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_RR_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_RR_score.txt"

    rm "${filename}_RR_score_temp.txt"
done

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/YY_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148+$7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$150)-($1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133+$143+$2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145);
    }' | \
    > "${filename}_YY_score_temp.txt"

    tail -n +2 SCORES/YY_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_YY_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_YY_score.txt"

    rm "${filename}_YY_score_temp.txt"
done

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/SS_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143);
    }' | \
    > "${filename}_SS_score_temp.txt"
    tail -n +2 SCORES/SS_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_SS_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_SS_score.txt"

    rm "${filename}_SS_score_temp.txt"
done

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/WW_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$11+$21+$31+$41+$51+$62+$72+$82+$92+$102+$113+$123+$133+$143)-($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148);
    }' | \
     > "${filename}_WW_score_temp.txt"
     tail -n +2 SCORES/WW_${filename}_150_sense_SCORES.out | awk '{ OFS = "\t"; if ($2 > 0) print $2; else print 1 }' > temp.out
    paste ${filename}_WW_score_temp.txt temp.out | awk '{
        OFS = "\t";
        print $1/$2 ;
    }' > "${filename}_WW_score.txt"

    rm "${filename}_WW_score_temp.txt"
done

rm temp.out
rm *.fa *.fa.fai

for file in SCORES/*_Adj+1Nuc_TSS.bed SCORES/*_noAdj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)

    paste ${filename}_YY_score.txt \
          ${filename}_RR_score.txt \
          ${filename}_WW_score.txt \
          ${filename}_SS_score.txt | \
    awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
        count++;
    } END {
        avg1 = sum1 / count;
        avg2 = sum2 / count;
        avg3 = sum3 / count;
        avg4 = sum4 / count;
        print avg1, avg2, avg3, avg4;
    }' > ${filename}_YY_RR_WW_SS_means.out

    paste ${filename}_YY_score.txt \
          ${filename}_RR_score.txt \
          ${filename}_WW_score.txt \
          ${filename}_SS_score.txt | \
    paste $file - > ${filename}_YY_RR_WW_SS_score.bed
    rm ${filename}_YY_score.txt
    rm ${filename}_RR_score.txt
    rm ${filename}_WW_score.txt
    rm ${filename}_SS_score.txt
done


mv *_means.out SCORES/

## merge all nucleosome together
cat  sameYRlowWS_Adj+1Nuc_TSS_*score.bed lowYRsameWS_Adj+1Nuc_TSS_*score.bed YRWS_Adj+1Nuc_TSS_*score.bed   antiYRsameWS_Adj+1Nuc_TSS_*score.bed antiYRantiWS_Adj+1Nuc_TSS_*score.bed sameYRantiWS_Adj+1Nuc_TSS_*score.bed lessDNAencode_Adj+1Nuc_TSS_*score.bed  YRWS_noAdj+1Nuc_TSS_*score.bed sameYRlowWS_noAdj+1Nuc_TSS_*score.bed lowYRsameWS_noAdj+1Nuc_TSS_*score.bed sameYRantiWS_noAdj+1Nuc_TSS_*score.bed  antiYRsameWS_noAdj+1Nuc_TSS_*score.bed antiYRantiWS_noAdj+1Nuc_TSS_*score.bed lessDNAencode_noAdj+1Nuc_TSS_*score.bed > Adj+1Nuc_TSS_all.bed
rm YRWS_Adj+1Nuc_TSS_*score.bed sameYRlowWS_Adj+1Nuc_TSS_*score.bed lowYRsameWS_Adj+1Nuc_TSS_*score.bed sameYRantiWS_Adj+1Nuc_TSS_*score.bed  antiYRsameWS_Adj+1Nuc_TSS_*score.bed antiYRantiWS_Adj+1Nuc_TSS_*score.bed lessDNAencode_Adj+1Nuc_TSS_*score.bed  YRWS_noAdj+1Nuc_TSS_*score.bed sameYRlowWS_noAdj+1Nuc_TSS_*score.bed lowYRsameWS_noAdj+1Nuc_TSS_*score.bed sameYRantiWS_noAdj+1Nuc_TSS_*score.bed  antiYRsameWS_noAdj+1Nuc_TSS_*score.bed antiYRantiWS_noAdj+1Nuc_TSS_*score.bed lessDNAencode_noAdj+1Nuc_TSS_*score.bed


## add phase 

    for file in Adj+1Nuc_TSS_all.bed  ; do
        filename=Adj+1Nuc_TSS
        awk -v filename="$filename" '{
            last_digit = substr($5, length($5), 1)
            if ($5 >= 0) {
                print $0 > (filename "_phase" last_digit ".bed")
            } else {
                new_digit = 10 - last_digit
                print $0 > (filename "_phase" new_digit ".bed")
            }
        }' "$file"
    done

wc -l Adj+1Nuc_TSS_phase*.bed

    4944 Adj+1Nuc_TSS_phase0.bed
    5039 Adj+1Nuc_TSS_phase1.bed
    4926 Adj+1Nuc_TSS_phase2.bed
    4749 Adj+1Nuc_TSS_phase3.bed
    4400 Adj+1Nuc_TSS_phase4.bed
    4160 Adj+1Nuc_TSS_phase5.bed
    4196 Adj+1Nuc_TSS_phase6.bed
    4106 Adj+1Nuc_TSS_phase7.bed
    4342 Adj+1Nuc_TSS_phase8.bed
    4686 Adj+1Nuc_TSS_phase9.bed
   45548 total

awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase1",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase1.bed | sort -k27,27n > TSS_phase1_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase2",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase2.bed | sort -k27,27n > TSS_phase2_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase3",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase3.bed | sort -k27,27n > TSS_phase3_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase4",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase4.bed | sort -k27,27n > TSS_phase4_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase5",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase5.bed | sort -k27,27n > TSS_phase5_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase6",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase6.bed | sort -k27,27n > TSS_phase6_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase7",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase7.bed | sort -k27,27n > TSS_phase7_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase8",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase8.bed | sort -k27,27n > TSS_phase8_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase9",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase9.bed | sort -k27,27n > TSS_phase9_adj+1Nuc.bed
awk '{OFS="\t"; print $9,$10,$11,$12,$13,$14,$15,$16"_phase0",$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$31,$32,$33,$34,$35}' Adj+1Nuc_TSS_phase0.bed | sort -k27,27n > TSS_phase0_adj+1Nuc.bed
rm Adj+1Nuc_TSS_phase*.bed

mkdir -p temp 
mv TSS_phase*_adj+1Nuc.bed temp/

cat  temp/TSS_phase*_adj+1Nuc.bed | bedtools sort -i | uniq > TSS_all_phase_adj+1Nuc_Di.bed
awk '{OFS="\t"; print $23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' TSS_all_phase_adj+1Nuc_Di.bed | bedtools sort -i | uniq | \
> Adj+1Nuc_Di_TSS_all_phase.bed

## seperate each group by phase

cat temp/TSS_phase9_adj+1Nuc.bed temp/TSS_phase0_adj+1Nuc.bed temp/TSS_phase1_adj+1Nuc.bed temp/TSS_phase2_adj+1Nuc.bed temp/TSS_phase3_adj+1Nuc.bed temp/TSS_phase4_adj+1Nuc.bed temp/TSS_phase5_adj+1Nuc.bed temp/TSS_phase6_adj+1Nuc.bed temp/TSS_phase7_adj+1Nuc.bed temp/TSS_phase8_adj+1Nuc.bed | \
awk  '{ if ( $27 >=40 && $27 <=160 && $29 != "0" ) { print $0 > ("TSS_allphase_adj+1Nuc.bed")
     }  }' 

awk '{OFS="\t"; print $23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$1,$2,$3,$4,$5,$6,$7,$8}' TSS_allphase_adj+1Nuc.bed > adj+1Nuc_allphase_TSS.bed
awk '{OFS="\t"; print $17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7,$8}' TSS_allphase_adj+1Nuc.bed > PIC_adj+1Nuc_allphase.bed
wc -l adj+1Nuc_allphase_TSS.bed
 27838 adj+1Nuc_allphase_TSS.bed

wc -l temp/TSS_phase*_adj+1Nuc.bed
    4944 temp/TSS_phase0_adj+1Nuc.bed
    5039 temp/TSS_phase1_adj+1Nuc.bed
    4926 temp/TSS_phase2_adj+1Nuc.bed
    4749 temp/TSS_phase3_adj+1Nuc.bed
    4400 temp/TSS_phase4_adj+1Nuc.bed
    4160 temp/TSS_phase5_adj+1Nuc.bed
    4196 temp/TSS_phase6_adj+1Nuc.bed
    4106 temp/TSS_phase7_adj+1Nuc.bed
    4342 temp/TSS_phase8_adj+1Nuc.bed
    4686 temp/TSS_phase9_adj+1Nuc.bed


## choose  phase 1 and phase 6 out, plot CoPRO and DNA encoding

awk  '{ if (  $21 ~ /phase1/ || $21 ~ /phase9/ $21 ~ /phase0/ || $21 ~ /phase1/ || $21 ~ /phase2/ || $21 ~ /phase3/ ) { print $0 > ("adj+1Nuc_phase93.bed")
     }  }'  adj+1Nuc_allphase_TSS.bed

bedtools intersect -v -a adj+1Nuc_allphase_TSS.bed -b adj+1Nuc_phase93.bed > adj+1Nuc_phase48.bed 

######## +2 nucleosome and random nucleosome
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $6}' +1Nuc_TSS_all.bed | bedtools sort -i - | uniq | bedtools closest -a - -b "$Nuc" -iu -io -d -D a -t first \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $6, $13}' \
| awk 'BEGIN{OFS="\t"} $8 != "-1" {print}' > "+1Nuc_+2Nuc.bed"

awk '{OFS="\t"; print $7,$8,$9,$10,$13,$6,$1,$2,$3,$4,$5,$6 }' +1Nuc_+2Nuc.bed | awk -v filename="+2Nuc_+1Nuc" '{
            if ($5 >= 70) {
                print $0 >> (filename "_70.bed")
            } else {
                print $0 >> (filename "_half.bed")
            }
        }'

mkdir -p Nuc2_SCORES

for file in +2Nuc_+1Nuc_70.bed +2Nuc_+1Nuc_half.bed; do
    filename=`basename $file ".bed"`
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $Input --cpu 4 -5 -1 -M Nuc2_SCORES/BZ_${filename}_300bp_read1_original
    rm ${filename}_300bp.bed
done

for file in +2Nuc_+1Nuc_half.bed ; do
    filename=$(basename "$file" ".bed")
    BAM="BZ"
    tail -n +2 Nuc2_SCORES/${BAM}_${filename}_300bp_read1_original_sense.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t";
        print  $56+$66+$76+$86+$96,
               $57+$67+$77+$87+$97,
               $58+$68+$78+$88+$98,
               $59+$69+$79+$89+$99,
               $60+$70+$80+$90+$100,
               $61+$71+$81+$91+$101,
               $52+$62+$72+$82+$92,
               $53+$63+$73+$83+$93,
               $54+$64+$74+$84+$94,
               +$55+$65+$75+$85+$95
    }' | \
    > ${filename}_read1_original_sense_peak.bed

    tail -n +2 Nuc2_SCORES/${BAM}_${filename}_300bp_read1_original_anti.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t"; print \
            $56+$66+$76+$86+$96+$107+$117+$127+$137+$147,
            $57+$67+$77+$87+$97+$108+$118+$128+$138+$148,
            $58+$68+$78+$88+$98+$109+$119+$129+$139+$149,
            $59+$69+$79+$89+$99+$110+$120+$130+$140+$150,
            $50+$60+$70+$80+$90+$101+$111+$121+$131+$141,
            $51+$61+$71+$81+$91+$102+$112+$122+$132+$142,
            $52+$62+$72+$82+$92+$103+$113+$123+$133+$143,
            $53+$63+$73+$83+$93+$104+$114+$124+$134+$144,
            $54+$64+$74+$84+$94+$105+$115+$125+$135+$145,
            $55+$65+$75+$85+$95+$106+$116+$126+$136+$146
    }'  | \
    paste ${filename}_read1_original_sense_peak.bed - | \
    awk '{OFS="\t"; print $1+$19,$2+$20,$3+$11,$4+$12,$5+$13,$6+$14,$7+$15,$8+$16,$9+$17,$10+$18 }' | \
    awk '{ OFS="\t"; total = 0; for (i=1;i<=10;i++) total += $i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,total }' | \
    awk '{ OFS="\t"; for (i=1;i<=10;i++) printf "%s\t", $i/($11+1); print "" }' | \
    paste "$file" - > ${filename}_read1_original_senseanti_peak.bed
    rm ${filename}_read1_original_sense_peak.bed


    awk -v filename="${filename}_read1_original" '
    {
        max_val = $13
        max_col = 13
        max_count = 1

        for (i = 14; i <= 22; i++) {
            if ($i > max_val) {
                max_val = $i
                max_col = i
                max_count = 1
            } else if ($i == max_val) {
                max_count++
            }
        }

        if (max_count == 1) {
            offset = max_col - 12
            if (offset == 10) offset = 0
            out_file = filename "_10x_+" offset ".bed"
            print $0 > out_file
        } else {
            print $0 > (filename "_nonunique_max.bed")
        }
    }'  ${filename}_read1_original_senseanti_peak.bed

    rm "${filename}_read1_original_senseanti_peak.bed"

done

for file in +2Nuc_+1Nuc_70.bed ; do
    filename=$(basename "$file" ".bed")
    BAM="BZ"
    tail -n +2 Nuc2_SCORES/${BAM}_${filename}_300bp_read1_original_sense.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t";
        print  $5+$15+$25+$35+$45+$56+$66+$76+$86+$96,
               $6+$16+$26+$36+$46+$57+$67+$77+$87+$97,
               $7+$17+$27+$37+$47+$58+$68+$78+$88+$98,
               $8+$18+$28+$38+$48+$59+$69+$79+$89+$99,
               $9+$19+$29+$39+$49+$60+$70+$80+$90+$100,
               $10+$20+$30+$40+$50+$61+$71+$81+$91+$101,
               $1+$11+$21+$31+$41+$52+$62+$72+$82+$92,
               $2+$12+$22+$32+$42+$53+$63+$73+$83+$93,
               $3+$13+$23+$33+$43+$54+$64+$74+$84+$94,
               $4+$14+$24+$34+$44+$55+$65+$75+$85+$95
    }'  | \
    > ${filename}_read1_original_sense_peak.bed

    tail -n +2 Nuc2_SCORES/${BAM}_${filename}_300bp_read1_original_anti.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t"; print \
            $56+$66+$76+$86+$96+$107+$117+$127+$137+$147,
            $57+$67+$77+$87+$97+$108+$118+$128+$138+$148,
            $58+$68+$78+$88+$98+$109+$119+$129+$139+$149,
            $59+$69+$79+$89+$99+$110+$120+$130+$140+$150,
            $50+$60+$70+$80+$90+$101+$111+$121+$131+$141,
            $51+$61+$71+$81+$91+$102+$112+$122+$132+$142,
            $52+$62+$72+$82+$92+$103+$113+$123+$133+$143,
            $53+$63+$73+$83+$93+$104+$114+$124+$134+$144,
            $54+$64+$74+$84+$94+$105+$115+$125+$135+$145,
            $55+$65+$75+$85+$95+$106+$116+$126+$136+$146
    }'  | \
    paste ${filename}_read1_original_sense_peak.bed - | \
    awk '{OFS="\t"; print $1+$19,$2+$20,$3+$11,$4+$12,$5+$13,$6+$14,$7+$15,$8+$16,$9+$17,$10+$18 }' | \
    awk '{ OFS="\t"; total = 0; for (i=1;i<=10;i++) total += $i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,total }' | \
    awk '{ OFS="\t"; for (i=1;i<=10;i++) printf "%s\t", $i/($11+1); print "" }' | \
    paste "$file" - > ${filename}_read1_original_senseanti_peak.bed
    rm ${filename}_read1_original_sense_peak.bed


    awk -v filename="${filename}_read1_original" '
    {
        max_val = $13
        max_col = 13
        max_count = 1

        for (i = 14; i <= 22; i++) {
            if ($i > max_val) {
                max_val = $i
                max_col = i
                max_count = 1
            } else if ($i == max_val) {
                max_count++
            }
        }

        if (max_count == 1) {
            offset = max_col - 12
            if (offset == 10) offset = 0
            out_file = filename "_10x_+" offset ".bed"
            print $0 > out_file
        } else {
            print $0 > (filename "_nonunique_max.bed")
        }
    }'  ${filename}_read1_original_senseanti_peak.bed

    rm "${filename}_read1_original_senseanti_peak.bed"
done

mkdir -p Nuc2_temp
mv  +2Nuc_+1Nuc_*_read1_original_10x_*.bed  +2Nuc_+1Nuc_*_read1_original_nonunique_max.bed Nuc2_temp/
 cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13-$18}' > +2Nuc_+1Nuc_0.bed
 cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+2.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +1 -m -1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+1,$6,$7,$8,$9,$10,$11,$12,$14-$19}' > +2Nuc_+1Nuc_down1.bed
  cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+3.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +2 -m -2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+2,$6,$7,$8,$9,$10,$11,$12,$15-$20}'  > +2Nuc_+1Nuc_down2.bed
 cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+4.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +3 -m -3 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+3,$6,$7,$8,$9,$10,$11,$12,$16-$21}' > +2Nuc_+1Nuc_down3.bed
  cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+5.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +4 -m -4 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+4,$6,$7,$8,$9,$10,$11,$12,$17-$22}' > +2Nuc_+1Nuc_down4.bed
 cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+6.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +5 -m -5 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+5,$6,$7,$8,$9,$10,$11,$12,$18-$13}' > +2Nuc_+1Nuc_down5.bed
cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+0.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -1 -m +1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-1,$6,$7,$8,$9,$10,$11,$12,$22-$17}' > +2Nuc_+1Nuc_up1.bed
  cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+9.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -2 -m +2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-2,$6,$7,$8,$9,$10,$11,$12,$21-$16}' > +2Nuc_+1Nuc_up2.bed
cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+8.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -3 -m +3 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-3,$6,$7,$8,$9,$10,$11,$12,$20-$15}' > +2Nuc_+1Nuc_up3.bed
  cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_10x_+7.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -4 -m +4 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-4,$6,$7,$8,$9,$10,$11,$12,$19-$14}' > +2Nuc_+1Nuc_up4.bed 
cat Nuc2_temp/+2Nuc_+1Nuc_*_read1_original_nonunique_max.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0"}' > +2Nuc_+1Nuc_no10x.bed

cat +2Nuc_+1Nuc_0.bed +2Nuc_+1Nuc_down*.bed +2Nuc_+1Nuc_up*.bed +2Nuc_+1Nuc_no10x.bed | bedtools sort | uniq > Adj+2Nuc_+1Nuc.bed

rm +2Nuc_+1Nuc_0.bed +2Nuc_+1Nuc_down*.bed +2Nuc_+1Nuc_up*.bed +2Nuc_+1Nuc_no10x.bed

########## 
# Other nucleosome
##########


cat TSS_all_+.bed TSS_all_-.bed | awk 'BEGIN{OFS="\t"} {print $1, $9, $10, $11, $12, $13}' | bedtools intersect -v -a $Nuc -b - | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' > random_Nuc.bed
awk 'BEGIN{srand()} {print rand()"\t"$0}' random_Nuc.bed | sort -k1,1n | cut -f2- | head -n 40000 > random_Nuc_40000.bed

java -jar $SCRIPTMANAGER peak-analysis filter-bed -e 200 -o random_Nuc random_Nuc_40000.bed
rm random_Nuc-CLUSTER.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1"_"$2"_"$3, $5, $6}' random_Nuc-FILTER.bed | awk '{
            if ( $1 !~ /K/ && $1 !~ /GL/ ) {
                print $0 > "rNuc.bed"
            } 
        }'  
rm random_Nuc-FILTER.bed


mkdir -p random_SCORES

for file in rNuc.bed ; do
    filename=`basename $file ".bed"`
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $Input --cpu 4 -5 -1 -M random_SCORES/BZ_${filename}_300bp_read1_original
    rm ${filename}_300bp.bed
done


for file in rNuc.bed ; do
    filename=$(basename "$file" ".bed")
    BAM="BZ"
    tail -n +2 random_SCORES/${BAM}_${filename}_300bp_read1_original_sense.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t";
        print  $5+$15+$25+$35+$45+$56+$66+$76+$86+$96,
               $6+$16+$26+$36+$46+$57+$67+$77+$87+$97,
               $7+$17+$27+$37+$47+$58+$68+$78+$88+$98,
               $8+$18+$28+$38+$48+$59+$69+$79+$89+$99,
               $9+$19+$29+$39+$49+$60+$70+$80+$90+$100,
               $10+$20+$30+$40+$50+$61+$71+$81+$91+$101,
               $1+$11+$21+$31+$41+$52+$62+$72+$82+$92,
               $2+$12+$22+$32+$42+$53+$63+$73+$83+$93,
               $3+$13+$23+$33+$43+$54+$64+$74+$84+$94,
               $4+$14+$24+$34+$44+$55+$65+$75+$85+$95
    }'  | \
    > ${filename}_read1_original_sense_peak.bed

    tail -n +2 random_SCORES/${BAM}_${filename}_300bp_read1_original_anti.cdt | cut -f 78-227 | \
    awk '{
        OFS="\t"; print \
            $56+$66+$76+$86+$96+$107+$117+$127+$137+$147,
            $57+$67+$77+$87+$97+$108+$118+$128+$138+$148,
            $58+$68+$78+$88+$98+$109+$119+$129+$139+$149,
            $59+$69+$79+$89+$99+$110+$120+$130+$140+$150,
            $50+$60+$70+$80+$90+$101+$111+$121+$131+$141,
            $51+$61+$71+$81+$91+$102+$112+$122+$132+$142,
            $52+$62+$72+$82+$92+$103+$113+$123+$133+$143,
            $53+$63+$73+$83+$93+$104+$114+$124+$134+$144,
            $54+$64+$74+$84+$94+$105+$115+$125+$135+$145,
            $55+$65+$75+$85+$95+$106+$116+$126+$136+$146
    }'  | \
    paste ${filename}_read1_original_sense_peak.bed - | \
    awk '{OFS="\t"; print $1+$19,$2+$20,$3+$11,$4+$12,$5+$13,$6+$14,$7+$15,$8+$16,$9+$17,$10+$18 }' | \
    awk '{ OFS="\t"; total = 0; for (i=1;i<=10;i++) total += $i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,total }' | \
    awk '{ OFS="\t"; for (i=1;i<=10;i++) printf "%s\t", $i/($11+1); print "" }' | \
    paste "$file" - > ${filename}_read1_original_senseanti_peak.bed
    rm ${filename}_read1_original_sense_peak.bed


    awk -v filename="${filename}_read1_original" '
    {
        max_val = $7
        max_col = 7
        max_count = 1

        for (i = 8; i <= 16; i++) {
            if ($i > max_val) {
                max_val = $i
                max_col = i
                max_count = 1
            } else if ($i == max_val) {
                max_count++
            }
        }

        if (max_count == 1) {
            offset = max_col - 6
            if (offset == 10) offset = 0
            out_file = filename "_10x_+" offset ".bed"
            print $0 > out_file
        } else {
            print $0 > (filename "_nonunique_max.bed")
        }
    }'  ${filename}_read1_original_senseanti_peak.bed

    rm "${filename}_read1_original_senseanti_peak.bed"
done
mkdir -p random_temp
mv  rNuc_read1_original_10x_*.bed  rNuc_read1_original_nonunique_max.bed random_temp/
 cat random_temp/rNuc_read1_original_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7-$12}' > rNuc_0.bed
 cat random_temp/rNuc_read1_original_10x_+2.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +1 -m -1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+1,$6,$8-$13}' > rNuc_down1.bed
  cat random_temp/rNuc_read1_original_10x_+3.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +2 -m -2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+2,$6,$9-$14}'  > rNuc_down2.bed
 cat random_temp/rNuc_read1_original_10x_+4.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +3 -m -3 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+3,$6,$10-$15}' > rNuc_down3.bed
  cat random_temp/rNuc_read1_original_10x_+5.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +4 -m -4 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+4,$6,$11-$16}' > rNuc_down4.bed
 cat random_temp/rNuc_read1_original_10x_+6.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p +5 -m -5 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5+5,$6,$12-$7}' > rNuc_down5.bed
cat random_temp/rNuc_read1_original_10x_+0.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -1 -m +1 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-1,$6,$16-$11}' > rNuc_up1.bed
  cat random_temp/rNuc_read1_original_10x_+9.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -2 -m +2 | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-2,$6,$15-$10}' > rNuc_up2.bed
cat random_temp/rNuc_read1_original_10x_+8.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -3 -m +3 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-3,$6,$14-$9}' > rNuc_up3.bed
  cat random_temp/rNuc_read1_original_10x_+7.bed | bedtools sort -i | uniq | bedtools shift -i - -g $Genome -p -4 -m +4 |  awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3,$5-4,$6,$13-$8}' > rNuc_up4.bed 
cat random_temp/rNuc_read1_original_nonunique_max.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,"0"}' > rNuc_no10x.bed
cat rNuc_0.bed rNuc_down*.bed rNuc_up*.bed rNuc_no10x.bed | awk '{
            if ( $2 >= 0 && $3 >= 0 ) {
                print $0 > "test.bed"
            } 
        }'  

bedtools sort -i test.bed | uniq > AdjrNuc.bed

rm  rNuc_0.bed rNuc_down*.bed rNuc_up*.bed rNuc_no10x.bed




