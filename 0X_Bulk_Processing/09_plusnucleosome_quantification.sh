module load anaconda3
source activate bioinfo
# Script for F2a and E5

### CHANGE ME
WRK=/Path/to/Title/
Reference=$WRK/X_Bulk_Processing/Reference
BAMDIR=$WRK/data/BAM
###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
[ -d $WRK/Library/E6 ] || mkdir -p $WRK/Library/E6
[ -d $WRK/Library/F3d ] || mkdir -p $WRK/Library/F3d
[ -d $WRK/Library/E8 ] || mkdir -p $WRK/Library/E8
[ -d $WRK/Library/F3e ] || mkdir -p $WRK/Library/F3e
[ -d $WRK/Library/F3g ] || mkdir -p $WRK/Library/F3g

BAMDIR=$WRK/data/BAM
cd $WRK/Library/E6


for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed $WRK/04_plusoneNucleosome/SCORES/*_noAdj+1Nuc_TSS.bed ; do
        filename=$(basename "$file" ".bed")
        cp $WRK/04_plusoneNucleosome/SCORES/${filename}_YY_RR_WW_SS_means.out $WRK/Library/E6
done

cd $WRK/Library/F3d

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed $WRK/04_plusoneNucleosome/SCORES/*_noAdj+1Nuc_TSS.bed  ; do
    filename=$(basename "$file" .bed)
    awk '{
        sum1 += $31;
        count++;
    } END {
        avg1 = sum1 / count;
        print avg1;
    }' > $WRK/Library/F3d/${filename}_Conservation_means.out
done

cd $WRK/Library/E8

## times conservation score to YRWS pattern

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed   ; do
        filename=$(basename "$file" ".bed")

## reorder
awk -F'\t' '
NR==1 {
    print
    next
}
{
    # Copy original values
    for (i = 1; i <= NF; i++) {
        orig[i] = $i
    }

    # Only check from column 4 to end (3rd index = data starts)
    for (i = 4; i <= NF; i++) {
        if (orig[i] == 1) {
            $ (i - 1) = 1
        }
    }

    print
}' OFS='\t' $WRK/Library/F3c/CDT/SS_${filename}_sense.cdt > $WRK/Library/E8/SS_${filename}_sense_reorder.cdt


head -n 1  $WRK/Library/F3c/CDT/SS_${filename}_sense.cdt  > $WRK/Library/E8/SSxphyloP30way_"${filename}"_sense.cdt
paste $WRK/Library/E8/SS_${filename}_sense_reorder.cdt $WRK/Library/E8/CDT/hg38.phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
  yorfn = $1
  name = $2
  total_fields = NF
  score_count = (total_fields / 2) - 2   # number of score columns
  line = yorfn OFS name

  # Loop through score1 .. scoreN
  for (i = 1; i <= score_count; i++) {
    c1 = 2 + i
    c2 = (total_fields / 2) + i + 1
    product = $(c1) * $(c2)
    if (product == 0) {
      sum = "-"
    } else {
      sum = product
    }
    line = line OFS sum
  }

  print line
}' >> $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt


head -n 1 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt  > $WRK/Library/E8/CDT/sensehead.cdt

tail -n +3 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt |  cat $WRK/Library/E8/CDT/sensehead.cdt - > $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt

rm $WRK/Library/E8/CDT/sensehead.cdt $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt $WRK/Library/E8/SS_${filename}_sense_reorder.cdt

perl $COMPOSITEFILTER $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/SSxphyloP30way_"${filename}"_sense_reorder.out

done

## times conservation score to YR pattern

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed   ; do
        filename=$(basename "$file" ".bed")

## reorder
awk -F'\t' '
NR==1 {
    print
    next
}
{
    # Copy original values
    for (i = 1; i <= NF; i++) {
        orig[i] = $i
    }

    # Only check from column 4 to end (3rd index = data starts)
    for (i = 4; i <= NF; i++) {
        if (orig[i] == 1) {
            $ (i - 1) = 1
        }
    }

    print
}' OFS='\t' $WRK/Library/F3c/CDT/WW_${filename}_sense.cdt > $WRK/Library/E8/WW_${filename}_sense_reorder.cdt


head -n 1  $WRK/Library/F3c/CDT/WW_${filename}_sense.cdt  > $WRK/Library/E8/SSxphyloP30way_"${filename}"_sense.cdt
paste $WRK/Library/E8/WW_${filename}_sense_reorder.cdt $WRK/Library/E8/CDT/hg38.phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
  yorfn = $1
  name = $2
  total_fields = NF
  score_count = (total_fields / 2) - 2   # number of score columns
  line = yorfn OFS name

  # Loop through score1 .. scoreN
  for (i = 1; i <= score_count; i++) {
    c1 = 2 + i
    c2 = (total_fields / 2) + i + 1
    product = $(c1) * $(c2)
    if (product == 0) {
      sum = "-"
    } else {
      sum = product
    }
    line = line OFS sum
  }

  print line
}' >> $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt


head -n 1 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt  > $WRK/Library/E8/CDT/sensehead.cdt

tail -n +3 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt |  cat $WRK/Library/E8/CDT/sensehead.cdt - > $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt

rm $WRK/Library/E8/CDT/sensehead.cdt $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt $WRK/Library/E8/WW_${filename}_sense_reorder.cdt

perl $COMPOSITEFILTER $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/SSxphyloP30way_"${filename}"_sense_reorder.out

done

## times conservation score to YRWS pattern

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed   ; do
        filename=$(basename "$file" ".bed")

## reorder
awk -F'\t' '
NR==1 {
    print
    next
}
{
    # Copy original values
    for (i = 1; i <= NF; i++) {
        orig[i] = $i
    }

    # Only check from column 4 to end (3rd index = data starts)
    for (i = 4; i <= NF; i++) {
        if (orig[i] == 1) {
            $ (i - 1) = 1
        }
    }

    print
}' OFS='\t' $WRK/Library/F3c/CDT/YY_${filename}_sense.cdt > $WRK/Library/E8/YY_${filename}_sense_reorder.cdt


head -n 1  $WRK/Library/F3c/CDT/YY_${filename}_sense.cdt  > $WRK/Library/E8/SSxphyloP30way_"${filename}"_sense.cdt
paste $WRK/Library/E8/YY_${filename}_sense_reorder.cdt $WRK/Library/E8/CDT/hg38.phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
  yorfn = $1
  name = $2
  total_fields = NF
  score_count = (total_fields / 2) - 2   # number of score columns
  line = yorfn OFS name

  # Loop through score1 .. scoreN
  for (i = 1; i <= score_count; i++) {
    c1 = 2 + i
    c2 = (total_fields / 2) + i + 1
    product = $(c1) * $(c2)
    if (product == 0) {
      sum = "-"
    } else {
      sum = product
    }
    line = line OFS sum
  }

  print line
}' >> $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt


head -n 1 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt  > $WRK/Library/E8/CDT/sensehead.cdt

tail -n +3 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt |  cat $WRK/Library/E8/CDT/sensehead.cdt - > $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt

rm $WRK/Library/E8/CDT/sensehead.cdt $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt $WRK/Library/E8/YY_${filename}_sense_reorder.cdt

perl $COMPOSITEFILTER $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/SSxphyloP30way_"${filename}"_sense_reorder.out

done

## times conservation score to YR pattern

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed   ; do
        filename=$(basename "$file" ".bed")

## reorder
awk -F'\t' '
NR==1 {
    print
    next
}
{
    # Copy original values
    for (i = 1; i <= NF; i++) {
        orig[i] = $i
    }

    # Only check from column 4 to end (3rd index = data starts)
    for (i = 4; i <= NF; i++) {
        if (orig[i] == 1) {
            $ (i - 1) = 1
        }
    }

    print
}' OFS='\t' $WRK/Library/F3c/CDT/RR_${filename}_sense.cdt > $WRK/Library/E8/RR_${filename}_sense_reorder.cdt


head -n 1  $WRK/Library/F3c/CDT/RR_${filename}_sense.cdt  > $WRK/Library/E8/SSxphyloP30way_"${filename}"_sense.cdt
paste $WRK/Library/E8/RR_${filename}_sense_reorder.cdt $WRK/Library/E8/CDT/hg38.phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
  yorfn = $1
  name = $2
  total_fields = NF
  score_count = (total_fields / 2) - 2   # number of score columns
  line = yorfn OFS name

  # Loop through score1 .. scoreN
  for (i = 1; i <= score_count; i++) {
    c1 = 2 + i
    c2 = (total_fields / 2) + i + 1
    product = $(c1) * $(c2)
    if (product == 0) {
      sum = "-"
    } else {
      sum = product
    }
    line = line OFS sum
  }

  print line
}' >> $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt


head -n 1 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt  > $WRK/Library/E8/CDT/sensehead.cdt

tail -n +3 $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt |  cat $WRK/Library/E8/CDT/sensehead.cdt - > $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt

rm $WRK/Library/E8/CDT/sensehead.cdt $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense.cdt $WRK/Library/E8/RR_${filename}_sense_reorder.cdt

perl $COMPOSITEFILTER $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/SSxphyloP30way_"${filename}"_sense_reorder.out

done



cd $WRK/Library/F3g
mkdir -p $WRK/Library/F3g/SCORES


for file in  $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed   ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/*_Pol2_*.bam --cpu 4 -5 -1  -M SCORES/Pol2_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $BAMDIR/*_Pol2_*.bam --cpu 4 -5 -1 -M SCORES/Pol2_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/Pol2_${filename}_150_read1_sense.cdt -o SCORES/Pol2_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/Pol2_${filename}_150_read1_anti.cdt -o SCORES/Pol2_${filename}_150_read1_anti_SCORES.out 
        rm SCORES/Pol2_${filename}_150_read1_sense.cdt SCORES/Pol2_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done

for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/Pol2_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "Pol2_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/Pol2_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "Pol2_${filename}_read1_anti_temp.txt"

   paste SCORES/Pol2_${filename}_150_read1_sense_SCORES.out SCORES/Pol2_${filename}_150_read1_anti_SCORES.out \
   | tail -n +2 | cut -f 2,4 | paste - "Pol2_${filename}_read1_sense_temp.txt" "Pol2_${filename}_read1_anti_temp.txt" \
   | awk '{
    OFS = "\t";
    total = $1 + $2;
    sum = $3 + $4;
    if (total == 0) total = 1;
    ratio = sum / total;
    print ratio
    }' | paste $file - > ${filename}_pol2score.bed
    paste SCORES/Pol2_${filename}_150_read1_sense_SCORES.out SCORES/Pol2_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "Pol2_${filename}_read1_sense_temp.txt" "Pol2_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        sametotal =  sum1;
        antitotal =  sum2;
        total =  sum1 + sum2;
        samepeak =  sum3;
        antipeak =  sum4;
        peaktotal =  sum3 + sum4;
        averagesame = sum3 / sum1;
        averageanti = sum4 / sum2;
        average = peaktotal / total;
        print sametotal, antitotal, samepeak, antipeak, averagesame, averageanti, average ;
    }' > ${filename}_Pol2phasescore.txt
    rm "Pol2_${filename}_read1_sense_temp.txt" "Pol2_${filename}_read1_anti_temp.txt" 
done

cat +1Nuc_YRWS_Pol2phasescore.txt +1Nuc_sameYR_lowWS_Pol2phasescore.txt +1Nuc_lowYR_sameWS_Pol2phasescore.txt +1Nuc_sameYR_antiWS_Pol2phasescore.txt +1Nuc_antiYR_sameWS_Pol2phasescore.txt +1Nuc_antiYR_antiWS_Pol2phasescore.txt +1Nuc_lessDNAencode_Pol2phasescore.txt >  Pol2Supphasescore.txt
rm +1Nuc_YRWS_Pol2phasescore.txt +1Nuc_sameYR_lowWS_Pol2phasescore.txt +1Nuc_lowYR_sameWS_Pol2phasescore.txt +1Nuc_sameYR_antiWS_Pol2phasescore.txt +1Nuc_antiYR_sameWS_Pol2phasescore.txt +1Nuc_antiYR_antiWS_Pol2phasescore.txt +1Nuc_lessDNAencode_Pol2phasescore.txt


for file in  $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam --cpu 4 -5 -1  -M SCORES/BI_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $Input --cpu 4 -5 -1 -M SCORES/BI_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/BI_${filename}_150_read1_sense.cdt -o SCORES/BI_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/BI_${filename}_150_read1_anti.cdt -o SCORES/BI_${filename}_150_read1_anti_SCORES.out
        rm SCORES/BI_${filename}_150_read1_sense.cdt SCORES/BI_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done

for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/BI_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "BI_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/BI_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "BI_${filename}_read1_anti_temp.txt"

     paste SCORES/BI_${filename}_150_read1_sense_SCORES.out SCORES/BI_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "BI_${filename}_read1_sense_temp.txt" "BI_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        sametotal =  sum1;
        antitotal =  sum2;
        total =  sum1 + sum2;
        samepeak =  sum3;
        antipeak =  sum4;
        peaktotal =  sum3 + sum4;
        averagesame = sum3 / sum1;
        averageanti = sum4 / sum2;
        average = peaktotal / total;
        print sametotal, antitotal, samepeak, antipeak, averagesame, averageanti, average ;
    }' > ${filename}_BIphasescore.txt
    rm "BI_${filename}_read1_sense_temp.txt" "BI_${filename}_read1_anti_temp.txt" 
done


cat +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt >  BIphasescore.txt
rm +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt

cd $WRK/Library/F3e
mkdir -p $WRK/Library/F3e/SCORES

for file in  $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam --cpu 4 -5 -1  -M SCORES/BI_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $Input --cpu 4 -5 -1 -M SCORES/BI_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/BI_${filename}_150_read1_sense.cdt -o SCORES/BI_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/BI_${filename}_150_read1_anti.cdt -o SCORES/BI_${filename}_150_read1_anti_SCORES.out
        rm SCORES/BI_${filename}_150_read1_sense.cdt SCORES/BI_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done

for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/BI_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "BI_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/BI_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "BI_${filename}_read1_anti_temp.txt"

     paste SCORES/BI_${filename}_150_read1_sense_SCORES.out SCORES/BI_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "BI_${filename}_read1_sense_temp.txt" "BI_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        sametotal =  sum1;
        antitotal =  sum2;
        total =  sum1 + sum2;
        samepeak =  sum3;
        antipeak =  sum4;
        peaktotal =  sum3 + sum4;
        averagesame = sum3 / sum1;
        averageanti = sum4 / sum2;
        average = peaktotal / total;
        print sametotal, antitotal, samepeak, antipeak, averagesame, averageanti, average ;
    }' > ${filename}_BIphasescore.txt
    rm "BI_${filename}_read1_sense_temp.txt" "BI_${filename}_read1_anti_temp.txt" 
done


cat +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt >  BIphasescore.txt
rm +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt

## salt

for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/K562_Input_Native100BI_rep1_hg38.bam --cpu 4 -5 -1  -M SCORES/100Sup_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $BAMDIR/K562_Input_Native100BI_rep1_hg38.bam --cpu 4 -5 -1 -M SCORES/100Sup_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/100Sup_${filename}_150_read1_sense.cdt -o SCORES/100Sup_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/100Sup_${filename}_150_read1_anti.cdt -o SCORES/100Sup_${filename}_150_read1_anti_SCORES.out 
        rm SCORES/100Sup_${filename}_150_read1_sense.cdt SCORES/100Sup_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done

for file in +1Nuc_antiYR_sameWS.bed  +1Nuc_lessDNAencode.bed  +1Nuc_YRWS.bed +1Nuc_sameYR_lowWS.bed  +1Nuc_lowYR_sameWS.bed +1Nuc_antiYR_antiWS.bed +1Nuc_sameYR_antiWS.bed ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/100Sup_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "100Sup_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/100Sup_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "100Sup_${filename}_read1_anti_temp.txt"

     paste SCORES/100Sup_${filename}_150_read1_sense_SCORES.out SCORES/100Sup_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "100Sup_${filename}_read1_sense_temp.txt" "100Sup_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        sametotal =  sum1;
        antitotal =  sum2;
        total =  sum1 + sum2;
        samepeak =  sum3;
        antipeak =  sum4;
        peaktotal =  sum3 + sum4;
        averagesame = sum3 / sum1;
        averageanti = sum4 / sum2;
        average = peaktotal / total;
        print sametotal, antitotal, samepeak, antipeak, averagesame, averageanti, average ;
    }' > ${filename}_100Supphasescore.txt
    rm "100Sup_${filename}_read1_sense_temp.txt" "100Sup_${filename}_read1_anti_temp.txt" 
done

cat +1Nuc_YRWS_100Supphasescore.txt +1Nuc_sameYR_lowWS_100Supphasescore.txt +1Nuc_lowYR_sameWS_100Supphasescore.txt +1Nuc_sameYR_antiWS_100Supphasescore.txt +1Nuc_antiYR_sameWS_100Supphasescore.txt +1Nuc_antiYR_antiWS_100Supphasescore.txt +1Nuc_lessDNAencode_100Supphasescore.txt >  100Supphasescore.txt
rm +1Nuc_YRWS_100Supphasescore.txt +1Nuc_sameYR_lowWS_100Supphasescore.txt +1Nuc_lowYR_sameWS_100Supphasescore.txt +1Nuc_sameYR_antiWS_100Supphasescore.txt +1Nuc_antiYR_sameWS_100Supphasescore.txt +1Nuc_antiYR_antiWS_100Supphasescore.txt +1Nuc_lessDNAencode_100Supphasescore.txt

for file in  $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/K562_Input_Native1000BI_rep1_hg38.bam    --cpu 4 -5 -1  -M SCORES/1000Sup_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $HisDIR/K562_1000Sup_NativeBI_merge_hg38.bam --cpu 4 -5 -1 -M SCORES/1000Sup_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/1000Sup_${filename}_150_read1_sense.cdt -o SCORES/1000Sup_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/1000Sup_${filename}_150_read1_anti.cdt -o SCORES/1000Sup_${filename}_150_read1_anti_SCORES.out 
        rm SCORES/1000Sup_${filename}_150_read1_sense.cdt SCORES/1000Sup_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done


for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/1000Sup_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "1000Sup-_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/1000Sup_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "1000Sup-_${filename}_read1_anti_temp.txt"

     paste SCORES/1000Sup_${filename}_150_read1_sense_SCORES.out SCORES/1000Sup_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "1000Sup-_${filename}_read1_sense_temp.txt" "1000Sup-_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        sametotal =  sum1;
        antitotal =  sum2;
        total =  sum1 + sum2;
        samepeak =  sum3;
        antipeak =  sum4;
        peaktotal =  sum3 + sum4;
        averagesame = sum3 / sum1;
        averageanti = sum4 / sum2;
        average = peaktotal / total;
        print sametotal, antitotal, samepeak, antipeak, averagesame, averageanti, average ;
    }' > ${filename}_1000Sup-phasescore.txt
    rm "1000Sup-_${filename}_read1_sense_temp.txt" "1000Sup-_${filename}_read1_anti_temp.txt" 
done


cat +1Nuc_YRWS_1000Sup-phasescore.txt +1Nuc_sameYR_lowWS_1000Sup-phasescore.txt +1Nuc_lowYR_sameWS_1000Sup-phasescore.txt +1Nuc_sameYR_antiWS_1000Sup-phasescore.txt +1Nuc_antiYR_sameWS_1000Sup-phasescore.txt +1Nuc_antiYR_antiWS_1000Sup-phasescore.txt +1Nuc_lessDNAencode_1000Sup-phasescore.txt >  1000Sup-phasescore.txt
rm +1Nuc_YRWS_1000Sup-phasescore.txt +1Nuc_sameYR_lowWS_1000Sup-phasescore.txt +1Nuc_lowYR_sameWS_1000Sup-phasescore.txt +1Nuc_sameYR_antiWS_1000Sup-phasescore.txt +1Nuc_antiYR_sameWS_1000Sup-phasescore.txt +1Nuc_antiYR_antiWS_1000Sup-phasescore.txt +1Nuc_lessDNAencode_1000Sup-phasescore.txt




