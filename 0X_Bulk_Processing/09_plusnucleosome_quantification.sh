module load anaconda3
source activate bioinfo


### CHANGE ME
WRK=/Path/to/Title/
Reference=$WRK/0X_Bulk_Processing/Reference
BAMDIR=$WRK/data/BAM
###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library ] || mkdir -p $WRK/Library
[ -d $WRK/Library/E6 ] || mkdir -p $WRK/Library/E6
[ -d $WRK/Library/F3d ] || mkdir -p $WRK/Library/F3d
[ -d $WRK/Library/F3e ] || mkdir -p $WRK/Library/F3e
[ -d $WRK/Library/F3g ] || mkdir -p $WRK/Library/F3g
[ -d $WRK/Library/E8 ] || mkdir -p $WRK/Library/E8
[ -d $WRK/Library/ET3 ] || mkdir -p $WRK/Library/ET3

BAMDIR=$WRK/data/BAM
cd $WRK/Library/E6


for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed  ; do
        filename=$(basename "$file" ".bed")
        cp $WRK/04_plusoneNucleosome/SCORES/${filename}_YY_RR_WW_SS_means.out $WRK/Library/E6
done

cd $WRK/Library/F3d

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 Composites/"phyloP30way_${filename}.out" | \
    cut -f 427-576 | \
    awk '{
        OFS = "\t";
        print ($3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147)-($8+$18+$28+$38+$48+$59+$69+$79+$89+$100+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142);
    }'  > ${filename}_Conservation_means.out
done

cd $WRK/Library/E8

## times conservation score to WW pattern

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
}' OFS='\t' $WRK/04_plusoneNucleosome/SCORES/WW_${filename}_sense.cdt > WW_${filename}_sense_reorder.cdt
head -n 1  $WRK/04_plusoneNucleosome/SCORES/WW_${filename}_sense.cdt  > WWxphyloP30way_"${filename}"_sense.cdt
paste WW_${filename}_sense_reorder.cdt $WRK/04_plusoneNucleosome/SCORES/phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
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
}' >> WWxphyloP30way_"${filename}"_sense.cdt
head -n 1 WWxphyloP30way_"${filename}"_sense.cdt  > sensehead.cdt
tail -n +3 WWxphyloP30way_"${filename}"_sense.cdt |  cat sensehead.cdt - > $WRK/Library/E8/CDT/WWxphyloP30way_"${filename}"_sense_reorder.cdt
rm sensehead.cdt WWSxphyloP30way_"${filename}"_sense.cdt WW_${filename}_sense_reorder.cdt
perl $COMPOSITEFILTER $WRK/Library/E8/CDT/WWxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/WWxphyloP30way_"${filename}"_sense_reorder.out
done
## times conservation score to SS pattern
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
}' OFS='\t' $WRK/04_plusoneNucleosome/SCORES/SS_${filename}_sense.cdt > SS_${filename}_sense_reorder.cdt
head -n 1  $WRK/04_plusoneNucleosome/SCORES/SS_${filename}_sense.cdt  > SSxphyloP30way_"${filename}"_sense.cdt
paste SS_${filename}_sense_reorder.cdt $WRK/04_plusoneNucleosome/SCORES/phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
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
}' >> SSxphyloP30way_"${filename}"_sense.cdt
head -n 1 SSxphyloP30way_"${filename}"_sense.cdt  > sensehead.cdt
tail -n +3 SSxphyloP30way_"${filename}"_sense.cdt |  cat sensehead.cdt - > $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt
rm sensehead.cdt SSxphyloP30way_"${filename}"_sense.cdt SS_${filename}_sense_reorder.cdt
perl $COMPOSITEFILTER $WRK/Library/E8/CDT/SSxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/SSxphyloP30way_"${filename}"_sense_reorder.out
done

## times conservation score to WW pattern

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
}' OFS='\t' $WRK/04_plusoneNucleosome/SCORES/YY_${filename}_sense.cdt > YY_${filename}_sense_reorder.cdt
head -n 1  $WRK/04_plusoneNucleosome/SCORES/YY_${filename}_sense.cdt  > YYxphyloP30way_"${filename}"_sense.cdt
paste YY_${filename}_sense_reorder.cdt $WRK/04_plusoneNucleosome/SCORES/phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
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
}' >> YYxphyloP30way_"${filename}"_sense.cdt
head -n 1 YYxphyloP30way_"${filename}"_sense.cdt  > sensehead.cdt
tail -n +3 YYxphyloP30way_"${filename}"_sense.cdt |  cat sensehead.cdt - > $WRK/Library/E8/CDT/YYxphyloP30way_"${filename}"_sense_reorder.cdt
rm sensehead.cdt YYSxphyloP30way_"${filename}"_sense.cdt YY_${filename}_sense_reorder.cdt
perl $COMPOSITEFILTER $WRK/Library/E8/CDT/YYxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/YYxphyloP30way_"${filename}"_sense_reorder.out
done

## times conservation score to WW pattern

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
}' OFS='\t' $WRK/04_plusoneNucleosome/SCORES/RR_${filename}_sense.cdt > RR_${filename}_sense_reorder.cdt
head -n 1  $WRK/04_plusoneNucleosome/SCORES/RR_${filename}_sense.cdt  > RRxphyloP30way_"${filename}"_sense.cdt
paste RR_${filename}_sense_reorder.cdt $WRK/04_plusoneNucleosome/SCORES/phyloP30way_"${filename}".cdt | awk -F'\t' 'BEGIN { OFS = "\t" } {
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
}' >> RRxphyloP30way_"${filename}"_sense.cdt
head -n 1 RRxphyloP30way_"${filename}"_sense.cdt  > sensehead.cdt
tail -n +3 RRxphyloP30way_"${filename}"_sense.cdt |  cat sensehead.cdt - > $WRK/Library/E8/CDT/RRxphyloP30way_"${filename}"_sense_reorder.cdt
rm sensehead.cdt RRxphyloP30way_"${filename}"_sense.cdt RR_${filename}_sense_reorder.cdt
perl $COMPOSITEFILTER $WRK/Library/E8/CDT/RRxphyloP30way_"${filename}"_sense_reorder.cdt $WRK/Library/E8/Composites/RRxphyloP30way_"${filename}"_sense_reorder.out
done


cd $WRK/Library/F3g
mkdir -p $WRK/Library/F3g/SCORES

for file in  $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed   ; do
        filename=$(basename "$file" ".bed")
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $file -o ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_300bp.bed $BAMDIR/K562_PolII_BX_rep1_hg38.bam --cpu 4 -5 -1  -M SCORES/PolII_${filename}_read1
        rm  ${filename}_300bp.bed
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $BAMDIR/K562_PolII_BX_rep1_hg38.bam  --cpu 4 -5 -1 -M SCORES/PolII_${filename}_150_read1
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/PolII_${filename}_150_read1_sense.cdt -o SCORES/PolII_${filename}_150_read1_sense_SCORES.out 
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/PolII_${filename}_150_read1_anti.cdt -o SCORES/PolII_${filename}_150_read1_anti_SCORES.out 
        rm SCORES/PolII_${filename}_150_read1_sense.cdt SCORES/PolII_${filename}_150_read1_anti.cdt  ${filename}_150bp.bed
done

for file in $Reference/+1Nuc_antiYR_sameWS.bed  $Reference/+1Nuc_lessDNAencode.bed  $Reference/+1Nuc_YRWS.bed $Reference/+1Nuc_sameYR_lowWS.bed  $Reference/+1Nuc_lowYR_sameWS.bed $Reference/+1Nuc_antiYR_antiWS.bed $Reference/+1Nuc_sameYR_antiWS.bed  ; do
    filename=$(basename "$file" .bed)
    tail -n +2 "SCORES/PolII_${filename}_read1_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146+$5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147+$6+$16+$26+$36+$46+$57+$67+$77+$87+$97+$108+$118+$128+$138+$148)-($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141+$10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142+$1+$11+$21+$31+$41+$52+$62+$72+$82+$92+$103+$113+$123+$133);
    }'  > "PolII_${filename}_read1_sense_temp.txt"
    tail -n +2 "SCORES/PolII_${filename}_read1_anti.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144+$3+$13+$23+$33+$43+$54+$64+$74+$84+$94+$105+$115+$125+$135+$145+$4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146)-($7+$17+$27+$37+$47+$58+$68+$78+$88+$99+$109+$119+$129+$139+$8+$18+$28+$38+$48+$59+$69+$79+$89+$99+$110+$120+$130+$140+$9+$19+$29+$39+$49+$60+$70+$80+$90+$100+$111+$121+$131+$141);
    }'  > "PolII_${filename}_read1_anti_temp.txt"

    paste SCORES/PolII_${filename}_150_read1_sense_SCORES.out SCORES/PolII_${filename}_150_read1_anti_SCORES.out | tail -n +2 | cut -f 2,4 |  paste - "PolII_${filename}_read1_sense_temp.txt" "PolII_${filename}_read1_anti_temp.txt" | awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4;
    } END {
        total =  sum1 + sum2;
        peaktotal =  sum3 + sum4;
        average = peaktotal / total;
        print peaktotal, total, average ;
    }' > ${filename}_PolIIphasescore.txt
    rm "PolII_${filename}_read1_sense_temp.txt" "PolII_${filename}_read1_anti_temp.txt" 
done

cat +1Nuc_YRWS_PolIIphasescore.txt +1Nuc_sameYR_lowWS_PolIIphasescore.txt +1Nuc_lowYR_sameWS_PolIIphasescore.txt +1Nuc_sameYR_antiWS_PolIIphasescore.txt +1Nuc_antiYR_sameWS_PolIIphasescore.txt +1Nuc_antiYR_antiWS_PolIIphasescore.txt +1Nuc_lessDNAencode_PolIIphasescore.txt >  PolIISupphasescore.txt
rm +1Nuc_YRWS_PolIIphasescore.txt +1Nuc_sameYR_lowWS_PolIIphasescore.txt +1Nuc_lowYR_sameWS_PolIIphasescore.txt +1Nuc_sameYR_antiWS_PolIIphasescore.txt +1Nuc_antiYR_sameWS_PolIIphasescore.txt +1Nuc_antiYR_antiWS_PolIIphasescore.txt +1Nuc_lessDNAencode_PolIIphasescore.txt


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
        total =  sum1 + sum2;
        peaktotal =  sum3 + sum4;
        average = peaktotal / total;
        print peaktotal, total, average ;
    }' > ${filename}_BIphasescore.txt
    rm "BI_${filename}_read1_sense_temp.txt" "BI_${filename}_read1_anti_temp.txt" 
done

cat +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt >  BIphasescore.txt
rm +1Nuc_YRWS_BIphasescore.txt +1Nuc_sameYR_lowWS_BIphasescore.txt +1Nuc_lowYR_sameWS_BIphasescore.txt +1Nuc_sameYR_antiWS_BIphasescore.txt +1Nuc_antiYR_sameWS_BIphasescore.txt +1Nuc_antiYR_antiWS_BIphasescore.txt +1Nuc_lessDNAencode_BIphasescore.txt

cd $WRK/Library/F3e
mkdir -p $WRK/Library/F3e/SCORES

## salt treatment

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
        total =  sum1 + sum2;
        peaktotal =  sum3 + sum4;
        average = peaktotal / total;
        print peaktotal, total, average ;
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
        total =  sum1 + sum2;
        peaktotal =  sum3 + sum4;
        average = peaktotal / total;
        print peaktotal, total, average ;
    }' > ${filename}_1000Sup-phasescore.txt
    rm "1000Sup-_${filename}_read1_sense_temp.txt" "1000Sup-_${filename}_read1_anti_temp.txt" 
done

cat +1Nuc_YRWS_1000Sup-phasescore.txt +1Nuc_sameYR_lowWS_1000Sup-phasescore.txt +1Nuc_lowYR_sameWS_1000Sup-phasescore.txt +1Nuc_sameYR_antiWS_1000Sup-phasescore.txt +1Nuc_antiYR_sameWS_1000Sup-phasescore.txt +1Nuc_antiYR_antiWS_1000Sup-phasescore.txt +1Nuc_lessDNAencode_1000Sup-phasescore.txt >  1000Sup-phasescore.txt
rm +1Nuc_YRWS_1000Sup-phasescore.txt +1Nuc_sameYR_lowWS_1000Sup-phasescore.txt +1Nuc_lowYR_sameWS_1000Sup-phasescore.txt +1Nuc_sameYR_antiWS_1000Sup-phasescore.txt +1Nuc_antiYR_sameWS_1000Sup-phasescore.txt +1Nuc_antiYR_antiWS_1000Sup-phasescore.txt +1Nuc_lessDNAencode_1000Sup-phasescore.txt

cp $WRK/Library/F3g/BIphasescore.txt  $WRK/Library/F3e

### Sup table : calculate dinucleotide count in entire region ###
cd $WRK/Library/ET2
for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/RR_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144);
    }' | 
    > "${filename}_RR_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/YY_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149);
    }' | \
    > "${filename}_YY_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/SS_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147);
    }' | \
    > "${filename}_SS_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/WW_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142);
    }' | \
     > "${filename}_WW_score_temp.txt"

    paste ${filename}_RR_score_temp.txt ${filename}_YY_score_temp.txt ${filename}_SS_score_temp.txt ${filename}_WW_score_temp.txt | awk '{
        OFS = "\t";
        print $1/15, $2/15, $3/15, $3/14, ($1+$2+$3+$4)/59 ;
    }' > "${filename}_dinucleotide_ratio.txt"

    rm ${filename}_RR_score_temp.txt ${filename}_YY_score_temp.txt ${filename}_SS_score_temp.txt ${filename}_WW_score_temp.txt
done


for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed ; do
    filename=$(basename "$file" .bed)
    
    tail -n +2 "SCORES/RR_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($7+$17+$27+$37+$47+$58+$68+$78+$88+$98+$109+$119+$129+$139+$149);
    }' | 
    > "${filename}_RR_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/YY_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($2+$12+$22+$32+$42+$53+$63+$73+$83+$93+$104+$114+$124+$134+$144);
    }' | \
    > "${filename}_YY_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/SS_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($10+$20+$30+$40+$50+$61+$71+$81+$91+$101+$112+$122+$132+$142);
    }' | \
    > "${filename}_SS_score_temp.txt"

    tail -n +2 "$WRK/04_plusoneNucleosome/SCORES/WW_${filename}_sense.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($5+$15+$25+$35+$45+$56+$66+$76+$86+$96+$107+$117+$127+$137+$147);
    }' | \
     > "${filename}_WW_score_temp.txt"

    paste ${filename}_RR_score_temp.txt ${filename}_YY_score_temp.txt ${filename}_SS_score_temp.txt ${filename}_WW_score_temp.txt | awk '{
        OFS = "\t";
        print $1/15, $2/15, $3/14, $3/15, ($1+$2+$3+$4)/59 ;
    }' > "${filename}_dinucleotide_outofphase_ratio.txt"

    rm ${filename}_RR_score_temp.txt ${filename}_YY_score_temp.txt ${filename}_SS_score_temp.txt ${filename}_WW_score_temp.txt
done

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed  ; do
    filename=$(basename "$file" .bed)

    paste ${filename}_dinucleotide_ratio.txt \
          ${filename}_dinucleotide_outofphase_ratio.txt | \
    awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4; sum5 += $5; sum6 += $6; sum7 += $7; sum8 += $8; sum9 += $9; sum10 += $10;
        count++;
    } END {
        avg1 = sum1 / count;
        avg2 = sum2 / count;
        avg3 = sum3 / count;
        avg4 = sum4 / count;
        avg5 = sum5 / count;
        avg6 = sum6 / count;
        avg7 = sum7 / count;
        avg8 = sum8 / count;
        avg9 = sum9 / count;
        avg10 = sum10 / count;
        print avg1, avg6, avg2, avg7, avg3, avg8,avg4, avg9, avg5, avg10;
    }' > ${filename}_RYSW_in_out.out

    rm ${filename}_dinucleotide_ratio.txt ${filename}_dinucleotide_outofphase_ratio.txt
done

### Sup table : calculate dinucleotide count in entire region ###
cd $WRK/Library/ET3

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed   ; do
        filename=$(basename "$file" ".bed")

tail -n +2 $WRK/Library/E8/CDT/"SSxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146);
    }'  > "SSxphyloP30way_${filename}_score_in_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"SSxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141);
    }'  > "SSxphyloP30way_${filename}_score_outof_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"WWxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146);
    }'  > "WWxphyloP30way_${filename}_score_in_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"WWxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141);
    }'  > "WWxphyloP30way_${filename}_score_outof_phase.txt"


tail -n +2 $WRK/Library/E8/CDT/"YYxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146);
    }'  > "YYxphyloP30way_${filename}_score_in_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"YYxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141);
    }'  > "YYxphyloP30way_${filename}_score_outof_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"RRxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($4+$14+$24+$34+$44+$55+$65+$75+$85+$95+$106+$116+$126+$136+$146);
    }'  > "RRxphyloP30way_${filename}_score_in_phase.txt"

tail -n +2 $WRK/Library/E8/CDT/"RRxphyloP30way_"${filename}"_sense_reorder.cdt" | \
    cut -f 78-227 | \
    awk '{
        OFS = "\t";
        print ($9+$19+$29+$39+$49+$60+$70+$80+$90+$101+$111+$121+$131+$141);
    }'  > "RRxphyloP30way_${filename}_score_outof_phase.txt"

done

for file in $WRK/04_plusoneNucleosome/SCORES/*_Adj+1Nuc_TSS.bed; do
    filename=$(basename "$file" .bed)

    paste RRxphyloP30way_${filename}_score_in_phase.txt \
          RRxphyloP30way_${filename}_score_outof_phase.txt \
          YYxphyloP30way_${filename}_score_in_phase.txt \
          YYxphyloP30way_${filename}_score_outof_phase.txt \
          SSxphyloP30way_${filename}_score_in_phase.txt \
          SSxphyloP30way_${filename}_score_outof_phase.txt \
          WWxphyloP30way_${filename}_score_in_phase.txt \
          WWxphyloP30way_${filename}_score_outof_phase.txt \
          ${filename}_phyloP30way_score_in_phase.txt \
          ${filename}_phyloP30way_score_outof_phase.txt | \
    awk '{
        sum1 += $1; sum2 += $2; sum3 += $3; sum4 += $4; sum5 += $5; sum6 += $6; sum7 += $7; sum8 += $8; sum9 += $9; sum10 += $10;
        count++;
    } END {
        avg1 = sum1 / count;
        avg2 = sum2 / count;
        avg3 = sum3 / count;
        avg4 = sum4 / count;
        avg5 = sum5 / count;
        avg6 = sum6 / count;
        avg7 = sum7 / count;
        avg8 = sum8 / count;
        avg9 = sum9 / count;
        avg10 = sum10 / count;
        print avg1, avg2, avg3, avg4, avg5, avg6, avg7, avg8, avg9, avg10;
    }' > ${filename}_RYSWXphyloP30way_in_out.out

    rm *xphyloP30way_${filename}_score_*_phase.txt ${filename}_phyloP30way_score_*_phase.txt 
done
