module load anaconda3
source activate virtualenv


### CHANGE ME
WRK=/Path/to/Title/
Reference=$WRK/0X_Bulk_Processing/Reference
###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/F3f ] || mkdir -p $WRK/Library/F3f
[ -d $WRK/Library/E10 ] || mkdir -p $WRK/Library/E10

cd $WRK/Library/F3f

awk -v file="+1Nuc_" '{ if ( $8 ~ /YRWS/  && $13> 0 && $12> 0 && $11> 0 && $10> 0 && $21 !~ /Divergent/  && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160    )  { print $0 > ( file "YRWS.bed")
     } else if ( $8 ~ /sameYR_lowWS/  && $11> 0 && $10> 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160  ) { print $0 > ( file "sameYR_lowWS.bed")
     } else if ( $8 ~ /sameYR_antiWS/ && $13< 0 && $12< 0 && $11> 0 && $10> 0 && $21 !~ /Divergent/  && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160   ) { print $0 > ( file "sameYR_antiWS.bed")
     } else if ( $8 ~ /antiYR_antiWS/ && $13< 0 && $12< 0 && $11< 0 && $10< 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160  ) { print $0 > ( file "antiYR_antiWS.bed")
     } else if ( $8 ~ /antiYR_sameWS/ && $13> 0 && $12> 0 && $11< 0 && $10< 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/ && $7 > 0 && $5 >=50 && $5 <=160   )  { print $0 > ( file "antiYR_sameWS.bed")
     } else if ( $8 ~ /lessDNAencode/ && $13 == "0" && $12 == "0"  && $11 == "0" && $10 == "0" && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/  && $7 > 0 && $5 >=50 && $5 <=160  )  { print $0 > ( file "lessDNAencode.bed")
     } else if ($8 ~ /lowYR_sameWS/ && $12> 0 && $13> 0 && $21 !~ /Divergent/ && $27 !~ /_noncodingTSS/  && $7 > 0 && $5 >=50 && $5 <=160  )  { print $0 > ( file "lowYR_sameWS.bed")
     }
        }' $WRK/04_plusoneNucleosome/Adj+1Nuc_Di_TSS_all_phase.bed
cat +1Nuc_YRWS.bed +1Nuc_sameYR_lowWS.bed +1Nuc_lowYR_sameWS.bed +1Nuc_sameYR_antiWS.bed +1Nuc_antiYR_sameWS.bed +1Nuc_antiYR_antiWS.bed  +1Nuc_lessDNAencode.bed   > phasescorethrethold.bed
rm +1Nuc_YRWS.bed +1Nuc_sameYR_lowWS.bed +1Nuc_lowYR_sameWS.bed +1Nuc_sameYR_antiWS.bed +1Nuc_antiYR_sameWS.bed  +1Nuc_antiYR_antiWS.bed +1Nuc_lessDNAencode.bed
echo -e "Site\tNuc_type\tPhase_score\tcore-Nuc\tOrientation\tcore-promoter\tTSS\tCpG_region\tHMM" > phasescorethrethold.csv
awk '{OFS="\t"} {print $4,$8,$7,$21,$27}' phasescorethrethold.bed > temp.txt

awk 'BEGIN{OFS="\t"}{
    # Print first two columns as is
    printf "%s\t%s\t", $1, $2;
    
    # Process remaining columns starting from column 3
    for(i=3;i<=NF;i++){
        n = split($i, arr, "_");
        for(j=1;j<=n;j++){
            # Use tab between fields, newline at end of row
            printf "%s%s", arr[j], (i==NF && j==n ? "\n" : "\t")
        }
    }
}' temp.txt > output.txt


awk '{OFS="\t"} {print $1,$2,$3,$15,$17,$20,$23,$24,$25}'  output.txt >> phasescorethrethold.csv
rm temp.txt  output.txt
python $WRK/bin/Correlation.py -i phasescorethrethold.csv

cp -r enrichment_tables  $WRK/Library/E10


