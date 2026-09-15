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
[ -d $WRK/Library/F5E ] || mkdir -p $WRK/Library/F5E
[ -d $WRK/Library/E14 ] || mkdir -p $WRK/Library/E14

cd $WRK/Library/F5E

#make to TF_M1_adjNuc.csv and TF_M1_nearestTSS.csv

echo -e "Site\tTFBS\tDistance_Nuc\tTFBS_phase\tTFBS_strand" > TF_M1_adjNuc.csv
awk '{OFS="\t"} {print $4,$7,$30,$31,$32}' $WRK/05_Call_RefPT/TF_M1_AdjNuc_annotated.bed >> TF_M1_adjNuc.csv

python $WRK/bin/TFBS_Nuc.py TF_M1_adjNuc.csv

echo -e "Site\tTFBS\tDistance_TSS\tTFBS_phase\tTFBS_strand" > TF_M1_nearestTSS.csv
awk '{OFS="\t"} {print $4,$7,$30,$31,$32}' $WRK/05_Call_RefPT/TF_M1_nearestTSS_annotated.bed >> TF_M1_nearestTSS.csv

python $WRK/bin/TFBS_TSS.py TF_M1_nearestTSS.csv

python $WRK/bin/plot_all_tf_phase_Nuc_pies.sh TFBS_phase_skew_ratios_Nuc.tsv TFBS_phase_counts_Nuc

python $WRK/bin/plot_all_tf_phase_TSS_pies.sh TFBS_phase_skew_ratios_TSS.tsv TFBS_phase_counts_TSS

## for TSS-Nuc, same logic with TFBS


cd $WRK/Library/E14
conda deactivate 
conda activate bioinfo
for file in $WRK/05_Call_RefPT/*_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    filename=$(basename "$file" "_1bp.bed")
    mkdir -p "${TF}"
    mkdir -p "${TF}"/plot
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  $WRK/05_Call_RefPT/${TF}/AdjNuc.bed -o ${TF}/AdjNuc_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  $WRK/05_Call_RefPT/${TF}/TSS.bed -o ${TF}/TSS_1000bp.bed
    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_+1Nuc $WRK/05_Call_RefPT/${TF}/${TF}_M1.bed ${TF}/AdjNuc_1000bp.bed
    perl $COMPOSITE ${TF}_M1_+1Nuc_sense.cdt ${TF}_M1_+1Nuc_sense
    perl $COMPOSITE ${TF}_M1_+1Nuc_anti.cdt ${TF}_M1_+1Nuc_anti
    tail -1 ${TF}_M1_+1Nuc_sense | cat ${TF}_M1_+1Nuc_anti - > ${TF}/plot/${TF}_M1_+1Nuc.out
    rm  ${TF}_M1_+1Nuc_sense ${TF}_M1_+1Nuc_anti ${TF}_M1_+1Nuc_*.cdt

    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_TSS $WRK/05_Call_RefPT/${TF}/${TF}_M1.bed ${TF}/TSS_1000bp.bed
    perl $COMPOSITE ${TF}_M1_TSS_sense.cdt ${TF}_M1_TSS_sense
    perl $COMPOSITE ${TF}_M1_TSS_anti.cdt ${TF}_M1_TSS_anti
    tail -1 ${TF}_M1_TSS_sense | cat ${TF}_M1_TSS_anti - > ${TF}/plot/${TF}_M1_TSS.out
    rm  ${TF}_M1_TSS_sense ${TF}_M1_TSS_anti ${TF}_M1_TSS_*.cdt

    rm ${TF}/AdjNuc_1000bp.bed
    rm ${TF}/TSS_1000bp.bed
done

cp -r $WRK/Library/F5E/TF-Nuc_pie_charts  $WRK/Library/E14
cp -r $WRK/Library/F5E/TF-TSS_pie_charts  $WRK/Library/E14
cp $WRK/Library/F4e/TFBS_strand_bias_ratios.tsv $WRK/Library/E14

cd $WRK/Library/E14a
mkdir -p Heatmap 
for file in WDR5_Occupancy_1bp.bed   YY1_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    Ref=nearestTSS
    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_${Ref} $Reference/${TF}_M1_20bp.bed $Reference/nearestTSS_${TF}_M1_same-oppo_500bp.bed
    java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue ${TF}_M1_${Ref}_sense.cdt -o ${TF}_M1_${Ref}_sense.png
    java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red ${TF}_M1_${Ref}_anti.cdt -o ${TF}_M1_${Ref}_anti.png
    java -jar $SCRIPTMANAGER figure-generation merge-heatmap ${TF}_M1_${Ref}_sense.png ${TF}_M1_${Ref}_anti.png -o ${TF}_M1_${Ref}_merge.png
    java -jar $SCRIPTMANAGER figure-generation label-heatmap ${TF}_M1_${Ref}_merge.png -f 20 -l -250 -m 0 -r 250 -o Heatmap/${TF}_M1_${Ref}_merge.svg
    rm ${TF}_M1_${Ref}_*.png ${TF}_M1_${Ref}_*.cdt
done



