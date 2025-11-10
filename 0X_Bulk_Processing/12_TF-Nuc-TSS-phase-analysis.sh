module load anaconda3
source activate bioinfo
# Script for F2a and E5

### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/F4e ] || mkdir -p $WRK/Library/F4e
[ -d $WRK/Library/E15 ] || mkdir -p $WRK/Library/E15

cd $WRK/Library/F4e


#make to TF_M1_adjNuc.csv and TF_M1_nearestTSS.csv

python TFBS_Nuc.py

python TFBS_TSS.py

python plot_all_tf_phase_Nuc_pies.sh

python plot_all_tf_phase_TSS_pies.sh

for file in *_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    filename=$(basename "$file" "_1bp.bed")
    mkdir -p "${TF}"/plot
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  ${TF}/AdjNuc.bed -o ${TF}/AdjNuc_1000bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  ${TF}/TSS.bed -o ${TF}/TSS_1000bp.bed
    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_+1Nuc ${TF}/${TF}_M1.bed ${TF}/AdjNuc_1000bp.bed
    perl $COMPOSITE ${TF}_M1_+1Nuc_sense.cdt ${TF}_M1_+1Nuc_sense
    perl $COMPOSITE ${TF}_M1_+1Nuc_anti.cdt ${TF}_M1_+1Nuc_anti
    tail -1 ${TF}_M1_+1Nuc_sense | cat ${TF}_M1_+1Nuc_anti - > ${TF}/plot/${TF}_M1_+1Nuc.out
    rm  ${TF}_M1_+1Nuc_sense ${TF}_M1_+1Nuc_anti ${TF}_M1_+1Nuc_*.cdt

    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_TSS ${TF}/${TF}_M1.bed ${TF}/TSS_1000bp.bed
    perl $COMPOSITE ${TF}_M1_TSS_sense.cdt ${TF}_M1_TSS_sense
    perl $COMPOSITE ${TF}_M1_TSS_anti.cdt ${TF}_M1_TSS_anti
    tail -1 ${TF}_M1_TSS_sense | cat ${TF}_M1_TSS_anti - > ${TF}/plot/${TF}_M1_TSS.out
    rm  ${TF}_M1_TSS_sense ${TF}_M1_TSS_anti ${TF}_M1_TSS_*.cdt

    rm ${TF}/AdjNuc_1000bp.bed
    rm ${TF}/TSS_1000bp.bed
done

for file in *_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    TFBAM=../BAM/K562_${TF}_BX_rep1_hg38.bam
    FACTOR=`grep 'Scaling factor' ../NormalizationFactors/K562_${TF}_BX_rep1_hg38_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`
    mkdir -p "$TF"/plot
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 ${TF}/nearestTSS_${TF}_M1_same-oppo.bed -o nearestTSS_${TF}_M1_same-oppo_1000bp.bed
    Ref=nearestTSS
    java -jar $SCRIPTMANAGER peak-analysis peak-align-ref --separate -o ${TF}_M1_${Ref} ${TF}/${TF}_M1.bed nearestTSS_${TF}_M1_same-oppo_1000bp.bed
    perl $COMPOSITE ${TF}_M1_${Ref}_sense.cdt ${TF}_M1_${Ref}_sense
    perl $COMPOSITE ${TF}_M1_${Ref}_anti.cdt ${TF}_M1_${Ref}_anti
    tail -1 ${TF}_M1_${Ref}_sense | cat ${TF}_M1_${Ref}_anti - > "$TF"/plot/${TF}_M1_${Ref}.out
    rm ${TF}_M1_${Ref}_sense.cdt
    rm ${TF}_M1_${Ref}_anti.cdt
    rm ${TF}_M1_${Ref}_anti
    rm ${TF}_M1_${Ref}_sense
    rm nearestTSS_${TF}_M1_same-oppo_1000bp.bed
    rm nearestTSS_${TF}_M1_same-oppo_500bp.bed
    rm ${TF}_M1_20bp.bed
done
cd $WRK/Library/E15

