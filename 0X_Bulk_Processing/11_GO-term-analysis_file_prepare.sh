WRK=/Path/to/Title/
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/S11B ] || mkdir -p $WRK/Library/S11B
[ -d $WRK/Library/S13B ] || mkdir -p $WRK/Library/S13B
cd $WRK/Library/S11B
cp $Reference/TSS_+1Nuc_*.bed $WRK/Library/S11B

cp $Reference/TSS_same_WDR5_M1.bed $WRK/Library/S13B
cp $Reference/TSS_all_phase_adj+1Nuc_Di.bed $WRK/Library/S13B
