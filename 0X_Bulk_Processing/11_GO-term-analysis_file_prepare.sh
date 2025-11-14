WRK=/Path/to/Title/
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/E11 ] || mkdir -p $WRK/Library/E11
[ -d $WRK/Library/E14 ] || mkdir -p $WRK/Library/E14
cd $WRK/Library/E11
cp $Reference/TSS_+1Nuc_*.bed $WRK/Library/E11

cp $Reference/TSS_same_WDR5_M1.bed $WRK/Library/E14
cp $Reference/TSS_all_phase_adj+1Nuc_Di.bed $WRK/Library/E14
