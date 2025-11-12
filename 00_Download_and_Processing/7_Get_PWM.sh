#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory
module load anaconda3
source activate meme

### CHANGE ME
WRK=/Path/to/Title/
###

# Inputs and outputs
MDIR=$WRK/data/sample-MEME
CALL_RefPT=$WRK/05_Call_RefPT/
PWM=$WRK/05_Call_RefPT/PWM

# Create output directories if they don't exist
[ -d $CALL_RefPT ] || mkdir $CALL_RefPT
[ -d $PWM ] || mkdir $PWM

# Hardcode move of MEME files from PEGR workflow to raname as TF.meme.txt file
cp $MDIR/34218_ZFP91_HPA065325_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/ZFP91.meme.txt
cp $MDIR/36714_GATA1_HPA000232_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/GATA1.meme.txt
cp $MDIR/33918_WDR5_HPA047182_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/WDR5.meme.txt
cp $MDIR/32124_YY1_HPA001119_K562_-_-_-_BX_hg38.meme.txt $PWM/YY1.meme.txt
cp $MDIR/40080_E2F7_HPA064866_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/E2F7.meme.txt
cp $MDIR/34364_NFYC_HPA055011_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/NFYC.meme.txt
cp $MDIR/32120_Sp1_HPA001853_K562_-_-_-_BX_hg38.meme.txt $PWM/SP1.meme.txt
cp $MDIR/32122_USF1_HPA036233_K562_-_-_-_BX_hg38.meme.txt $PWM/USF1.meme.txt
cp $MDIR/32114_GABPA_HPA003258_K562_-_-_-_BX_hg38.meme.txt $PWM/GABPA.meme.txt
cp $MDIR/38477_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.meme.txt $PWM/FOXA1.meme.txt
cp $MDIR/32117_NRF1_HPA029329_K562_-_-_-_BX_hg38.meme.txt $PWM/NRF1.meme.txt