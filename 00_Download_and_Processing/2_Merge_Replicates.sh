#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/Path/to/Title/00_Download_and_Preprocessing
###

module load anaconda3
source activate bioinfo
PICARD=$WRK/../bin/picard.jar
# Script shortcuts

cd $WRK
# Dowload K562_IgG_BX_merge_hg38.bam from GEO GSE267711 
# Dowload K562_FOXA1_BX_merge_hg38.bam from GEO GSE267711
# Download Benzonase-seq BNase-seq_50U-10min_merge_hg38.bam and H3K4me3 ChIP-seq BNase-ChIP_H3K4me3_merge_hg38.bam from GSE266547

[ -d $WRK/sample-BAM ] || mkdir $WRK/sample-BAM
[ -d $WRK/../data ] || mkdir $WRK/../data
 [-d $WRK/../data/BAM ] || mkdir $WRK/../data/BAM

mv *.bam $WRK/sample-BAM
# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# ChIP-exo rep1

cd $WRK/sample-BAM
cp 34601_RBBP5_A300-109A_K562_-_-_-_BX_hg38.bam K562_RBBP5_BX_rep1_hg38.bam
cp 33918_WDR5_HPA047182_K562_-_IMDM_-_BX_hg38.bam K562_WDR5_BX_rep1_hg38.bam
cp 36714_GATA1_HPA000232_K562_-_IMDM_-_BX_hg38.bam K562_GATA1_BX_rep1_hg38.bam
cp 33911_DR1_HPA055308_K562_-_IMDM_-_BX_hg38.bam K562_DR1_BX_rep1_hg38.bam

 java -jar $PICARD MergeSamFiles -I 33943_Med1_HPA052818_K562_-_IMDM_-_BX_hg38.bam\
                                -I 34204_Med1_HPA052818_K562_-_IMDM_-_BX_hg38.bam      \
                                -O K562_MED1_BX_rep1_hg38.bam 
                               
java -jar $PICARD MergeSamFiles -I 33945_MED12_HPA003185_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34206_MED12_HPA003185_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_MED12_BX_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 40049_E2F7_HPA064866_K562_-_IMDM_-_BX_hg38.bam \
                                -I 40080_E2F7_HPA064866_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_E2F7_BX_rep1_hg38.bam
java -jar $PICARD MergeSamFiles -I 32064_EP300_HPA003128_K562_-_-_-_BX_hg38.bam \
                                -I 32113_EP300_HPA003128_K562_-_-_-_BX_hg38.bam \
                                -O K562_EP300_BX_rep1_hg38.bam
java -jar $PICARD MergeSamFiles -I 33907_ERCC3_HPA046077_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34076_ERCC3_HPA046077_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_ERCC3_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32065_GABPA_HPA003258_K562_-_-_-_BX_hg38.bam \
                                -I 32114_GABPA_HPA003258_K562_-_-_-_BX_hg38.bam \
                                -O K562_GABPA_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 33902_GTF2A1_HPA000869_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34072_GTF2A1_HPA000869_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_GTF2A1_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32066_GTF2B_HPA061626_K562_-_-_-_BX_hg38.bam \
                                -I 32115_GTF2B_HPA061626_K562_-_-_-_BX_hg38.bam \
                                -O K562_GTF2B_BX_rep1_hg38.bam


java -jar $PICARD MergeSamFiles -I 33941_NELFA_HPA043931_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34202_NELFA_HPA043931_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_NELFA_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 35656_NFYC_HPA055011_K562_-_-_-_BX_hg38.bam \
                                -I 34364_NFYC_HPA055011_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34166_NFYC_HPA055011_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_NFYC_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32068_NRF1_HPA029329_K562_-_-_-_BX_hg38.bam \
                                -I 32117_NRF1_HPA029329_K562_-_-_-_BX_hg38.bam \
                                -O K562_NRF1_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32868_PolII_ab76123_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37479_PolII_ab76123_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_PolII_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32071_Sp1_HPA001853_K562_-_-_-_BX_hg38.bam \
                                -I 32120_Sp1_HPA001853_K562_-_-_-_BX_hg38.bam \
                                -O K562_SP1_BX_rep1_hg38.bam  

cp 33910_Taf1_HPA001075_K562_-_IMDM_-_BX_hg38.bam K562_TAF1_BX_rep1_hg38.bam                           

java -jar $PICARD MergeSamFiles -I 37571_Taf3_HPA066184_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37424_Taf3_HPA066184_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF3_BX_rep1_hg38.bam  


java -jar $PICARD MergeSamFiles -I 38045_TAF4B_HPA028937_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37539_TAF4B_HPA028937_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37354_TAF4B_HPA028937_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF4B_BX_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 40083_Taf9_HPA072658_K562_-_IMDM_-_BX_hg38.bam \
                                -I 40057_Taf9_HPA072658_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF9_BX_rep1_hg38.bam


java -jar $PICARD MergeSamFiles -I 38648_Taf11_HPA049127_K562_-_IMDM_-_BX_hg38.bam \
                                -I 38459_Taf11_HPA049127_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF11_BX_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 36982_Taf10_HPA004148_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37329_Taf10_HPA004148_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF10_BX_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 42801_Taf13_HPA044492_K562_-_IMDM_-_BX_hg38.bam \
                                -I 42672_Taf13_HPA044492_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF13_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 32662_TBP_hTBPserum_K562_-_-_-_BX_hg38.bam \
                                -I 32755_TBP_hTBPserum_K562_-_-_-_BX_hg38.bam \
                                -O K562_TBP_BX_rep1_hg38.bam  

java -jar $PICARD MergeSamFiles -I 32073_USF1_HPA036233_K562_-_-_-_BX.bam \
                                -I 32122_USF1_HPA036233_K562_-_-_-_BX.bam \
                                -O K562_USF1_BX_rep1_hg38.bam  

java -jar $PICARD MergeSamFiles -I 32075_YY1_HPA001119_K562_-_-_-_BX_hg38.bam \
                                -I 32124_YY1_HPA001119_K562_-_-_-_BX_hg38.bam \
                                -O K562_YY1_BX_rep1_hg38.bam  

java -jar $PICARD MergeSamFiles -I 34218_ZFP91_HPA065325_K562_-_IMDM_-_BX_hg38.bam \
                                -I 33961_ZFP91_HPA065325_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_ZFP91_BX_rep1_hg38.bam  

java -jar $PICARD MergeSamFiles -I 42132_Input_-_K562_-_IMDM_-_BI_hg38.bam \
                                -I 41941_Input_-_K562_-_IMDM_-_BI_hg38.bam \
                                -O K562_Input_Native100BI_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 41943_Input_-_K562_-_IMDM_-_BI_hg38.bam \
                                -I 42133_Input_-_K562_-_IMDM_-_BI_hg38.bam \
                                -O K562_Input_Native1000BI_rep1_hg38.bam 

java -jar $PICARD MergeSamFiles -I 33788_IgG_i5006_K562_-_-_Triptolide_BX_hg38.bam  \
                                -I 34427_IgG_i5006_K562_-_-_Triptolide_BX_hg38.bam \
                                -I 35225_IgG_i5006_K562_-_-_Triptolide_BX_hg38.bam  \
                                -I 36555_IgG_i5006_K562_-_IMDM_Triptolide_BX_hg38.bam \
                                -O TriptolideK562_IgG_BX_merge_hg38.bam

java -jar $PICARD MergeSamFiles -I 36557_IgG_i5006_K562_-_IMDM_DMSO_BX_hg38.bam  \
                                -I 34429_IgG_i5006_K562_-_-_DMSO_BX_hg38.bam \
                                -I 35227_IgG_i5006_K562_-_-_DMSO_BX_hg38.bam   \
                                -I 33821_IgG_i5006_K562_-_-_DMSO_BX_hg38.bam \
                                -I 33853_IgG_i5006_K562_-_-_DMSO_BX_hg38.bam  \
                                -O DMSOK562_IgG_BX_merge_hg38.bam
                                
cp 40520_GTF2A1_HPA000869_K562_-_IMDM_Triptolide_BX_hg38.bam TriptolideK562_GTF2A1_BX_rep1_hg38.bam
cp 40522_GTF2A1_HPA000869_K562_-_IMDM_DMSO_BX_hg38.bam DMSOK562_GTF2A1_BX_rep1_hg38.bam
cp 33780_TBP_hTBP_K562_-_-_Triptolide_BX_hg38.bam    TriptolideK562_TBP_BX_rep1_hg38.bam      
cp 33782_TBP_hTBP_K562_-_-_DMSO_BX_hg38.bam  DMSOK562_TBP_BX_rep1_hg38.bam
cp 33776_PolII_ab76123_K562_-_-_Triptolide_BX_hg38.bam    TriptolideK562_PolII_BX_rep1_hg38.bam
cp 33778_PolII_ab76123_K562_-_-_DMSO_BX_hg38.bam    DMSOK562_PolII_BX_rep1_hg38.bam
java -jar $PICARD MergeSamFiles -I 35295_TFIIB_HPA061626_K562_-_-_Triptolide_BX_hg38.bam    \
                                -I 33997_TFIIB_HPA061626_K562_-_-_Triptolide_BX_hg38.bam    \
                                -O TriptolideK562_TFIIB_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 35297_TFIIB_HPA061626_K562_-_-_DMSO_BX_hg38.bam      \
                                -I 33999_TFIIB_HPA061626_K562_-_-_DMSO_BX_hg38.bam    \
                                -O DMSOK562_TFIIB_BX_rep1_hg38.bam
# ChIP-exo rep2
## GTF2A2 and GTF2A1 both are subunits of TFIIA
cp 33903_GTF2A2_HPA056239_K562_-_IMDM_-_BX_hg38.bam  K562_GTF2A2_BX_rep1_hg38.bam
cp 40517_GTF2A1_HPA000869_K562_-_IMDM_-_BX_hg38.bam K562_GTF2A1_BX_rep2_hg38.bam

## DR1 and DRAP1 are same complex
cp 33912_DRAP1_HPA006790_K562_-_IMDM_-_BX_hg38.bam K562_DRAP1_BX_rep1_hg38.bam
 ## ERCC3 and GTF2H1 GTF2H2 are components of the TFIIH complex                      
java -jar $PICARD MergeSamFiles -I 33905_GTF2H1_HPA046660_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34074_GTF2H1_HPA046660_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_GTF2H1_BX_rep1_hg38.bam
java -jar $PICARD MergeSamFiles -I 33906_GTF2H2_HPA047001_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34075_GTF2H2_HPA047001_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_GTF2H2_BX_rep1_hg38.bam  

 ## NELFA and NELFB NELFE NELF are components of the NELF complex    
java -jar $PICARD MergeSamFiles -I 33942_NELFE_HPA046502_K562_-_IMDM_-_BX_hg38.bam \
                                -I 34203_NELFE_HPA046502_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_NELFE_BX_rep1_hg38.bam

java -jar $PICARD MergeSamFiles -I 40833_NELFB_HPA020259_K562_-_IMDM_-_BX_hg38.bam \
                                -I 40931_NELFB_HPA020259_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_NELFB_BX_rep1_hg38.bam 
## TAF4B and TAF4 from same subunit
java -jar $PICARD MergeSamFiles -I 36756_Taf4_HPA008599_K562_-_IMDM_-_BX_hg38.bam \
                                -I 37188_Taf4_HPA008599_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF4_BX_rep1_hg38.bam 
## TAF9B and TAF9 from the same subunit
java -jar $PICARD MergeSamFiles -I 38857_TAF9B_HPA045275_K562_-_IMDM_-_BX_hg38.bam \
                                -I 38552_TAF9B_HPA045275_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_TAF9B_BX_rep1_hg38.bam 
## real reaplicate 2 

java -jar $PICARD MergeSamFiles -I 34959_EP300_A300-358A_K562_-_IMDM_-_BX_hg38.bam \
                                -I 35104_EP300_A300-358A_K562_-_IMDM_-_BX_hg38.bam \
                                -O K562_EP300_BX_rep2_hg38.bam
cp 34504_GTF2H1_Antibody-1A10_K562_-_IMDM_-_BX_hg38.bam K562_GTF2H1_BX_rep2_hg38.bam
cp 34599_NRF1_3H1-s_K562_-_-_-_BX_hg38.bam K562_NRF1_BX_rep2_hg38.bam
cp 32869_PolII_ab76123_K562_-_IMDM_-_BX_hg38.bam  K562_PolII_BX_rep2_hg38.bam
cp 34600_RBBP5_A300-109A_K562_-_-_-_BX_hg38.bam K562_RBBP5_BX_rep2_hg38.bam
cp 34604_Sp1_HPA001853_K562_-_-_-_BX_hg38.bam K562_SP1_BX_rep2_hg38.bam
cp 32661_TBP_hTBPpurified_K562_-_-_-_BX_hg38.bam K562_TBP_BX_rep2_hg38.bam
cp 32099_USF1_1B8_K562_-_-_-_BX_hg38.bam K562_USF1_BX_rep2_hg38.bam
cp 34606_WDR5_HPA047182_K562_-_-_-_BX_hg38.bam K562_WDR5_BX_rep2_hg38.bam
cp 32097_YY1_1B2_K562_-_-_-_BX_hg38.bam K562_YY1_BX_rep2_hg38.bam
cp 32245_ZFP91_A303-245A_K562_-_-_-_BX_hg38.bam K562_ZFP91_BX_rep2_hg38.bam
cp 30699_GATA1_HPA000233_K562_-_-_Lysis-50Unuclease-10cycSonic-splintOligoswithout3primeddC-OnlyPost2ndLigationAmpure-52degreeAnnealing_BX_hg38.bam K562_GATA1_BX_rep2_hg38.bam
cp 34928_E2F7_A303-037A_K562_-_IMDM_-_BX_hg38.bam K562_E2F7_BX_rep2_hg38.bam
cp 34607_GABPA_HPA003258_K562_-_-_-_BX_hg38.bam K562_GABPA_BX_rep2_hg38.bam
cp 34617_GTF2B_HPA061626_K562_-_-_-_BX_hg38.bam K562_GTF2B_BX_rep2_hg38.bam
cp 35371_NELFE_sc-377052_K562_-_IMDM_-_BX_hg38.bam K562_NELFE_BX_rep2_hg38.bam

cp 41892_Input_-_K562_-_IMDM_-_BI_hg38.bam K562_Input_Native100BI_rep2_hg38.bam 
cp 41945_Input_-_K562_-_IMDM_-_BI_hg38.bam K562_Input_Native1000BI_rep2_hg38.bam

## will add TAF1 TAF3 NFYC replicate 2
#cp 34531_NFYC_NFYC-1A11_K562_-_IMDM_-_BX_hg38.bam K562_NFYC_BX_rep2_hg38.bam

java -jar $PICARD MergeSamFiles -I 33993_PolII_ab76123_K562_-_-_Triptolide_BX_hg38.bam \
                                -I 35291_PolII_ab76123_K562_-_-_Triptolide_BX_hg38.bam \
                                -O TriptolideK562_PolII_BX_rep2_hg38.bam
  
java -jar $PICARD MergeSamFiles -I 33995_PolII_ab76123_K562_-_-_DMSO_BX_hg38.bam   \
                                -I 35293_PolII_ab76123_K562_-_-_DMSO_BX_hg38.bam   \
                                -O DMSOK562_PolII_BX_rep2_hg38.bam


java -jar $PICARD MergeSamFiles -I 36545_TFIIB_HPA061626_K562_-_IMDM_Triptolide_BX_hg38.bam   \
                                -I 33807_TFIIB_HPA061626_K562_-_-_Triptolide_BX_hg38.bam      \
                                -O TriptolideK562_TFIIB_BX_rep2_hg38.bam

java -jar $PICARD MergeSamFiles -I 36547_TFIIB_HPA061626_K562_-_IMDM_DMSO_BX_hg38.bam      \
                                -I 33809_TFIIB_HPA061626_K562_-_-_DMSO_BX_hg38.bam               \
                                -O DMSOK562_TFIIB_BX_rep2_hg38.bam


cd $WRK/../data
mv sample-BAM/K562_*.bam BAM/
mv sample-BAM/TriptolideK562_*.bam sample-BAM/DMSOK562_*.bam BAM/

# Index set of BAM files
for FILE in  $WRK/../data/BAM/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done
