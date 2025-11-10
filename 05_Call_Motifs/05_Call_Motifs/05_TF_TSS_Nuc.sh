module load gcc
#module load samtools
module load anaconda3
source activate bioinfo
# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/Path/to/Title/
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data/hg38_files/hg38.info.txt
BLACKLIST=$WRK/hg38_files/hg38-blacklist.bed

SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
MOTIFSCAN=$WRK/bin/scan_FASTA_for_motif_as_binary_string.py
PILEUPBW=$WRK//bin/pileup_BigWig_on_RefPT.py
COMPOSITE=$WRK/bin/sum_Col_CDT.pl 
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
BAMDIR=$WRK/data/BAM
NormDIR=$WRK/data/NormalizationFactors
Call_Motifs=$WRK/05_Call_Motifs/
cd $Call_Motifs


## collect motif, include TATA
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/SP1/SP1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/SP1/SP1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "SP1_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/YY1/YY1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/YY1/YY1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "YY1_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/GABPA/GABPA_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/GABPA/GABPA_experiment.bed | awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "GABPA_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/WDR5/WDR5_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/WDR5/WDR5_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "WDR5_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/E2F7_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/E2F7/E2F7_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "E2F7_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/FOXA1/FOXA1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/FOXA1/FOXA1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "FOXA1_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/GATA1/GATA1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/GATA1/GATA1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "GATA1_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/NRF1/NRF1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/NRF1/NRF1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "NRF1_Occupancy_1bp.bed" } }'
awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/USF1/USF1_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/USF1/USF1_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "USF1_Occupancy_1bp.bed" } }'

awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/NFYC/NFYC_motif1_Occupancy.bed | bedtools intersect -u -a - -b Chexmix/SP1/NFYC_experiment.bed |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "NFYC_Occupancy_originalstrand.bed" } }'
 awk '{ if ($6 == "+") { print $0 >> "NFY_originalstrand_+.bed" } else { print $0 >> "NFY_originalstrand_-.bed" } }'  NFYC_Occupancy_originalstrand.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"-",$7}' NFY_originalstrand_+.bed >  NFY_originalstrand_+_flip.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"+",$7}' NFY_originalstrand_-.bed  >  NFY_originalstrand_-_flip.bed
cat NFY_originalstrand_+_flip.bed  NFY_originalstrand_-_flip.bed | sort -k5,5nr  >  NFYC_Occupancy_1bp.bed
rm NFY_originalstrand_+.bed NFY_originalstrand_-.bed NFY_originalstrand_+_flip.bed  NFY_originalstrand_-_flip.bed  NFYC_Occupancy_originalstrand.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$8,$6,$7}' FIMO/ZFP91/ZFP91_motif1_Occupancy.bed | bedtools shift -g $Genome -p 8 -m -8  | bedtools intersect -u -a - -b Chexmix/ZFP91/ZFP91_experiment.bed  |  awk '{ if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ && $1 !~/chrM/) { print $0 > "ZFP91_Occupancy_originalstrand.bed" } }'
 awk '{ if ($6 == "+") { print $0 >> "ZFP91_Occupancy_1bp_originalstrand_+.bed" } else { print $0 >> "ZFP91_Occupancy_1bp_originalstrand_-.bed" } }' ZFP91_Occupancy_originalstrand.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"-",$7}' ZFP91_Occupancy_1bp_originalstrand_+.bed >  ZFP91_Occupancy_1bp_originalstrand_+_flip.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"+",$7}' ZFP91_Occupancy_1bp_originalstrand_-.bed >  ZFP91_Occupancy_1bp_originalstrand_-_flip.bed
cat ZFP91_Occupancy_1bp_originalstrand_+_flip.bed ZFP91_Occupancy_1bp_originalstrand_-_flip.bed | bedtools sort -i | uniq | sort -k5,5nr > ZFP91_Occupancy_1bp.bed
rm ZFP91_Occupancy_1bp_originalstrand_+.bed ZFP91_Occupancy_1bp_originalstrand_-.bed ZFP91_Occupancy_1bp_originalstrand_+_flip.bed ZFP91_Occupancy_1bp_originalstrand_-_flip.bed

cat ../03_core-promoter/FixedTATA_TSS.bed | \
awk '{OFS="\t"; print $1,$2,$3,$4,"0",$6,"TATA_M1"}' >  TATA_Occupancy_1bp.bed


for file in *_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    mkdir -p "$TF"

    # Sort and find closest TSS
    bedtools sort -i "${TF}_Occupancy_1bp.bed" | uniq | \
    bedtools closest -a - \
                     -b ../04_plusoneNucleosome/TSS_all_phase_adj+1Nuc_Di.bed \
                     -d -D a -t first | \
    sort -k43,43n > "${TF}/${TF}_M1_nearestTSS.bed"

    # Classify TSS proximity and orientation
    awk -v TF="$TF" '{
        if ($6 == $13 && $2 >= ($16 - 100) && $2 <= ($17 + 100)) {
            print $0 > (TF "_M1_same_TSS.bed");
        } else if ($2 >= ($16 - 100) && $2 <= ($17 + 100)) {
            print $0 > (TF "_M1_oppo_TSS.bed");
        } else if ($6 == $13) {
            print $0 > (TF "_M1_same_TSS_away.bed");
        } else {
            print $0 > (TF "_M1_oppo_TSS_away.bed");
        }
    }' "${TF}/${TF}_M1_nearestTSS.bed"

    # Generate TSS and TFBS and Nuc BED files
    cat "${TF}_M1_same_TSS.bed" "${TF}_M1_oppo_TSS.bed" | \
    awk 'BEGIN{OFS="\t"} { print $8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$30,$31,$32,$33,$34,$35,$36,$37 }' | \
    bedtools sort -i - | uniq > "${TF}/TSS.bed"

    cat "${TF}_M1_same_TSS.bed" "${TF}_M1_oppo_TSS.bed" | \
    awk 'BEGIN{OFS="\t"} { print $30,$31,$32,$33,$34,$35,$36,$37,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21 }' | \
    bedtools sort -i - | uniq > "${TF}/AdjNuc.bed"

    cat "${TF}_M1_same_TSS.bed" "${TF}_M1_oppo_TSS.bed" | \
    awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7 }' | \
    bedtools sort -i - | uniq > "${TF}/${TF}_M1.bed"

    # Closest TSS to TFBS
    bedtools closest -a "${TF}/TSS.bed" \
                     -b "${TF}/${TF}_M1.bed" \
                     -d -D a -t first | \
    sort -k30,30n > "${TF}/nearestTSS_${TF}_M1_sort.bed"

    # Closest Nuc to TFBS
    bedtools closest -a "${TF}/AdjNuc.bed" \
                     -b "${TF}/${TF}_M1.bed" \
                     -d -D a -t first | \
    sort -k30,30n > "${TF}/AdjNuc_${TF}_M1_sort.bed"

    # Classify same/opposite strand in closest TSS-TFBS pairs
    awk -v TF="$TF" '{
        if ($6 == $28) {
            print $0 > ("TSS_" TF "_same.bed");
        } else {
            print $0 > ("TSS_" TF "_oppo.bed");
        }
    }' "${TF}/nearestTSS_${TF}_M1_sort.bed"

    # Combine and move classified TSS files
    cat "TSS_${TF}_same.bed" "TSS_${TF}_oppo.bed" > "${TF}/nearestTSS_${TF}_M1_same-oppo.bed"
    mv "TSS_${TF}_same.bed" "TSS_${TF}_oppo.bed" "$TF/"
    mv "${TF}_M1_same_TSS.bed" "${TF}_M1_oppo_TSS.bed" "$TF/"
    # Clean up temporary away files
    rm -f "${TF}_M1_same_TSS_away.bed" "${TF}_M1_oppo_TSS_away.bed"
done


## annotation

cat */nearestTSS_*_M1_sort.bed | awk 'BEGIN{OFS="\t"} { print  $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > TF_M1_nearestTSS_sort.bed 
cat */AdjNuc_*_M1_sort.bed | awk 'BEGIN{OFS="\t"} { print  $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > TF_M1_AdjNuc_sort.bed 

    for file in TF_M1_nearestTSS_sort.bed  ; do
        filename=TF_M1_nearestTSS
        awk -v filename="$filename" '{
            last_digit = substr($30, length($30), 1)
            if ($30 >= 0) {
                print $0 > (filename "_phase" last_digit ".bed")
            } else {
                new_digit = 10 - last_digit
                print $0 > (filename "_phase" new_digit ".bed")
            }
        }' "$file"
    done


cat TF_M1_nearestTSS_phase10.bed TF_M1_nearestTSS_phase0.bed > TF_M1_nearestTSS_phase0_new.bed
rm TF_M1_nearestTSS_phase10.bed TF_M1_nearestTSS_phase0.bed
mv TF_M1_nearestTSS_phase0_new.bed TF_M1_nearestTSS_phase0.bed


for file in TF_M1_nearestTSS_phase*.bed; do
    filename=$(basename "$file" .bed)
    phase=$(echo "$filename" | cut -d "_" -f 4)

    # Split into "same" and "oppo"
    awk -v prefix="$filename" '{
        if ($6 == $13) {
            print > (prefix "_same.bed")
        } else {
            print > (prefix "_oppo.bed")
        }
    }' "$file"

    # Add metadata columns: phase and strand relationship
    awk -v phase="$phase" 'BEGIN{OFS="\t"} {
        print $0, phase, "same"
    }' "${filename}_same.bed" > "${filename}_same_annotated.bed"

    awk -v phase="$phase" 'BEGIN{OFS="\t"} {
        print $0, phase, "oppo"
    }' "${filename}_oppo.bed" > "${filename}_oppo_annotated.bed"

    rm ${filename}_oppo.bed ${filename}_same.bed
done


cat TF_M1_nearestTSS_phase*_*_annotated.bed > TF_M1_nearestTSS_annotated.bed
rm TF_M1_nearestTSS_phase*_*_annotated.bed
rm TF_M1_nearestTSS_phase*.bed


    for file in  TF_M1_AdjNuc_sort.bed  ; do
        filename=TF_M1_AdjNuc
        awk -v filename="$filename" '{
            last_digit = substr($30, length($30), 1)
            if ($30 >= 0) {
                print $0 > (filename "_phase" last_digit ".bed")
            } else {
                new_digit = 10 - last_digit
                print $0 > (filename "_phase" new_digit ".bed")
            }
        }' "$file"
    done

cat TF_M1_AdjNuc_phase10.bed TF_M1_AdjNuc_phase0.bed > TF_M1_AdjNuc_phase0_new.bed
rm TF_M1_AdjNuc_phase10.bed TF_M1_AdjNuc_phase0.bed
mv TF_M1_AdjNuc_phase0_new.bed TF_M1_AdjNuc_phase0.bed


for file in TF_M1_AdjNuc_phase*.bed; do
    filename=$(basename "$file" .bed)
    phase=$(echo "$filename" | cut -d "_" -f 4)

    # Split into "same" and "oppo"
    awk -v prefix="$filename" '{
        if ($6 == $13) {
            print > (prefix "_same.bed")
        } else {
            print > (prefix "_oppo.bed")
        }
    }' "$file"

    # Add metadata columns: phase and strand relationship
    awk -v phase="$phase" 'BEGIN{OFS="\t"} {
        print $0, phase, "same"
    }' "${filename}_same.bed" > "${filename}_same_annotated.bed"

    awk -v phase="$phase" 'BEGIN{OFS="\t"} {
        print $0, phase, "oppo"
    }' "${filename}_oppo.bed" > "${filename}_oppo_annotated.bed"

    rm ${filename}_oppo.bed ${filename}_same.bed
done

cat TF_M1_AdjNuc_phase*_*_annotated.bed > TF_M1_AdjNuc_annotated.bed
rm TF_M1_AdjNuc_phase*_*_annotated.bed
rm TF_M1_AdjNuc_phase*.bed

###########

## TF-specific analysis

### WDR5

cd WDR5

awk '{
    if ($6 == $28) {
        print $0 > ("TSS_same_WDR5_M1.bed");
    }
    else  {
        print $0 > ("TSS_oppo_WDR5_M1.bed");
    }
}' nearestTSS_WDR5_M1_sort.bed 

## +1Nuc

awk '{OFS="\t"} { print $15,$16,$17,$18,$19,$20,$21,$22}' TSS_same_WDR5_M1.bed > +1Nuc_TSS_same_WDR5_M1.bed

cd ..

## SP1 NFYC GABPA 

for file in SP1_Occupancy_1bp.bed NFYC_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    filename=$(basename "$file" "_1bp.bed")
    awk '{OFS="\t"} { print $15,$16,$17,$18,$19,$20,$21,$22}' ${TF}/nearestTSS_${TF}_M1_same-oppo.bed > ${TF}/+1Nuc_TSS_${TF}_M1.bed
done


for file in SP1_Occupancy_1bp.bed NFYC_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed  ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)

    awk -v TF=${TF} '{
    if ($6 == $28 && $30 < 0 && $30 > -100 ) {
        print $0 > ("TSS_same_" TF "_M1_up.bed");
    }
    else if ($37 < 0 && $30 > -100) {
        print $0 > ("TSS_oppo_" TF "_M1_up.bed");
    }}' ${TF}/nearestTSS_${TF}_M1_sort.bed

    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_oppo_${TF}_M1_up.bed > ${TF}/${TF}_M1_TSS_oppo_up100.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_up.bed > ${TF}/${TF}_M1_TSS_same_up100.bed

    mv TSS_oppo_${TF}_M1_up.bed TSS_same_${TF}_M1_up.bed ${TF}/
done



for file in  GABPA_Occupancy_1bp.bed NFYC_Occupancy_1bp.bed SP1_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    filename=$(basename "$file" "_1bp.bed")
    awk -v TF=${TF} '{
    if ($8 ~ /_TATA/) {
        print $0 > ("TSS_same_" TF "_M1_TATA.bed");
    }
    else if ($14 ~ /_noTATA/) {
        print $0 > ("TSS_same_" TF "_M1_noTATA.bed");
    }}'  ${TF}/TSS_same_${TF}_M1_up.bed

    mv TSS_same_${TF}_M1_TATA.bed TSS_same_${TF}_M1_noTATA.bed ${TF}/

    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' ${TF}/TSS_same_${TF}_M1_noTATA.bed > ${TF}/${TF}_M1_TSS_same_noTATA.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' ${TF}/TSS_same_${TF}_M1_TATA.bed > ${TF}/${TF}_M1_TSS_same_TATA.bed
done
 

## take core-promoter location and TF sites

for file in NFYC_Occupancy_1bp.bed  SP1_Occupancy_1bp.bed GABPA_Occupancy_1bp.bed ; do
    TF=$(basename "$file" ".bed" | cut -d "_" -f 1)
    bedtools sort -i ${TF}/TSS_same_${TF}_M1_TATA.bed > TSS_same_${TF}_M1_TATA_sort.bed
    bedtools intersect -u -a ../Nuc1/TSS_all_phase_adj+1Nuc_Di.bed -b  TSS_same_${TF}_M1_TATA_sort.bed | \
    awk '{OFS="\t"} { print $17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ${TF}/TATA_TSS_same_${TF}_M1_sort.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_TATA_sort.bed > ${TF}/${TF}_M1_TSS_same_TATA_sort.bed
    rm TSS_same_${TF}_M1_TATA_sort.bed
    bedtools sort -i ${TF}/TSS_same_${TF}_M1_noTATA.bed > TSS_same_${TF}_M1_noTATA_sort.bed
    bedtools intersect -u -a ../Nuc1/TSS_all_phase_adj+1Nuc_Di.bed -b  TSS_same_${TF}_M1_noTATA_sort.bed | \
    awk '{OFS="\t"} { print $17,$18,$19,$20,$21,$22,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ${TF}/noTATA_TSS_same_${TF}_M1_sort.bed
    awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_noTATA_sort.bed > ${TF}/${TF}_M1_TSS_same_noTATA_sort.bed
    rm TSS_same_${TF}_M1_noTATA_sort.bed
done


##combine SP1 and GABPA NFYC
cat GABPA/TSS_GABPA_same.bed GABPA/TSS_GABPA_oppo.bed  > TSS_GABPA.bed
cat SP1/TSS_SP1_same.bed SP1/TSS_SP1_oppo.bed  > TSS_SP1.bed
cat NFYC/TSS_NFYC_same.bed NFYC/TSS_NFYC_oppo.bed  > TSS_NFYC.bed

awk '{OFS="\t"} { print $1,$9,$10,$11,$12,$13}' TSS_GABPA.bed > NFR_TSS_GABPA.bed
awk '{OFS="\t"} { print $1,$9,$10,$11,$12,$13}' TSS_SP1.bed > NFR_TSS_SP1.bed
awk '{OFS="\t"} { print $1,$9,$10,$11,$12,$13}' TSS_NFYC.bed > NFR_TSS_NFYC.bed

bedtools intersect -u -a TSS_GABPA.bed -b NFR_TSS_SP1.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > GABPA_TSS_withSP1.bed
bedtools intersect -u -a TSS_NFYC.bed -b NFR_TSS_SP1.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > NFYC_TSS_withSP1.bed

bedtools intersect -v -a TSS_GABPA.bed -b NFR_TSS_SP1.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > GABPA_TSS_withoutSP1.bed
bedtools intersect -v -a TSS_NFYC.bed -b NFR_TSS_SP1.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > NFYC_TSS_withoutSP1.bed


bedtools intersect -u -a TSS_GABPA.bed -b NFR_TSS_NFYC.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > GABPA_TSS_withNFYC.bed
bedtools intersect -u -a TSS_SP1.bed -b NFR_TSS_NFYC.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > SP1_TSS_withNFYC.bed

bedtools intersect -v -a TSS_GABPA.bed -b NFR_TSS_NFYC.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > GABPA_TSS_withoutNFYC.bed
bedtools intersect -v -a TSS_SP1.bed -b NFR_TSS_NFYC.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > SP1_TSS_withoutNFYC.bed


bedtools intersect -u -a TSS_NFYC.bed -b NFR_TSS_GABPA.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > NFYC_TSS_withGABPA.bed
bedtools intersect -u -a TSS_SP1.bed -b NFR_TSS_GABPA.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > SP1_TSS_withGABPA.bed

bedtools intersect -v -a TSS_NFYC.bed -b NFR_TSS_GABPA.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > NFYC_TSS_withoutGABPA.bed
bedtools intersect -v -a TSS_SP1.bed -b NFR_TSS_GABPA.bed |awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' > SP1_TSS_withoutGABPA.bed

rm TSS_GABPA.bed TSS_NFYC.bed  TSS_SP1.bed  NFR_TSS_SP1.bed NFR_TSS_NFYC.bed NFR_TSS_GABPA.bed 

mkdir -p 3MOTIF 

mv *with* 3MOTIF 

cd YY1
TF=YY1

awk '{
    if ($6 == $28 && $30 <= 50 && $30 >= -50 ) {
        print $0 > "TSS_same_YY1_M1_overlap.bed";
    }
    else if ($6 != $28 && $30 <= 50 && $30 >= -50) {
        print $0 > "TSS_oppo_YY1_M1_overlap.bed";
    }
}' nearestTSS_${TF}_M1_sort.bed

wc -l nearestTSS_${TF}_M1_sort.bed
wc -l TSS_same_YY1_M1_overlap.bed
wc -l TSS_oppo_YY1_M1_overlap.bed
    1108 nearestTSS_YY1_M1_sort.bed
     531 TSS_same_YY1_M1_overlap.bed
     281 TSS_oppo_YY1_M1_overlap.bed
awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_oppo_${TF}_M1_overlap.bed > ${TF}_M1_TSS_oppo_overlap.bed
awk '{OFS="\t"} { print $23,$24,$25,$26,$27,$28,$29,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$30}' TSS_same_${TF}_M1_overlap.bed > ${TF}_M1_TSS_same_overlap.bed

awk '{
    if ($30 <= 5 && $30 >= 0 ) {
        print $0 > "YY1_M1_TSS_same_1.bed";
    }
    else if ($30 <= 18 && $30 >= 5 ) {
        print $0 > "YY1_M1_TSS_same_2.bed";
    }
    else if ($30 >= 18) {
        print $0 > "YY1_M1_TSS_same_3.bed";
    }
}' ${TF}_M1_TSS_same_overlap.bed 






























