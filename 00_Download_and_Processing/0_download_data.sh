#!/bin/bash

# Pull FASTQ, BAM, and MEME files from PEGR
# This will need to be updated to run off of SRR samples before manuscript submission.

### CHANGE ME
WRK=/path/to/Title/00_Download_and_Preprocessing
###

# Dependencies
# - python
# - bioblend
# - pug

MYSAMPLES=$WRK/sampleIDs.txt
MEMESAMPLES=$WRK/MEMEsampleIDs.txt
FDIR=../data/FASTQ
BDIR=../data/sample-BAM
MDIR=../data/sample-MEME

[ -d $FDIR ] || mkdir $FDIR
[ -d $BDIR ] || mkdir $BDIR
[ -d $MDIR ] || mkdir -p $MDIR

# Download FASTQ
pug-fastq $MYSAMPLES
mv *.fastq.gz $FDIR

# Download Filtered BAM
 pug-bam $MYSAMPLES -b hg38
 mv *.bam $BDIR

# Download MEME files
pug-meme $MEMESAMPLES -b hg38
mv *.meme.txt $MDIR

# Download previous Benzonase-seq and ChIP-exo v6 data:
#!/usr/bin/env bash
set -uo pipefail
OUT=fastq; mkdir -p "$OUT"

while read -r GSM LABEL; do
  [[ -z "${GSM:-}" ]] && continue

  # GSM -> SRA internal id(s) -> SRR run accession(s)
  IDS=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${GSM}" | grep -oP '(?<=<Id>)\d+(?=</Id>)'); sleep 0.4
  for ID in $IDS; do
    for SRR in $(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${ID}&rettype=runinfo&retmode=csv" | tail -n +2 | cut -d, -f1 | grep '^[SED]RR'); do
      sleep 0.4
      # SRR -> ENA FASTQ url(s) (works for single- and paired-end)
      for URL in $(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | tr ';' '\n'); do
        FN=$(basename "$URL")
        wget -c -O "${OUT}/${LABEL}_${FN}" "https://${URL}" || echo "FAILED: $GSM $LABEL $SRR" >> failed.txt
      done
    done
  done
done <<'EOF'
GSM8607769	38477_K562_FOXA1_rep1-1
GSM8607773 33925_K562_IgG_rep1-1
GSM8607774 33963_K562_IgG_rep2-1
GSM8607775 34031_K562_IgG_rep3-1
GSM8607776 34055_K562_IgG_rep4-1
GSM8607777 34175_K562_IgG_rep5-1
GSM8607778 34469_K562_IgG_rep6-1
GSM8607779 34668_K562_IgG_rep7-1
GSM8607780 36713_K562_IgG_rep8-1
GSM8607781 36900_K562_IgG_rep9-1
GSM8607782 37449_K562_IgG_rep10-1
GSM8248142 28385_50U-3min_1-1
GSM8248143 28386_50U-3min_2-1
GSM8248144 28393_50U-30min_1-1
GSM8248145 28394_50U-30min_2-1
GSM8248146 25769_50U-10min_1-1
GSM8248147 25976_50U-10min_2-1
GSM8248148 28326_50U-10min_3-1
GSM8248149 28389_50U-10min_4-1
GSM8248150 28390_50U-10min_5-1
GSM8248151 28468_50U-10min_6-1
GSM8248152 28469_50U-10min_7-1
GSM8248153 28471_50U-10min_8-1
GSM8248154 28472_50U-10min_9-1
GSM8248155 28473_50U-10min_10-1
GSM8248156 28820_50U-10min_11-1
GSM8248157 33023_50U-10min_12-1
GSM8248158 33159_50U-10min_12-2
GSM8248159 33024_50U-10min_13-1
GSM8248160 33160_50U-10min_13-2
GSM8248161 33025_50U-10min_14-1
GSM8248162 33161_50U-10min_14-2
GSM8248163 33026_50U-10min_15-1
GSM8248164 33162_50U-10min_15-2
GSM8248165 33027_50U-10min_16-1
GSM8248166 33163_50U-10min_16-2
GSM8248167 33028_50U-10min_17-1
GSM8248168 33164_50U-10min_17-2
GSM8248169 33646_50U-10min_18-1
GSM8248170 33647_50U-10min_19-1
GSM8248171 33648_50U-10min_20-1
GSM8248206 25861_H3K4me3_1-1
GSM8248207	25869_H3K4me3_2-1
GSM8248208	25965_H3K4me3_3-1
GSM8248209	25972_H3K4me3_4-1
GSM8248210	28805_H3K4me3_5-1
GSM8248211	28809_H3K4me3_6-1

mv *.fastq.gz $FDIR
