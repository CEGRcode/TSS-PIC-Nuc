#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/GROcap_download-ENCODE.log.out
#SBATCH -e logs/GROcap_download-ENCODE.log.err

# Download capped GRO-cap data from ENCODE

# data/BAM
#   |--GRO-cap_ENCODE_merge_hg38.bam

### CHANGE ME
WRK=/Path/to/Title/00_Download_and_Preprocessing


THREADS=8
###

# Dependencies
# - java
# - samtools
# - wget

set -exo
module load samtools
cd $WRK

# Inputs and outputs
Annotation=$WRK/../Annotation

# Download ref file from other github repo
wget https://github.com/CEGRcode/NucleosomePhasing/data/RefPT-Other/hg38_TSS-ALL_GENCODEV47_SORT.bed
wget https://github.com/CEGRcode/NucleosomePhasing/data/RefPT-Krebs/BNase-Nucleosomes.bed
wget https://github.com/CEGRcode/NucleosomePhasing/data/RefPT-Other/CpGIslands.bed
svn export https://github.com/CEGRcode/NucleosomePhasing/trunk/hg38_files
mv hg38_files $WRK/../
mv *.bed $Annotation

## download BroadHMM annotation and lift over to hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz
gunzip wgEncodeBroadHmmK562HMM.bed.gz
./liftOver wgEncodeBroadHmmK562HMM.bed hg19ToHg38.over.chain wgEncodeBroadHmmK562HMM_hg38.bed unmapped.bed
## seperate into Promoter/ Enhancer and other

awk '{if ($4 ~ /Promoter/) print $0 > "EncodeBroadHmmK562_Promoter.bed" ; else if ($4 ~ /Enhancer/) print $0 > "EncodeBroadHmmK562_Enhancer.bed"; else if ($4 ~ /Insulator/) print $0 > "EncodeBroadHmmK562_Insulator.bed"; else if ( ($4 ~ /Heterochrom/) || ($4 ~ /Repetitive/) || ($4 ~ /Repressed/)) print $0 > "EncodeBroadHmmK562_Repressive.bed"; else if ($4 ~ /Txn/) print $0 > "EncodeBroadHmmK562_Transcription.bed"; else print $0 > "EncodeBroadHmmK562_Other.bed" }' wgEncodeBroadHmmK562HMM_hg38.bed

rm wgEncodeBroadHmmK562HMM_hg38.bed  unmapped.bed wgEncodeBroadHmmK562HMM.bed

mv *.bed $Annotation

# Download GENCODE v47 GTF (or GFF) — adjust URL to match the official GENCODE v47 file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gunzip gencode.v47.annotation.gtf.gz

# Filter for long noncoding RNA entries (this depends on how they’re annotated in the GTF),
# then extract TSS as a BED-like format:
awk '$3 == "transcript" && /long_noncoding_RNA/ {
    # fields: $1 = chrom, $4 = start, $5 = end, $7 = strand
    chrom=$1;
    strand=$7;
    # For + strand, TSS = start; for – strand, TSS = end
    if (strand == "+") {
        tss = $4;
    } else {
        tss = $5;
    }
    # name: transcript_id or gene_id (extract from attribute column)
    match($0, /transcript_id "[^"]+"/);
    name = substr($0, RSTART, RLENGTH);
    # remove the `transcript_id ` prefix quotes
    gsub(/transcript_id /, "", name);
    gsub(/"/, "", name);
    # BED is 0-based, half-open — adjust: start = tss-1, end = tss
    bedstart = tss - 1;
    bedend = tss;
    print chrom "\t" bedstart "\t" bedend "\t" name "\t" 0 "\t" strand;
}' gencode.v47.annotation.gtf > gencode_v47.long_noncoding_RNAs_TSS.bed


# Assume gencode.v47.annotation.gtf is uncompressed

awk '$3 == "gene" && /gene_type "tRNA"/ {
    chrom = $1
    start = $4 - 1     # convert to 0-based for BED
    end = $5           # BED end is noninclusive, but for genes it is typically fine
    strand = $7

    # get gene_id (or any ID you prefer) from the attributes column
    attr = $9
    gene_id = "."

    # extract gene_id
    match(attr, /gene_id "[^"]+"/)
    if (RSTART > 0) {
      gene_id = substr(attr, RSTART+9, RLENGTH-10)
    }

    print chrom "\t" start "\t" end "\t" gene_id "\t0\t" strand
}' gencode.v47.annotation.gtf > gencode.v47.tRNAs.bed


mv *.bed $Annotation


