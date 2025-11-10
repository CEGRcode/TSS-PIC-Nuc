
# data/BAM
#   |--U2OS_YY1KD_scRNA_rep1_hg38.bam
#   |--U2OS_Control_scRNA_rep1_hg38.bam

## creat STAR genome index

### CHANGE ME
WRK=/Path/to/Title/
STAR_INDEX=$WRK/data/STAR_INDEX
cd $WRK
module load anaconda3
conda install -c bioconda star
conda install -c bioconda/label/cf201901 star
## download index
cd $WRK/data

wget http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

gunzip *.gz

STAR  --runMode genomeGenerate --genomeDIR STAR_INDEX/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.115.gtf --runThreadN 16

rm Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
rm Homo_sapiens.GRCh38.115.gtf.gz

# Dependencies
# - bwa
# - java
# - samtools
# - SRA toolkit (fasterq-dump)

set -exo
module load anaconda3
module load samtools

# Inputs and outputs
BAMDIR=$WRK/data/BAM
GENOME=$WRK/data/hg38_files/hg38.fa

# Script shortcuts
PICARD=$WRK/bin/picard.jar

cd BAMDIR=$WRK/data/

[ -d FQ ] || mkdir FQ
[ -d BAM ] || mkdir BAM

# Download raw FASTQ and process

for SRR in "SRR18485977" "SRR18486002" ; do
    mkdir -p FQ
    fasterq-dump $SRR -O FQ --threads 8

    # ---- Filter and clip reads ----
    echo "Filtering and trimming reads with fastp..."
    fastp \
        -i FQ/${SRR}.fastq \
        -o FQ/${SRR}_filtered.fastq \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        --thread 8 \
        --html FQ/${SRR}_fastp_report.html \
        --json FQ/${SRR}_fastp_report.json
     rm FQ/${SRR}.fastq

    STAR \
    --runMode alignReads \
    --runThreadN 16 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$WRK/data/FQ/${SRR}_filtered.fastq" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$BAMDIR/${SRR}_" \

done

## rm

cp SRR18485977_*.bam U2OS_YY1KD_scRNA_rep1_hg38.bam
cp SRR18486002_*.bam U2OS_Control_scRNA_rep1_hg38.bam

rm SRR18485977_*.bam SRR18486002_*.bam
rm *.out
rm *.tab

