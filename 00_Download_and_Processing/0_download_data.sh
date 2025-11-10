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
