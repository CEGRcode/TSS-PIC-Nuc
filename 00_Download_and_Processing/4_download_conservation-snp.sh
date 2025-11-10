#!/bin/bash

### CHANGE ME
WRK=/Path/to/Title/00_Download_and_Preprocessing
###
bigWigToWig=$WRK/../bin/bigWigToWig
BWTOBG=$WRK/../bin/convert_wig_to_bedgraph.py
BBTOBED=$WRK/../bin/bigBedToBed
BGTOBW=$WRK/../bin/bedGraphToBigWig

CDIR=$WRK/../data/Conservation-SNP
CHRSZ=$WRK/../data/hg38_files/hg38.chrom.sizes

mkdir -p $CDIR
cd $CDIR

# Conservation processing - Download
wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP30way/hg38.phyloP30way.bw
mv hg38.phyloP30way.bw $CDIR/
#ls -l ~/Downloads/bigWigToWig
#chmod +x ~/Downloads/bigWigToWig 
$bigWigToWig hg38.phyloP30way.bw hg38.phyloP30way.wig

#python $BWTOBG -i hg38.phyloP30way.wig -o phg38.phyloP30way.bedgraph
# SNP processing - download and expand to standard BED format
wget -c http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb
#ls -l ~/Downloads/bigBedToBed
#chmod +x ~/Downloads/bigBedToBed    
$BBTOBED dbSnp153.bb dbSnp153.bed

# SNP processsing - split by snv and others (ins,del,delins,mnv)
awk '{OFS="\t"}{FS="\t"}{
 if ($14=="snv") {
   print $0 > "dbSnp153_snv.bed";
 } 
}' dbSnp153.bed


# Rebuild each filtered BED file as a BigWig file (for BigWig-style pileup)
for FILE in dbSnp153_snv.bed;
do
	BASE=`basename $FILE ".bed"`
	echo $BASE
	continue

	# Make a smaller BED
	# bedtools intersect -sorted -a $FILE -b <(bedtools sorted -i )

	# Assume all notated for "+" strand, also stripping non-standard chr
	awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3}' $FILE \
		| grep -v '_fix' |grep -v '_random' | grep -v 'chrUn_' | grep -v '_alt' \
		| uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$1}' > $BASE.bedGraph
	$BGTOBW $BASE.bedGraph $CHRSZ $BASE.bw

	# Clean-up
	 rm $BASE.bedGraph
done
