# Genome-wide rotational and translational setting of transcription factors with nucleosomes

### Haining Chen<sup>&1</sup>, Olivia W. Lang<sup>1</sup>, William K. M. Lai<sup>1</sup>, B. Franklin Pugh<sup>1</sup>

<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

### Correspondence:fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [GSE266547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266547) [GSEXXXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX)

## Abstract
How transcription initiation complexes assemble on native promoters remains poorly understood. Here we establish global rules governing assembly of the RNA polymerase (Pol) II initiation complex in human cells via integration of transcription start sites, TATA boxes, +1 nucleosome rotational phasing, transcription factor binding. Our findings support an integrated model in which three promoter features—upstream bidirectional sequence-specific transcription factors, directional core promoter elements, and a DNA-encoded rotationally phased +1 nucleosome—work together to concentrate and orient general transcription factors (TFIIA and TFIID) near core promoters. This arrangement directs TFIIA and TBP to the −30 region through a mechanism consistent with established biochemical data. This positions Pol II to conduct a tightly focused search for the optimal initiator dinucleotide. Pol II then initiates transcription and pauses at the +1 nucleosome, without disrupting its rotational unless the nucleosome itself is unstable. Together, these findings define key organizational principles that govern preinitiation complex assembly, transcription initiation, and pausing within the native chromatin landscape of promoters.

## Directions
To recreate the figures for this manuscript, please execute the scripts in each directory in numerical order. Each directory's README includes more specific details on execution. To be more explicit, run the scripts in each directory in the following order: `00_Download_and_Preprocessing`, `01_Run_GenoPipe`, `02_TSS_NFR`, `03_core-promoter`, `04_plusonenucleosome`, `05_Call_Motifs`, `0X_Bulk_Processing`, and then finally `Library`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n bx -c bioconda -c conda-forge bedtools bowtie2 bwa cutadapt meme opencv pandas samtools scipy sra-tools wget pybigwig
```

For genetrack-executing script, a python2 environment needed to be created. The create command for that env is as follows:

```
conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy
```
For motif scanning and other python script, The create command for that env is as follows:

```
conda create -n virtualenv
pip install certifi contourpy cycler fonttools kiwisolver kneed matplotlib numpy packaging pandas patsy Pillow ply pyparsing PyQt5-sip pysam python-dateutil pytz scipy seaborn setuptools sip six statsmodels toml tornado tzdata wheel

```

## Table of Contents

### 00_Download_and_Preprocessing
Perform the preprocessing steps including alignment of raw sequencing data from both novel and previously published data.

### 01_Run_GenoPipe
Perform quality control for genetic background on these data by running GenoPipe on the aligned BAMs.

### 02_TSS_NFR
Call TSS sites based one PRO-cap RNA capped sites, define transcription activate region by determining +1 and -1 nucleosome relative to each TSS.

### 03_core-promoter
Define core-promoter region -- TSS upstream 30bp region -- is TATA or TATA-less

### 04_plusoneNucleosome
Call phased-aligned +1 nucleosome and group by dinucleotide encoding

### 05_Call_Motifs
Call TF binding motif 

### 0X_Bulk_Processing
With the BAM and BED files built from the scripts in the above directories, perform bulk read pileups for heatmaps and composites. Perform data quantification.

### Z_Figures
Copy/organize results from bulk processing into figure-specific directories corresponding to subfigures in the manuscript. 

### AI_files
all figures in paper

### data
Store large files to be globally accessed by the scripts in each directory

### bin
Generalized scripts and executables for global access by each of the numbered directories.
