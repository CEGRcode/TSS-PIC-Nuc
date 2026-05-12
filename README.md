# Fine architecture and mechanism of PIC assembly at natural genomic promoters


### Haining Chen<sup>1</sup>, Olivia W. Lang<sup>1</sup>, William K. M. Lai<sup>1,2</sup>, B. Franklin Pugh<sup>1,2,3</sup>

<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA
<sup>2</sup>Co-senior author
<sup>3</sup>Correspondence:fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [GSE318610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318610)

## Abstract
How human transcription pre-initiation complexes (PICs) assemble at promoters and enhancers within their natural genomic context remains poorly understood. Little is known about the role of the +1 nucleosome including whether it is rotationally phased and whether such phasing is DNA-encoded. Little is known about how sequence-specific transcription factors (ssTFs) orchestrate PIC assembly through TFIID, TFIIA, the +1 nucleosome, and whether the helical DNA structure imparts constraints on assembly. Here we use single-bp resolution ChIP-exo, coupled to same-molecule measurements of nucleosome rotational phasing, to reveal the molecular architecture and mechanistic steps governing natural PIC assembly. We find that +1 nucleosomes have robust DNA-encoded rotational phasing that engages TFIID (TAF3) and other ssTFs. This, plus ssTF (SP1, GABPA and NFYC) interactions with TFIID (TAF4) and TFIIA concentrate and orient the PIC within selected DNA gyres as predicted by cryoEM structures. Such placement positions Pol II to conduct a tightly focused search for the optimal initiator. When Pol II initiates and pauses at the +1 nucleosome it disrupts rotational phasing but only where the +1 nucleosome is biochemically unstable. Together, these finding reveal how promoters naturally recruit TFIID to +1 nucleosomes and deliver TBP to the core promoter via activator-guided intra-TFIID hand-off of TBP.


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
