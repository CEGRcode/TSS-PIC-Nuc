# Genome-wide rotational and translational setting of transcription factors with nucleosomes

### Jordan E. Krebs<sup>&#x2020; 1,2</sup>, Haining Chen<sup>&#x2020; 2</sup>, Olivia W. Lang<sup>2</sup>, William K. M. Lai<sup>2</sup>, B. Franklin Pugh<sup>2</sup>

&#x2020; Co-first author

<sup>1</sup>MD/PhD Medical Scientist Training Program, Penn State College of Medicine, Hershey, PA, USA.

<sup>2</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

### Correspondence:fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [GSE266547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266547)

## Abstract
How transcription factors (TFs) access DNA that is packed into chromatin has been challenging to decipher on a genomic scale due to inherent resolution limits of assays. Accessibility differs depending on whether a TF binding site (TFBS) is translationally and/or rotationally phased on or next to a nucleosome. Rotational phasing occurs where a DNA sequence in the DNA helix has a predominant rotational direction on the nucleosome surface, thereby consistently facing outward (accessible) or inward (inaccessible). While rotational phasing is DNA-encoded by dinucleotide periodicities at yeast TFBSs, such encoding has not been found in humans. Here we develop a genome-wide Benzonase nuclease-based assay to measure translational and rotational phasing, and a second assay to measure such phasing on the same DNA molecule to which a TF is bound. The latter uniquely allows phasing to be measured with transient TFs. We show that many types of human TFBSs have distinct translational and/or rotational phasing depending on whether they are bound by specific TFs. For example, unbound CTCF sites are nucleosomal with local but not global translational phasing. They possess DNA-encoded rotational phasing, with the predominant phase being in an accessible orientation. CTCF-bound sites have adjacent nucleosomes possessing a global translational phase and distinct DNA-encoded rotational phases. Similar themes recur for other TFs including NFIA and pioneer factor FoxA, but their activities are highly transient. Our assays further allow subnucleosomal structures to be examined along with their relationship to a transiting RNA polymerase. Together, these findings reveal an intimate relationship between nucleosome phasing and transcription factors.


## Directions
To recreate the figures for this manuscript, please execute the scripts in each directory in numerical order. Each directory's README includes more specific details on execution. To be more explicit, run the scripts in each directory in the following order: `00_Download_and_Preprocessing`, `01_Run_GenoPipe`, `02_Call_Nucleosomes`, `03_Call_JASPAR`, `04_Call_Motifs`, `X_Bulk_Processing`, and then finally `Z_Figures`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n bx -c bioconda -c conda-forge bedtools bowtie2 bwa cutadapt meme opencv pandas samtools scipy sra-tools wget pybigwig
```

For genetrack-executing script, a python2 environment needed to be created. The create command for that env is as follows:

```
conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy
```


## Table of Contents

### 00_Download_and_Preprocessing
Perform the preprocessing steps including alignment of raw sequencing data from both novel and previously published data.

### 01_Run_GenoPipe
Perform quality control for genetic background on these data by running GenoPipe on the aligned BAMs.

### 02_Call_Nucleosomes
Call nucleosome positions and identify TSS and +1 nucleosome reference points with different sorts.

### 03_Call_JASPAR
Call JASPAR motifs and subset to "bound" sites using ENCODE peak data.

### 04_Call_Motifs
Build de novo sequence-specific transcription factor (ssTF) motif reference points using Benzonase ChIP-exo data.

### X_Bulk_Processing
With the BAM and BED files built from the scripts in the above directories, perform bulk read pileups for heatmaps and composites.

### Z_Figures
Copy/organize results from bulk processing into figure-specific directories corresponding to subfigures in the manuscript. Also includes custom/one-off scripts for analysis that didn't need bulk-style execution.

### AI_files
all figures in paper

### data
Store large files to be globally accessed by the scripts in each directory

### bin
Generalized scripts and executables for global access by each of the numbered directories.
