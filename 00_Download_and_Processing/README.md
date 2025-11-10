

This directory includes the shell and SLURM scripts for downloading and preprocessing the sequencing samples for analysis throughout the rest of this repo.

The following scripts should be executed in numerical order. The user should update the working directory (`$WRK`) filepath variable within every script before executing.

### `sampleIDs.txt`
The set of ids for samples from PEGR to analyze

### `MEMEsampleIDs.txt`
The set of ids for samples from PEGR to pull MEME motif analysis (PWMs) for.

### 0_download_data.sh
TODESCRIBE

```
bash 0_download_data.sh
```

### 1_align_data.sbatch
Script that mimics Galaxy core pipeline to perform initial alignment off of raw FASTQ data and remove duplicates.

```
# ^change the number of BAM files samples (SBATCH --array)
# To execute, type
sbatch 1_align_data.sbatch
##
# To check status, type
sbatch -u <myusername> -t
```

TODO: mimic core pipeline MEME motif calling

### 2_Merge_Replicates.sh
Merge BAM files for technical replicates (and IgG biological replicates) and rename BAM files to use standard naming system that includes experiment metadata (Cell line, IP target, antibody, assay-type, replicate, and genome build).

```
bash 2_Merge_Replicates.sh
```

### 3_normalize_samples.sbatch
Both TotalTag and NCIS normalization factors are calculated for each `data/BAM/SAMPLE.bam` and saved to `data/BAM/NormalizationFactors/` with the name `SAMPLE_NCISb.out`. For the NCIS, a blacklist reference and IgG control BAMs that are cell-line and assay-specific are input based on the standard BAM filename structure (parse `_` delimited tokens for assay and cell line info).

```
# ^change the number of BAM files samples (SBATCH --array)
# To execute, type
sbatch 3_normalize_samples.sbatch
##
# To check status, type
sbatch -u <myusername> -t
```
### 4_download_conservation-snp.sh
Download conservation and SNPs from USCS browser.
make sure proper scripts in bin
```
bash 4_download_conservation-snp.sh
```
