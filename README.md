# Bluk RNA-seq pipeline: from fastq processing to differential gene expression analysis

* [What new features does this pipeline introduce?](#what-new-features-does-this-pipeline-introduce)
* [Getting data from GEO](#getting-data-from-geo)
* Trimming adapters with Cutadapt
* Mapping reads with STAR
  * Building a genome index
  * Align the FASTQ with the genome index
* Differential Expression Analysis

### What new features does this pipeline introduce?

RNA-sequencing provides an averaged transcriptomic profile, which can be analysed to reveal the differences in gene expression between groups of samples. It can be used for:

* **Differential Gene Expression (DGE) Analysis**: compare gene expression between different conditions, such as disease vs. healthy tissue, treated vs. untreated.
* **Functional Genomics and Pathway Analysis**: identify gene regulatory networks and co-expression patterns, as well as biological pathways and processes.
* **Biomarker Discovery**: find potential diagnostic, prognostic, or therapeutic biomarkers in diseases like cancer.
* **Cell Type Composition Estimation**: use deconvolution methods to estimate the proportion of cell populations in heterogeneous tissues.

However, [scientific evidence](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4) reveals that the mainstream methods (Deseq2 and EdgeR) of differential gene expression analysis suffers from an exaggerated number of false positives detected for differentially expressed genes.

Accounting this problem, this pipeline uses a costumized Wilcoxon rank-sum pipeline to differential gene expression analysis in bulk data, comparing both groups of samples and paired patient samples (Adjacent vs. Tumor for each patient).

### Getting data from GEO

With the bash script "1-download_sra.sh", we can call It in command line using a .txt file as argument, containing one SRR code for each line, which will be downloaded and converted to FASTQ and concatenated.

### Trimming adapters with Cutadapt

First of all, we can check the qulity of the FASTQ files with fastqc. It can be run with:

`fastqc --noextract --nogroup -o fastqc *.fastq.gz`

Then, we can run cutadapt to trim adapter sequences and increase the quality of our data, with this template:

`cutadapt -a ADAPTER_SEQUENCE -A ADAPTER_SEQUENCE_R2 -o trimmed_R1.fastq -p trimmed_R2.fastq R1.fastq R2.fastq`

For sanger data, for example, we use the following adapter sequences.:

`cutadapt -a TGTAAAACGACGGCCAGT -A CAGGAAACAGCTATGAC -o trimmed_R1.fastq -p trimmed_R2.fastq R1.fastq R2.fastq`
