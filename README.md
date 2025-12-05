[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17829597.svg)](https://doi.org/10.5281/zenodo.17829597)

# SkeletalVis-Transcriptomics

## Introduction

**SkeletalVis-Transcriptomics** is a bioinformatics pipeline for reproducible analyses of microarray and RNA-seq data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a portable workflow tool to run tasks across multiple compute infrastructures. This pipeline uses a singularity container containing all the software needed to run the analysis, making installation simple and the results reproducible.

## Pipeline summary

The **SkeletalVis-Transcriptomics** pipeline takes a sample table and a parameter file defining the experiment as input. If not provided microarray data and fastq files are automatically downloaded using the provided accession numbers/sample identifiers.

### Features:
(**a**) Download of raw microarray data from GEO, fastq files either directly from ENA, via conversion of sra files from SRA<br/>
(**b**) Microarray with affyQCReport and RNA-seq read quality trimming with trimmomatic, QC reports with fastqc and multiQC<br/>
(**c**)	RNA-seq Quantification using [`kallisto`](https://pachterlab.github.io/kallisto/) and processing with tximport to produce a sample x gene expression table<br/>
(**d**) Differential expression analysis with limma, DESeq2 and Characteristic Direction<br/>
(**e**) Pathway and gene ontology enrichment analysis with goseq[`goseq`](https://bioconductor.org/packages/release/bioc/html/goseq.html)<br/>
(**f**) Active subnetwork identification with [`GIGA`](https://pubmed.ncbi.nlm.nih.gov/15272936/)<br/>
(**i**) Identify transcription factors potenitally driving differential expression with [`CHEA3`](https://pubmed.ncbi.nlm.nih.gov/31114921/)<br/>

Analyses are run in parallel and in result of error you can resume with the `-resume` parameter to re-run the pipeline starting from the previous fault.

## Quick Start

### Analyse an example dataset

Try the pipeline on an example dataset (all inputs will be automatically downloaded): -

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

3. Download the pipeline
   ```console
     nextflow clone CBFLivUni/SkeletalVis-Transcriptomics
   ```
4. [`Configure`](https://www.nextflow.io/docs/latest/config.html) the resource profile for your HPC or local computer. A template for slurm schedulers is provided as an example in `nextflow.config`

There is a utility function provided to help replace paths within the config text files:

    ```console
     bash scripts/install/replacePath.sh nextflow.config /mnt/hc-storage/groups/cbf/Nextflow/SkeletalVis-Transcriptomics `pwd -P`
    ```
5. Test on the example dataset:

    ```console
      nextflow run main.nf -profile slurm -params-file params/GSE152805.yaml -with-singularity library://jsoul/default/skeletalvis-transcriptomics:latest
    ```
### Analyse your own data

1. Define the sampleTable for RNA-seq data

Create a tab seperated table with unique Sample names, SRR accession numbers (if download is needed) and any additional metadata e.g

|Sample|File|Condition|
| ---|---|---|
|Control_1|SRRXXX	|Control|
|Control_2|	SRRXXX	|Control|
|Treated_1|	SRRXXX	|Treated|
|Treated_2|	SRRXXX	|Treated|

Note for microarray data the metadata is retrived directly from GEO and we instead just need to specify the columns of interest that define variables to compare.

2. Define the configuration

Most parameters are set to sensible defaults within the main nextflow script, with only a few parameters required to be altered with typical use. Note the use of Groovy, python and R booleans.

|Parameter|Description|Options|
| ---|---|---|
|accession|A unique identifier for the experiment to be analysed e.g the GEO accession of the data - used to name output data and download fastq files||
|species|The species the reads originate from - used to create the kallisto index	|Human, Mouse, Rat, Cow, Pig|
|single|Is the data single ended RNA-seq?	|true, false|
|batchCorrect|Should batch effect correction (sva) be used?	|TRUE, FALSE|
|skipTrimming|Should read trimming be skipped?|false (default), true|

Parameters should be defined within a yaml file. See `params/GSE152805.yaml` for an example.

The `accession` parameter defines the default search path for fastq.gz files (data/`accession`/fastqFiles/). Trimmed unpaired reads e.g "*_R0.fastq.gz" are skipped by default. If fastq files are not found locally the data will be downloaded using the provided `accession` number.

3. Run the pipeline with your own parameters

    ```console
     nextflow run soulj/SkeletalVis-Transcriptomics -profile slurm -params-file ownData.yaml -with-singularity library://jsoul/default/skeletalvis-transcriptomics
    ```

### Testing modules
Modules can be tested using the [`pytest-workflow`](https://pypi.org/project/pytest-workflow/) framework. Module test directories within the `tests` folder contain a nextflow script and a configuration yaml file defining the test for each module.

1. Install pytest-workflow

    ```console
	conda install pytest-workflow
    ```

2. Run the tests - e.g to test the pathway enrichment module

    ```console
	pytest --symlink --kwdof --tag pathwayEnrichment
    ```


