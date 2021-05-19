## Project Overview
Create an automated pipeline for single-cell RNA sequencing analysis that has UMAP visualization, automatic cell type identification (vi scSorter), differential gene expression analysis, and RNA velocity using the package velocyto. We also aim to link the databases Metascape and  genecards, and integrate the two datasets. The goal is to see both datasets individually and visualize similariteis and differences between them. For example, _what cell types are shared between these data? what genes are up- or down-regulated?_, etc. 

**Working Methodology-** use nf-core framework for development of this pipeline (general information below)

**Automated Analysis-** cell sorter and RNA velocity

## Introduction
RNA sequencing (RNA-seq) uses next-generation sequencing to examine the quantity of sequences in a biological sample. RNA-seq analyzes the transcriptome of gene expression patterns encoded within RNA during different physiological or pathological states. As such, RNA-seq can provide important information on fundamental biological processes as well as diseased states. Single-cell RNA sequencing (scRNA-seq), as the name suggests, is an RNA-seq approach that provides information on individual cells, allowing researchers to examine cell-to-cell differences an identify cell subtypes. Doing so provides more specific insight on cellular function within different physiological states. To ascertain such information from RNA-seq experiments, biologist must perform data analysis. 

## Methods

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/teamrna** is a bioinformatics best-practice analysis pipeline for Single Cell RNA-seq analysis pipeline.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/teamrna/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
Pepline dependency:
Seurat,
Monocle3,
scSORTER,
RNA Velocity


1. Download data from NCBI GEO dataset
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
4. Trajectory analysis uses monocle 3.
5. Reduce dimensionality, clustering and visualize the cells.
6. Find marker genes expressed by each cluster.
7. Trajectory assignment.
7. Pseudotime assign and plot gene along the pseudotime.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/teamrna -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. It is also highly recommended to use the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) settings to store the images in a central location for future pipeline runs.
    * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/teamrna -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

See [usage docs](https://nf-co.re/teamrna/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/teamrna pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/teamrna/usage) and [output](https://nf-co.re/teamrna/output).

## Example
To demonstrate our work's effectiveness, we will analyze the following ARDS dataset: GSE166766
(Reference: https://pubmed.ncbi.nlm.nih.gov/33730024/)

**Biological example-** compare normal (healthy patients) and damaged (COVID-19 afflicted patients) lung tissue.

**Data-** single cell RNA seq data from human bronchial epithelial cells infected with SARS-CoV-2. The data provided in this example are in the most common scRNA-seq file format (matrix/gene/barcode) and the file sizes are a reasonable size for quick analysis.

**Data collection-** differentiated human bronchial epithelial cells were infected with SARS-CoV-2 at 1, 2, and 3 days post-infection. Cells were then harvested and processed for single cell RNA-seq using the 10X Genomic platform. Generated libraries were sequenced on NovaSeq 6000 system using HiSeq 100 base pair reads and dual indexing.

## Credits

nf-core/teamrna was originally written by Edmund Miller. Edmund Miller, Kaitlyn Saunders, Yan Fang, and Alexa M. Salsbury contributed to the development of the pipeline and documentation throughout the NCBI North Texas Codeathon event (2021). 

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#teamrna` channel](https://nfcore.slack.com/channels/teamrna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations (NF-CORE) 

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/teamrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

## Other works of interest to our project
https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7406263/
https://www.biorxiv.org/content/10.1101/416719v2


> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
