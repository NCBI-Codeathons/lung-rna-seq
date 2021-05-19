## Project Overview
Create an automated pipeline for single-cell RNA sequencing analysis that has UMAP visualization, automatic cell type identification (vi scSorter), differential gene expression analysis, and RNA velocity using the package velocyto. We also aim to link the databases Metascape and  genecards, and integrate the two datasets. The goal is to see both datasets individually and visualize similariteis and differences between them. For example, _what cell types are shared between these data? what genes are up- or down-regulated?_, etc. 

**Working Methodology-** use nf-core framework for development of this pipeline (general information below)

**Automated Analysis-** cell sorter and RNA velocity

## Introduction
RNA sequencing (RNA-seq) uses next-generation sequencing to examine the quantity of messenger RNA molecules in a biological sample and uses this transcriptomic information to extrapolate expression levels and changes in expression at the genomic level.[1] The ability to quantify genetic expression in different physiological or pathological states allows researchers to identify potential therapeutic targets and better understand the pathways underlying transitional processes from one state to another. Single-cell RNA sequencing (scRNA-seq), as the name suggests, is an RNA-seq approach that captures the transcriptome of individual cells, allowing researchers to identify cell subtypes and examine differences and similarities among these subtypes at the cellular level.[1]  Doing so provides increased resolution into cellular function within different physiological states. Both RNA- and single-cell RNA sequencing are computationally driven processes; as such, robust analysis pipelines with sound mathematical and biological are needed to drive data analysis and research.

With over a decade of use, there are many publicly available resources to help researchers store, analyze, visualize, and share scRNA-seq data, such as Seurat,[2–4]  Monocle,[5–9] and Scanpy.[10] Many of these resources require domain specific biological knowledge as well as coding experience. To cater to biologists without previous coding experience, resources like PIVOT[11] and CellRanger[12] have been developed. While the availability of such programs improves scRNA-seq data analysis, researchers often must use a combination of several different programs per dataset to complete analysis. This workflow also leaves room for pre-existing biases to shape analyses at critical steps. For example, cell type identification is usually performed manually, meaning annotation is based on pre-existing knowledge of marker genes. These groups are then used for downstream analyses, so the possibility of manual error can shift the entire analysis and may prevent researchers from identifying cell types that they do not already know. 

To address these issues, we have merged Seurat, Monocle, and velocyto into a single automated pipeline that is user-friendly. We are also appending the automatic cell identification package scSorter to hep eliminate cell-identification bias and compare against a library of marker genes, making cell identification more expansive than manual identification. 
 

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
To demonstrate our work's effectiveness, we will analyze the following ARDS dataset: GSE128033
(Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8025672/)

**Biological example-** compare normal (healthy patients) and damaged (COVID-19 afflicted patients) lung tissue.

**Data-** single cell RNA seq data from human bronchial epithelial cells infected with SARS-CoV-2. The data provided in this example are in the most common scRNA-seq file format (matrix/gene/barcode) and the file sizes are a reasonable size for quick analysis.

**Data collection-** differentiated human bronchial epithelial cells were infected with SARS-CoV-2 at 1, 2, and 3 days post-infection. Cells were then harvested and processed for single cell RNA-seq using the 10X Genomic platform. Generated libraries were sequenced on NovaSeq 6000 system using HiSeq 100 base pair reads and dual indexing.

## Credits

nf-core/teamrna was originally written by Edmund Miller. Edmund Miller, Kaitlyn Saunders, Yan Fang, and Alexa M. Salsbury contributed to the development of the pipeline and documentation throughout the NCBI North Texas Codeathon event (2021). 

## Citations
[1] 	Kukurba, K. R.; Montgomery, S. B. RNA Sequencing and Analysis. Cold Spring Harb. Protoc. 2015, 2015 (11), 951–969.

[2] 	Satija, R.; Farrell, J. A.; Gennert, D.; Schier, A. F.; Regev, A. Spatial Reconstruction of Single-Cell Gene Expression Data. Nat. Biotechnol. 2015, 33 (5), 495–502.

[3] 	Butler, A.; Hoffman, P.; Smibert, P.; Papalexi, E.; Satija, R. Integrating Single-Cell Transcriptomic Data across Different Conditions, Technologies, and Species. Nat. Biotechnol. 2018, 36 (5), 411–420.

[4] 	Hao, Y.; Hao, S.; Andersen-Nissen, E.; Mauck, W. M.; Zheng, S.; Butler, A.; Lee, M. J.; Wilk, A. J.; Darby, C.; Zagar, M.; et al. Integrated Analysis of Multimodal Single-Cell Data. bioRxiv. bioRxiv October 12, 2020.

[5] 	Monocle 3 https://cole-trapnell-lab.github.io/monocle3/docs/citations/ (accessed May 19, 2021).

[6] 	Trapnell, C.; Cacchiarelli, D.; Grimsby, J.; Pokharel, P.; Li, S.; Morse, M.; Lennon, N. J.; Livak, K. J.; Mikkelsen, T. S.; Rinn, J. L. The Dynamics and Regulators of Cell Fate Decisions Are Revealed by Pseudotemporal Ordering of Single Cells. Nat. Biotechnol. 2014, 32 (4), 381–386.

[7] 	Qiu, X.; Mao, Q.; Tang, Y.; Wang, L.; Chawla, R.; Pliner, H. A.; Trapnell, C. Reversed Graph Embedding Resolves Complex Single-Cell Trajectories. Nat. Methods 2017, 14 (10), 979–982.

[8] 	Qiu, X.; Hill, A.; Packer, J.; Lin, D.; Ma, Y. A.; Trapnell, C. Single-Cell MRNA Quantification and Differential Analysis with Census. Nat. Methods 2017, 14 (3), 309–315.

[9] 	Cao, J.; Spielmann, M.; Qiu, X.; Huang, X.; Ibrahim, D. M.; Hill, A. J.; Zhang, F.; Mundlos, S.; Christiansen, L.; Steemers, F. J.; et al. The Single-Cell Transcriptional Landscape of Mammalian Organogenesis. Nature 2019, 566 (7745), 496–502.

[10] 	Wolf, F. A.; Angerer, P.; Theis, F. J. SCANPY: Large-Scale Single-Cell Gene Expression Data Analysis. Genome Biol. 2018, 19 (1), 15.

[11] 	Zhu, Q.; Fisher, S. A.; Dueck, H.; Middleton, S.; Khaladkar, M.; Kim, J. PIVOT: Platform for Interactive Analysis and Visualization of Transcriptomics Data. BMC Bioinformatics 2018, 19 (6), 1–8.

[12] 	10x Genomics Cell Ranger 3.0.0. 2021.



CITE NF-CORE <!-- If you use  nf-core/teamrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->



> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
