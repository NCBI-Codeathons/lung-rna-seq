// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GEO_DOWNLOAD {
    tag "$id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    conda (params.enable_conda ? "conda-forge::r-base=3.4.1 bioconda::bioconductor-geoquery==2.46.3"  : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-geoquery:2.58.0--r40hdfd78af_1"
    } else {
        container "quay.io/biocontainers/bioconductor-geoquery:2.58.0--r40hdfd78af_1"
    }

    input:
    val id

    output:
    path "*.mtx"            , optional:true    , emit: matrix
    path "*.genes.tsv"      , optional:true    , emit: genes
    path "*.barcodes.tsv"   , optional:true    , emit: barcodes
    // path "*.RData"       , optional:true    , emit: rdata
    // path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    echo $id > id.txt
    GEO_download.R ${id} ./
    """

}