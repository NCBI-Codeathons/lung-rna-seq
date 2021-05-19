library(GEOquery)

args = commandArgs(trailingOnly=TRUE)

filePaths <- getGEOSuppFiles(args[1])
# TODO Turn these into their own nextflow process
GEO_ID_path <- paste0(args[2], args[1], "/")
tarF <- list.files(path = GEO_ID_path, pattern = "*.tar", full.names = TRUE)
untar(tarF, exdir = GEO_ID_path)
gzipF <- list.files(path = GEO_ID_path, pattern = "*.gz", full.names = TRUE)
ldply(.data = gzipF, .fun = gunzip)
