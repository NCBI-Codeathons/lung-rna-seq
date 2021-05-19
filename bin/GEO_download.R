library(GEOquery)
library(optparse)

option_list <- list(
    make_option(c("-d", "--geo_id" ), type="character", default='GSE166766' , metavar="string" , help="GEO Id"            )
    make_option(c("-o", "--outdir" ), type="character", default='./'        , metavar="path"   , help="Output directory." ),
)

filePaths <- getGEOSuppFiles(opt$geo_id)
# TODO Turn these into their own nextflow process
GEO_ID_path <- paste0("./", opt$geo_id, "/")
tarF <- list.files(path = GEO_ID_path, pattern = "*.tar", full.names = TRUE)
untar(tarF, exdir = GEO_ID_path)
gzipF <- list.files(path = GEO_ID_path, pattern = "*.gz", full.names = TRUE)
ldply(.data = gzipF, .fun = gunzip)
