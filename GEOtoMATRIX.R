#Download GSE Data from NCBI to make a list with [1]matrix and [2]sample information
GEO_ID = "GSE132771"
a = GEOtoMATRIX(GEO_ID) 
#It's hard to run such big data in local

GEOtoMATRIX <- function(GEO_ID){
  filePaths = getGEOSuppFiles(GEO_ID)  
  GEO_ID_path = paste0("./",GEO_ID,"/")
  tarF <- list.files(path = GEO_ID_path, pattern = "*.tar", full.names = TRUE)  
  untar(tarF, exdir = GEO_ID_path)  
  gzipF <- list.files(path = GEO_ID_path, pattern = "*.gz", full.names = TRUE)   
  ldply(.data = gzipF, .fun = gunzip)  
  matrix_file = sort(list.files(path = GEO_ID_path, pattern = "\\.mtx$",full.names = TRUE) )
  gene_file = sort(list.files(path = GEO_ID_path, pattern = "*.genes.tsv$", full.names = TRUE) ) 
  barcodes_file = sort(list.files(path = GEO_ID_path, pattern = "*.barcodes.tsv$", full.names = TRUE))
  for (i in 1:length(matrix_file)){
    sc_maxtrix = readMM(file = matrix_file[i]) 
    sc_maxtrix = as.matrix(sc_maxtrix)
    sc_gene = read.table(file = gene_file[i],  
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
    rownames(sc_maxtrix) <- sc_gene[,2] 
    sc_barcodes <- read.table(file = barcodes_file[i],  
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
    colnames(sc_maxtrix) <- sc_barcodes[,1] 
    colnames(sc_maxtrix) <- paste(colnames(sc_maxtrix), barcodes_file[i], sep = "_") 
    sc_pdat <- data.frame("samples" = colnames(sc_maxtrix), "treatment" = barcodes_file[i])
    if (i == 1){
      matrix_store = sc_maxtrix
      pdat_store = sc_pdat
    }else{
      matrix_store = cbind(matrix_store,sc_maxtrix)
      pdat_store = rbind(pdat_store,sc_pdat)
    }
    
  }
  return(list(matrix_store, pdat_store))
}

