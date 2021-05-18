#PACKAGE INSTALLATIONS =====================================================

install.packages('BiocManager') 
BiocManager::install("GEOquery") 
install.packages("plyr") 
install.packages("dplyr")
install.packages("Matrix") 
BiocManager::install('multtest') 
install.packages('metap')
install.packages("Seurat") 
install.packages('ggplot2') 
install.packages("cowplot") 
install.packages("devtools") 
devtools::install_github("mohuangx/SAVER") 
BiocManager::install("DESeq2")
install.packages("msigdbr")
BiocManager::install("fgsea")


#note, and I quote, "Monocle3 is hell to install". But here are their instructions:
https://cole-trapnell-lab.github.io/monocle3/docs/installation/ 
  
  
  #installing velocyto.R (package for RNA velocity analysis) instructions
  https://github.com/velocyto-team/velocyto.R

BiocManager::install("pcaMethods")
library(devtools)
library(pcaMethods)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)


#install LoomR (file type input into velocyto.R for RNA velocity analysis)
(library(devtools)
  devtools::install_github(repo = "hhoeflin/hdf5r")
  install.packages("stringi")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  
  
  #Install SeuratWrappers if you need it (will need for RNA velocity)
  install.packages("remotes")
  install.packages("vctrs")
  install.packages("lifecycle")
  rtools <- "C:\\Rtools\\bin"
  gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
  path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
  new_path <- c(rtools, gcc, path)
  new_path <- new_path[!duplicated(tolower(new_path))]
  Sys.setenv(PATH = paste(new_path, collapse = ";"))
  install.packages("ellipse")       # Install the dependent package ellipse
  install.packages('/tmp/RtmpCAMbyw/filedabd3e2e7efe/SeuratWrappers_0.3.0.tar.gz',  type="source", repos=NULL)
  
  
#LIBRARIES ===============================================================================================================
  
  library(BiocManager)
  library(GEOquery) 
  library(plyr)
  library(dplyr) 
  library(Matrix) 
  library(devtools)
  library(Seurat) 
  library(ggplot2) 
  library(cowplot) 
  library(SAVER) 
  library(metap)
  library(multtest)
  library(DESeq2)
  library(msigdbr)
  library(fgsea)
  library(monocle3)
  library(velocyto.R)
  library(loomR)
  
  
  
#LOADING IN DATA =========================================================================================================
  
#Put NCBI GEO# into the inside of the parentheses. (ex. GSE12345)  
filePaths = getGEOSuppFiles("GSE132771")  
  
#Copy paste the GEO12345 number between the "path = / /" sections   
tarF <- list.files(path = "./GSE132771/", pattern = "*.tar", full.names = TRUE)  
  
#This will extract that single raw data file and unzip the files within it.  
#Copy past the GEO12345 number between the "path = / /" sections  
  
untar(tarF, exdir = "./GSE132771/")  
gzipF <- list.files(path = "./GSE132771/", pattern = "*.gz", full.names = TRUE)  
ldply(.data = gzipF, .fun = gunzip)  
  
  
list.files(path = "./GSE132771/", pattern = "\\.mtx$",full.names = TRUE) 
list.files(path = "./GSE132771/", pattern = "*.genes.tsv$", full.names = TRUE)  
list.files(path = "./GSE132771/", pattern = "*.barcodes.tsv$", full.names = TRUE)  
  
  
#MATRIX CREATION (MAIN PART TO BE AUTOMATED essentially it is just just copy-paste the file paths from lines 94-96 into each "file =" =========================================================================================================
  
#control group______________________________________________________________________________________  
control_matrix <- readMM(file = './GSE119352/GSM3371684_Control_matrix.mtx') #copy-paste .mtx file path here
control_matrix <- as.matrix(control_matrix)
  
control_genes <- read.table(file = './GSE119352/GSM3371684_Control_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE) #copy-paste .genes.tsv file path here
rownames(control_matrix) <- control_genes[,2]
  
control_barcodes <- read.table(file = './GSE119352/GSM3371684_Control_barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE) #copy-paste .barcodes.tsv file path here
colnames(control_matrix) <- control_barcodes[,1]
  
colnames(control_matrix) <- paste(colnames(control_matrix), "control", sep = "_")
control_pdat <- data.frame("samples" = colnames(control_matrix), "treatment" = "control")

#antiPD1 group ________________________________________________________________________________________
aPD1_matrix <- readMM(file = './GSE119352/GSM3371685_aPD1_matrix.mtx')
aPD1_matrix <- as.matrix(aPD1_matrix)
aPD1_matrix[1:5,1:5]
#dim(aPD1_matrix)

aPD1_genes <- read.table(file = './GSE119352/GSM3371685_aPD1_genes.tsv', 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(aPD1_genes)
rownames(aPD1_matrix) <- aPD1_genes[,2]
aPD1_matrix[1:5,1:5]

aPD1_barcodes <- read.table(file = './GSE119352/GSM3371685_aPD1_barcodes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(aPD1_barcodes)
colnames(aPD1_matrix) <- aPD1_barcodes[,1]
aPD1_matrix[1:5,1:5]
colnames(aPD1_matrix) <- paste(colnames(aPD1_matrix), "aPD1", sep = "_")

aPD1_pdat <- data.frame("samples" = colnames(aPD1_matrix), "treatment" = "aPD1")


# MERGING MATRICES FOR SEURAT OBJECT CREATION ===========================================================================================

joined <- cbind(control_matrix,aPD1_matrix)
#dim(joined)
pdat <- rbind(control_pdat, aPD1_pdat)

rownames(pdat) <- pdat$samples
fdat <- toupper(as.matrix(control_genes))

rownames(fdat) <- fdat[,2]
fdat <- data.frame(fdat)
common_colnames <- c("ensembl_id", "gene_short_name")
colnames(fdat) <- common_colnames

rownames(joined) <- rownames(fdat)

#Seurat Object Creation______________________________________________________
sobj<- CreateSeuratObject(counts = joined)
sobj<-AddMetaData(sobj,metadata=pdat)
#head(sobj@meta.data, 5)

sobj[["RNA"]]@meta.features<-fdat
#head(sobj[["RNA"]]@meta.features)
slotNames(sobj[["RNA"]])


#PRE-PROCESSING =============================================================================================================

mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
#print(mito.genes)
percent.mito <- Matrix::colSums(sobj@assays[["RNA"]][mito.genes, ])/Matrix::colSums(sobj@assays[["RNA"]])
sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito") 

VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.001)

par(mfrow = c(1, 2))
FeatureScatter(object = sobj, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 0.001)
FeatureScatter(object = sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.001)

#__________________________________________________________________________________________

quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
quantile(sobj@meta.data$percent.mito, 0.05)
quantile(sobj@meta.data$percent.mito, 0.95)

#combo
sobj <- subset(x = sobj, subset = nFeature_RNA > 1085 & nFeature_RNA < 4456)
sobj <- subset(x = sobj, subset = percent.mito > 0.003 & percent.mito < 0.015)
sobj <- subset(x = sobj, subset = nCount_RNA > 2432 & nCount_RNA < 21076)


sobj <- NormalizeData(object = sobj, normalization.method = "LogNormalize", scale.factor = 10000)
#make the features number the same number of variable genes you want to use
sobj <- FindVariableFeatures(object = sobj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)


#scSORTER stuff ==================================================================================================================

# Taken from: https://cran.r-project.org/web//packages/scSorter/vignettes/scSorter.html
#expr = GetAssayData(expr_obj)
#topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
#topgenes = topgenes[topgene_filter]

#picked_genes = unique(c(anno$Marker, topgenes))
#expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, anno)


#if you don't want to wait all day, remove the "feature = rownames(sobj)" parameter
sobj <- ScaleData(object = sobj, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(sobj), block.size = 2000)


#PCA AND UMAP ===============================================================================================================================

sobj <- RunPCA(sobj, npcs = 100, ndims.print = 1:10, nfeatures.print = 5)
ElbowPlot(sobj, ndims = 100) #based on where base of "elbow" is, chooses number of principal components to use in dimension reduction

#alternatively, could use this to automate selection of number of principal components but may want to double check its actual effectiveness:
pct <- sobj[["pca"]]@stdev / sum(sobj[["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
# Determine which PC exhibits cumulative percent greater than 
#90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs #output of "pcs" is number need to put here: "dims = 1:pcs"

sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:25, nn.eps = 0.5) 
# for more clusters, increase resolution, and for fewer, decrease it
sobj <- FindClusters(sobj, resolution = 0.4, n.start = 10)
#changing seed number rotates clusters
sobj <- RunUMAP(object = sobj, reduction = "pca", dims = 1:9,min.dist = 0.5, seed.use = 123)

p1 <- DimPlot(sobj, reduction = "umap", pt.size = 0.1,label=TRUE) + ggtitle(label = "Plot")
p1
DimPlot(sobj, reduction = "umap", split.by = "treatment")


#CLUSTER IDENTIFICATION (MANUAL) =================================================================================================

#knowing which cluster is which_______________________________________________________________________________________
preknownmarkerlist <- c("HAVCR2", "GZMB", "VCAM1", "PRF1",  
             "KLRC1", "CCL4", "LAG3", "CCL3")
sobj <- AddModuleScore(object = sobj, features = cd8trmlist, name = "CD8TRM_List")  
FeaturePlot(object = sobj, features = "CD8TRM_List1") #MUST add 1 to the end of the name as that's how Seurat works

#Renaming clusters______________________________________________________________________________________________________
sobj<- RenameIdents(sobj, `3` = "CD8+ T(RM)", `8` = "CD8+ T(RM) mitotic", `6` = "CD8+ T(EM)",
                    `2` = "CD8+ gamma-delta", `4` = "CD4+ CD103+", 
                    `7` = "CD4+ CXCL13+",`1` = "CD4+ FOXP3+", `0`="CD4+ IL-7R+",
                    `5`= "CD4+ RGCC+",`9`= "CD3- Mono")
DimPlot(sobj, reduction = "umap", pt.size = 0.1,label=TRUE) + ggtitle(label = "Named Plot")
DimPlot(sobj, split.by = "treatment")


#DIFFERENTIAL GENE EXPRESSION ANALYSIS ===========================================================================================

markers <- FindAllMarkers(object = sobj, min.pct = 0.25, thresh.use = 0.25) 
#markers <- markers[ markers$p_val_adj < 0.01, ] #can set p-value cutoff
top.markers <- do.call(rbind, lapply(split(markers, markers$cluster), head))
DoHeatmap(sobj, features = top.markers$gene, group.bar = TRUE)


#RNA Velocity Code ================================================================================================================

#See "Converting to/from loom" : https://satijalab.org/seurat/articles/conversion_vignette.html 
#https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html 

# If you don't have velocyto's example mouse bone marrow dataset, download with the CURL command
# curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile
# = '~/Downloads/SCG71.loom')
ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)
