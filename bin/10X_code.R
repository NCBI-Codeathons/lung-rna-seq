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
 
  library(BiocManager)
  library(GEOquery) 
 
 #dataset to use :D https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766 
 
#Put NCBI GEO# into the inside of the parentheses. (ex. GSE12345)  
filePaths = getGEOSuppFiles("GSE166766")  
  
#Copy paste the GEO12345 number between the "path = / /" sections   
tarF <- list.files(path = "./GSE166766/", pattern = "*.tar", full.names = TRUE)  
  
#This will extract that single raw data file and unzip the files within it.  
#Copy past the GEO12345 number between the "path = / /" sections  
  
untar(tarF, exdir = "./GSE166766/")  
gzipF <- list.files(path = "./GSE166766/", pattern = "*.gz", full.names = TRUE)  
ldply(.data = gzipF, .fun = gunzip)  
  
  
list.files(path = "./GSE166766/", pattern = "\\.mtx$",full.names = TRUE) 
list.files(path = "./GSE166766/", pattern = "*.genes.tsv$", full.names = TRUE)  
list.files(path = "./GSE166766/", pattern = "*.barcodes.tsv$", full.names = TRUE)  
  
  
#MATRIX CREATION (MAIN PART TO BE AUTOMATED essentially it is just just copy-paste the file paths from lines 94-96 into each "file =" =========================================================================================================
 
  library(plyr)
  library(dplyr) 
  library(Matrix) 
  library(Seurat) 
  library(cowplot) 
  library(SAVER) 
 
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

  library(plyr)
  library(dplyr) 
  library(Matrix) 
  library(Seurat) 
  library(cowplot) 
  library(SAVER) 
 
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

  library(plyr)
  library(dplyr) 
  library(Matrix) 
  library(Seurat) 
  library(cowplot) 
  library(SAVER)
 
mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
#print(mito.genes)
percent.mito <- Matrix::colSums(sobj@assays[["RNA"]][mito.genes, ])/Matrix::colSums(sobj@assays[["RNA"]])
sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito") 

VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.001)

par(mfrow = c(1, 2))
FeatureScatter(object = sobj, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 0.001)
FeatureScatter(object = sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.001)

#__________________________________________________________________________________________

ncountQup <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
ncountQlow <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nfeatQup <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nfeatQlow <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
mitoQup <- quantile(sobj@meta.data$percent.mito, 0.05)
mitoQlow <- quantile(sobj@meta.data$percent.mito, 0.95)

#combo
sobj <- subset(x = sobj, subset = nFeature_RNA > nfeatQlow & nFeature_RNA < nfeatQup)
sobj <- subset(x = sobj, subset = percent.mito > mitoQlow & percent.mito < mitoQup)
sobj <- subset(x = sobj, subset = nCount_RNA > ncountQlow & nCount_RNA < ncountQup)


sobj <- NormalizeData(object = sobj, normalization.method = "LogNormalize", scale.factor = 10000)
#make the features number the same number of variable genes you want to use
sobj <- FindVariableFeatures(object = sobj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)


#scSORTER stuff ==================================================================================================================

  library(scSorter)
  library(Seurat)
 
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

  library(plyr)
  library(dplyr) 
  library(Matrix) 
  library(Seurat) 
  library(cowplot) 
  library(SAVER) 
 
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

library(plyr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot)
library(metap)
library(multtest)
 
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

 
library(plyr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot)
library(metap)
library(multtest)
library(DESeq2)
library(msigdbr)
library(fgsea)
 
markers <- FindAllMarkers(object = sobj, min.pct = 0.25, thresh.use = 0.25) 
#markers <- markers[ markers$p_val_adj < 0.01, ] #can set p-value cutoff
top.markers <- do.call(rbind, lapply(split(markers, markers$cluster), head))
DoHeatmap(sobj, features = top.markers$gene, group.bar = TRUE)


#RNA Velocity Code ================================================================================================================

library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
 
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

 
 # monocle 3=========================================================================
BiocManager::install()
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('cole-trapnell-lab/monocle3',ref='develop')
library(monocle3)
library(dplyr)

# convert the Seurat object to formmat used in monocle 3
cds <- as.cell_data_set(sobj) # use seurat project sobj here
cds <- cluster_cells(cds = sobj.cds, reduction_method = "UMAP")
cds <- learn_graph(sobj.cds, use_partition = TRUE)
cds <- order_cells(sobj.cds, reduction_method = "UMAP", root_cells = hsc)
# Generate a cell_data_set from 10X output
cds <- load_mm_data(mat_path = matrix_file, 
                    feature_anno_path = gene_file , 
                    cell_anno_path =  barcodes_file)
 
#Pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
#Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)
rowData(cds)$gene_short_name <- rowData(cds)$V2
#Clusterring
cds = cluster_cells(cds)
## Learn a graph
cds <- learn_graph(cds)
## Order cells
cds <- order_cells(cds)
plot_cells(cds)

#Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(8, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    ordering_type="cluster_row_col",
                    max.size=3)
library(stringr)
colData(cds)$sample =  str_sub(rownames(colData(cds)),-1,-1)
colData(cds_subset)$sample =  str_sub(rownames(colData(cds_subset)),-1,-1)
#Subset cells
cds_subset <- choose_cells(cds)
cds_subset <- reduce_dimension(cds_subset)
#dentify genes that are differentially expressed in different subsets of
# cells from this partition:
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)
plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
cds_subset = cluster_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="cluster")

#Trajectory
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
#pseudotime
cds_subset <- order_cells(cds_subset)                  
plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

AFD_genes <- c("Ly6c2", "Mki67", "Pltp", "Arg1", "Adgre1")
AFD_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% AFD_genes, ]
plot_genes_in_pseudotime(AFD_lineage_cds, min_expr=0.5)


# find the genes that are differentially expressed on the different paths through the trajectory
ciliated_cds_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
a = subset(ciliated_cds_pr_test_res, q_value < 0.05)
b = a[order(a$p_value), ]

AFD_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% AFD_genes, ]
plot_genes_in_pseudotime(AFD_lineage_cds, min_expr=0.5)
