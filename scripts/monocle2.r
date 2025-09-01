# --------- Prog Info ---------
# R version 3.6.3
# Author: Songming Liu
# E-Maile: smliu-bio@qq.com
# Date: 2021/04/01 ~ 2021/04/02
# Purpose: Trajectory Inference by using monocle2
# Attention: version 1.0, use given genes or high dispersion genes for trajectory inference

# Global Functions --------------------------------------------------------
checkPackage <- function(...) {
# Check whether package exists or not
    if (!suppressWarnings(require(...))) {
        stop(paste("Package:",...,"not exists, Please check!"))
    }
}
createDir <- function(dir){
# Check and Create dirs
    if (!dir.exists(dir)) {
        dir.create(path = dir,recursive = T)
    }
}

# Parse cmd and library packages ------------------------------------------
checkPackage("optparse")
option_ls <- list(
    make_option(opt_str = c("-i","--infile"),action = "store",dest = "infile",
                help = "UMI matrix, rds file of Seurat object or 10x-mtx dir."),
    make_option(opt_str = c("-o","--outdir"),action = "store",dest = "outdir",
                help = "directory to store results."),
    make_option(opt_str = "--intype",action = "store",dest = "intype",default = "seurat",
                help = "infile type [seurat, matrix, 10x] (default seurat)."),
    make_option(opt_str = "--gfile",action = "store",dest = "gfile",
                help = "gene file (one column, no header), use these genes to infer trajectory.")
)
obj_parse <- OptionParser(option_list = option_ls) # instance
arg_ls <- parse_args(object = obj_parse) # parse cmd
## check args
if (is.null(arg_ls$infile) || is.null(arg_ls$outdir)){
    stop("-i/-o are required and expected one argument, please check!")
}
if (!file.exists(arg_ls$infile)) {
    stop("infile not exist, please check!")
}
arg_ls$intype <- tolower(arg_ls$intype)
if (!(arg_ls$intype %in% c("seurat","matrix", '10x'))) {
    stop("--intype should be one of 'seurat', 'matrix'")
}
if ((!is.null(arg_ls$gfile)) && (!file.exists(arg_ls$gfile))) {
    stop("gene file (--gfile) not exist, please check!")
}
## loading required packages
checkPackage("Matrix")
checkPackage("monocle")
checkPackage("ggplot2")
checkPackage("tidyverse")
## print args
writeLines("--------------------- Parameters List ---------------------")
for (para in names(arg_ls)) {
    writeLines(paste0(para,": ",arg_ls[[para]]))
}
## create outdir
createDir(dir = arg_ls$outdir)
setwd(arg_ls$outdir)

# 1. Create CellDataSet Object --------------------------------------------
writeLines("\n--------------------- Step1. Create CellDataSet Object ---------------------")
if (arg_ls$intype == 'matrix') {
    writeLines("Creating (from UMI counts file) ...")
    raw_mat <- read_csv(arg_ls$infile) %>%
        column_to_rownames("...1")
    raw_mat <- as(as.matrix(raw_mat),"sparseMatrix")
    pd <- new("AnnotatedDataFrame",data = data.frame(orig.ident = rep(1,ncol(raw_mat)),
                                                     row.names = colnames(raw_mat)))
    fd <- new("AnnotatedDataFrame",data = data.frame(gene_short_name = rownames(raw_mat),
                                                     row.names = row.names(raw_mat)))
} else if (arg_ls$intype == '10x') {
    writeLines("Creating (from 10x format) ...")
    raw_mat <- Seurat::Read10X(arg_ls$infile)
    pd <- new("AnnotatedDataFrame",data = data.frame(orig.ident = rep(1,ncol(raw_mat)),
                                                     row.names = colnames(raw_mat)))
    fd <- new("AnnotatedDataFrame",data = data.frame(gene_short_name = rownames(raw_mat),
                                                     row.names = row.names(raw_mat)))
} else {
    writeLines("Creating (from Seurat counts slot) ...")
    scrna <- readRDS(arg_ls$infile)
    raw_mat <- scrna@assays$RNA@counts
    pd <- new("AnnotatedDataFrame",data = scrna@meta.data)
    fd <- new("AnnotatedDataFrame",data = data.frame(gene_short_name = rownames(scrna),
                                                     row.names = row.names(scrna)))
    rm(scrna)
    gc()  ## garbage collection
}
cds <- newCellDataSet(cellData = raw_mat,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())

# 2. Estimate size factors and dispersions --------------------------------
writeLines("\n--------------------- Step2. Estimate size factors and dispersions ---------------------")
cds <- estimateSizeFactors(cds) # useful for normalizing across cells
cds <- estimateDispersions(cds) # useful for DEA

# 3. Construct Single Cell Trajectories -----------------------------------
writeLines("\n--------------------- Step3. Construct Single Cell Trajectories ---------------------")
## select genes to order cells
if (is.null(arg_ls$gfile)) {
    writeLines("High Dispersion Genes Searching ...")
    disp_table <- dispersionTable(cds)
    ordGene_v <- subset(disp_table, mean_expression > 0.1 & dispersion_empirical >= dispersion_fit)$gene_id
} else {
    writeLines("Given Genes Reading ...")
    ordGene_v <- read.table(file = arg_ls$gfile,header = F,sep = "\t",stringsAsFactors = F)$V1
}
cds <- setOrderingFilter(cds,ordGene_v)
png("vlz01_orderingGenes.png",width = 6,height = 5,units = 'in',res = 300)
plot_ordering_genes(cds)
dev.off()
## reduce data dimensionality
writeLines("Dimensionality Reducing ...")
cds <- reduceDimension(cds,max_components = 2,reduction_method = "DDRTree")
## cell ordering
writeLines("Cell Ordering ...")
cds <- orderCells(cds)
## visualization
p <- plot_cell_trajectory(cds,show_cell_names = F,color_by = "Pseudotime") + 
    scale_color_viridis_c()
ggsave(filename = "vlz02_trajectory_Pseudotime.png",width = 6,height = 5,plot = p)
## save
writeLines("Environment Saving ...")
saveRDS(object = cds,file = "monocle.rds")
save.image(file = "monocle.RData")
