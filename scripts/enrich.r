# -------- Prog Info -------
# R version 4.1.3 (2022-03-10)
# Author: Songming Liu
# Date: 2022-10-14
# Desc: GO/KEGG Enrichment of DEGs using clusterProfiler & msigdbr
# Major packages: msigdbr (7.5.1), clusterProfiler (4.2.0)

# Global Functions --------------------------------------------------------
CheckPackage <- function(...) {
  # Check whether package exists or not
  if (!suppressWarnings(require(...))) {
    stop(paste("Package:", ..., "not exists, Please check!"))
  }
}

CreateDir <- function(dir) {
  # Check and Create dirs
  if (!dir.exists(dir)) {
    dir.create(path = dir, recursive = T)
  }
}

VlzDot <- function(data, outfile, title, 
                   showNum = 10,
                   col = 'pvalue',
                   maxLen = 50) {
    # Visualize GO/KEGG results
    # Args:
    #     data: Enrichment results for visualization
    #     outfile: output file
    #     title: plot title
    #     showNum: number of pathway with min pvalue to show in graph
    #     col: column to color by
    #     maxLen: max length of pathway names
    data$GeneRatio <- sapply(as.character(data$GeneRatio),
                             function(x) {
                               return(eval(parse(text = x)))
                             })
    pathway_v <- data %>%
      top_n(showNum, wt = -log10(pvalue)) %>%
      top_n(showNum, wt = GeneRatio) %>%
      pull(ID)
    df <-
      as.data.frame(data[data$ID %in% pathway_v, ], stringsAsFactors = F)
    # select data using for plotting
    if (col == 'p.adjust') {
      df[['color']] <- -log10(df$p.adjust)
      colorname <- "-log10(p.adjust)"
    } else {
      df[['color']] <- -log10(df$pvalue)
      colorname <- "-log10(pvalue)"
    }
    # max length of y lab
    string_cut <- function(x, len) {
      if (nchar(x) > len) {
        return(paste(str_sub(x, 1, len), "..."))
      } else {
        return(x)
      }
    }
    # plot
    df$Description <-
      sapply(df$Description, string_cut, maxLen) # max length
    p <- ggplot(data = df, aes(GeneRatio, Description)) +
      geom_point(aes(size = Count, color = color)) +
      labs(
        title = title,
        x = 'Gene Ratio',
        y = '',
        color = colorname
      ) +
      scale_color_gradient(low = "blue", high = "red") +
      scale_size(range = c(2, 7)) +
      guides(size = guide_legend(
        reverse = T,
        nrow = 5,
        order = 1
      ),
      color = guide_colorbar(order = 2)) +
      theme_bw()
    ggsave(
      filename = outfile,
      plot = p,
      width = 8,
      height = 6
    )
  }

EnrichFlow <- function(genelist, txtFile, graphFile, graphTi, ...) {
  # Enrich for given genelist, then show results using dotplot
  # Args:
  #     genelist: given gene list
  #     txtFile: file name of txt results
  #     graphFile: file name of graph
  #     graphTi: graph title
  #     ...: other parameters for 'enricher' function
  res <- enricher(gene = genelist, ...)
  if (is.null(res)) {
    writeLines("no results!")
    return(NULL)
  }
  write.table(
    x = res@result,
    file = txtFile,
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  VlzDot(data = res@result,
         outfile = graphFile,
         title = graphTi)
}

# Parse cmd and library packages ---------------------------------------------------------------
CheckPackage("optparse")
option_ls <- list(
  make_option(
    opt_str = c("-i", "--infile"),
    action = "store",
    dest = "infile",
    help = "input file contains gene information."
  ),
  make_option(
    opt_str = c("-o", "--outdir"),
    action = "store",
    dest = "outdir",
    help = "directory to store results."
  ),
  make_option(
    opt_str = "--geneCol",
    action = "store",
    dest = "geneCol",
    default = "gene",
    help = "column names of genes (default 'gene')."
  ),
  make_option(
    opt_str = "--typeCol",
    action = "store",
    dest = "typeCol",
    default = "cluster",
    help = "column names of gene types (default 'cluster')."
  ),
  make_option(
    opt_str = "--geneType",
    action = "store",
    dest = "geneType",
    default = "SYMBOL",
    help = "gene name type of input file, can be ENTREZID or SYMBOL (default SYMBOL)."
  ),
  make_option(
    opt_str = "--organism",
    action = "store",
    dest = "org",
    default = "Homo-sapiens",
    help = "organism (default 'Homo-sapiens')."
  )
)
obj_parse <- OptionParser(option_list = option_ls) # instance
arg_ls <- parse_args(object = obj_parse) # parse cmd

## check args
if (is.null(arg_ls$infile) || is.null(arg_ls$outdir)) {
  stop("-i/-o are required and expected one argument, please check!")
}
if (!file.exists(arg_ls$infile)) {
  stop("infile not exist, please check!")
}
if (!arg_ls$geneType %in% c("SYMBOL", "ENTREZID")) {
  stop("--geneType should be one of [SYMBOL, ENTREZID]")
}
arg_ls$org <- gsub(pattern = "-", replacement = " ", arg_ls$org)
## loading required packages
CheckPackage("clusterProfiler")
CheckPackage("msigdbr")
CheckPackage("dplyr")
CheckPackage("stringr")
CheckPackage("ggplot2")
## print args
writeLines("--------------------- Parameters List ---------------------")
for (para in names(arg_ls)) {
  writeLines(paste0(para, ": ", arg_ls[[para]]))
}
## create outdir
CreateDir(dir = arg_ls$outdir)
setwd(arg_ls$outdir)

# 1. Preparation for Enrichment ------------------------------------------------
writeLines("\n--------------------- Step1 Preparation ---------------------")
## get gene list
writeLines("Infile Reading ...")
infile_df <- read.table(
  file = arg_ls$infile,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = "#"
)
gene_ls <- split(x = infile_df[[arg_ls$geneCol]],
                 f = infile_df[[arg_ls$typeCol]])

## prepare for GO
writeLines("GO Terms Auquiring ....")
go <- msigdbr(species = arg_ls$org, category = "C5")
if (arg_ls$geneType == "SYMBOL") {
  go <- go[, c("gs_id", "gs_name", "gs_subcat", "gene_symbol")]
} else {
  go <- go[, c("gs_id", "gs_name", "gs_subcat", "entrez_gene")]
}
colnames(go) <- c("id", "name", "subcat", "gene")

## prepare for KEGG
writeLines("KEGG Pathways Auquiring ....")
kegg <- msigdbr(species = arg_ls$org,
                category = "C2",
                subcategory = "CP:KEGG")
if (arg_ls$geneType == "SYMBOL") {
  kegg <- kegg[, c("gs_id", "gs_name", "gene_symbol")]
} else {
  kegg <- kegg[, c("gs_id", "gs_name", "entrez_gene")]
}
colnames(kegg) <- c("id", "name", "gene")

# 2. Enrichment & Visualization -------------------------------------------------------------------
writeLines("\n--------------------- Step2 Enrichment & Visualization ---------------------")
for (type in names(gene_ls)) {
  writeLines(paste(arg_ls$typeCol, type, "Enriching ..."))
  if (length(gene_ls[[type]]) < 5) {
    writeLines("too few genes (< 5), skip !")
    next
  }
  # GO BP
  EnrichFlow(
    genelist = gene_ls[[type]],
    txtFile = paste0(type, "_GOBP.txt"),
    graphFile = paste0(type, "_GOBP.png"),
    graphTi = "GO BP",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    TERM2GENE = go[go$subcat == "GO:BP", c('id', 'gene')],
    TERM2NAME = go[go$subcat == "GO:BP", c('id', 'name')]
  )
  # GO MF
  EnrichFlow(
    genelist = gene_ls[[type]],
    txtFile = paste0(type, "_GOMF.txt"),
    graphFile = paste0(type, "_GOMF.png"),
    graphTi = "GO MF",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    TERM2GENE = go[go$subcat == "GO:MF", c('id', 'gene')],
    TERM2NAME = go[go$subcat == "GO:MF", c('id', 'name')]
  )
  # GO CC
  EnrichFlow(
    genelist = gene_ls[[type]],
    txtFile = paste0(type, "_GOCC.txt"),
    graphFile = paste0(type, "_GOCC.png"),
    graphTi = "GO CC",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    TERM2GENE = go[go$subcat == "GO:CC", c('id', 'gene')],
    TERM2NAME = go[go$subcat == "GO:CC", c('id', 'name')]
  )
  # KEGG
  EnrichFlow(
    genelist = gene_ls[[type]],
    txtFile = paste0(type, "_KEGG.txt"),
    graphFile = paste0(type, "_KEGG.png"),
    graphTi = "KEGG",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    TERM2GENE = kegg[, c('id', 'gene')],
    TERM2NAME = kegg[, c('id', 'name')]
  )
}
