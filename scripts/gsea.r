# -------- Prog Info -------
# R version 4.1.3 (2022-03-10)
# Author: Songming Liu
# Date: 2022-12-21
# Desc: GSEA using clusterProfiler & msigdbr.
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


# Parse cmd and library packages ---------------------------------------------------------------
CheckPackage(argparse)
parser <- ArgumentParser(formatter_class='argparse.ArgumentDefaultsHelpFormatter',
                         description='GSEA using clusterProfiler & msigdbr.')

parser$add_argument('--infile', 
                    dest='infile', 
                    required=TRUE,
                    help='DEGs file, csv format.')
parser$add_argument('--outdir',
                    dest='outdir',
                    required=TRUE,
                    help='output directory.')
parser$add_argument('--gene-col',
                    dest='geneCol',
                    default='gene',
                    help='column name in DEGs file defining genes.')
parser$add_argument('--rank-col',
                    dest='rankCol',
                    default='log2fc',
                    help='column name in DEGs file used to rank genes.')
parser$add_argument('--gene-type',
                    dest='geneType',
                    default='SYMBOL',
                    choices=c('SYMBOL', 'ENTREZID'),
                    help='type of gene identifiers in DEGs file.')
parser$add_argument('--organism',
                    dest='org',
                    default='Homo-sapiens',
                    help='organism.')
arg_ls <- parser$parse_args()

## check args
if (!file.exists(arg_ls$infile)) {
  stop("infile not exist, please check!")
}

## loading required packages
CheckPackage("clusterProfiler")
CheckPackage("msigdbr")
CheckPackage("tidyverse")

## print args
writeLines("--------------------- Parameters List ---------------------")
arg_ls$org <- gsub(pattern = "-", replacement = " ", arg_ls$org)
for (para in names(arg_ls)) {
  writeLines(paste0(para, ": ", arg_ls[[para]]))
}

## create outdir
CreateDir(dir = arg_ls$outdir)

# 1. Preparation for Enrichment ------------------------------------------------
writeLines("\n--------------------- Step1 Preparation ---------------------")
## get gene list
writeLines("Infile Reading ...")
df <- read_csv(file = arg_ls$infile) %>% 
  select(gene = arg_ls$geneCol, rank_val = arg_ls$rankCol) %>% 
  arrange(-rank_val)
gene_list <- df$rank_val
names(gene_list) <- df$gene

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

## GO BP
writeLines('GO BP GSEA ...')
res <- GSEA(geneList = gene_list, 
            pvalueCutoff = 1, 
            nPermSimple = 10000,
            TERM2GENE = go[go$subcat == 'GO:BP', c('id', 'gene')],
            TERM2NAME = go[go$subcat == 'GO:BP', c('id', 'name')])
res@result %>% 
  arrange(-NES) %>% 
  write_tsv(file = sprintf('%s/gsea_go_bp.tsv', arg_ls$outdir))

## GO CC
writeLines('GO CC GSEA ...')
res <- GSEA(geneList = gene_list, 
            pvalueCutoff = 1, 
            nPermSimple = 10000,
            TERM2GENE = go[go$subcat == 'GO:CC', c('id', 'gene')],
            TERM2NAME = go[go$subcat == 'GO:CC', c('id', 'name')])
res@result %>% 
  arrange(-NES) %>% 
  write_tsv(file = sprintf('%s/gsea_go_cc.tsv', arg_ls$outdir))

## GO MF
writeLines('GO MF GSEA ...')
res <- GSEA(geneList = gene_list, 
            pvalueCutoff = 1, 
            nPermSimple = 10000,
            TERM2GENE = go[go$subcat == 'GO:MF', c('id', 'gene')],
            TERM2NAME = go[go$subcat == 'GO:MF', c('id', 'name')])
res@result %>% 
  arrange(-NES) %>% 
  write_tsv(file = sprintf('%s/gsea_go_mf.tsv', arg_ls$outdir))

## KEGG
writeLines('KEGG GSEA ...')
res <- GSEA(geneList = gene_list, 
            pvalueCutoff = 1, 
            nPermSimple = 10000,
            TERM2GENE = kegg[, c('id', 'gene')],
            TERM2NAME = kegg[, c('id', 'name')])
res@result %>% 
  arrange(-NES) %>% 
  write_tsv(file = sprintf('%s/gsea_kegg.tsv', arg_ls$outdir))

