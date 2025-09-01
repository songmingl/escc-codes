# -------- Prog Info -------
# R version 4.1.3 (2022-03-10)
# Author: Songming Liu
# Date: 2023-01-10
# Desc: R functions frequently used in ESCC projects.


# global params -----------------------------------------------------------------------------------------
prog_group_comp <- list(
  c('R-Baseline', 'R-Treat'),
  c('NR-Baseline', 'NR-Treat'),
  c('R-Baseline', 'NR-Baseline'),
  c('R-Treat', 'NR-Treat')
)

mandard_group_comp <- list(
  c('good-Baseline', 'good-Treat'),
  c('poor-Baseline', 'poor-Treat'),
  c('good-Baseline', 'poor-Baseline'),
  c('good-Treat', 'poor-Treat')
)

response_comp <- list(
  c('Responsive-Baseline', 'Responsive-Treat'),
  c('Mild-Baseline', 'Mild-Treat'),
  c('Responsive-Baseline', 'Mild-Baseline'),
  c('Responsive-Treat', 'Mild-Treat')
)
response_degree_comp <- list(
  c('severe-Baseline', 'severe-Treat'),
  c('moderate-Baseline', 'moderate-Treat'),
  c('mild-Baseline', 'mild-Treat'),
  c('severe-Baseline', 'mild-Baseline'),
  c('severe-Treat', 'mild-Treat'),
  c('moderate-Baseline', 'mild-Baseline'),
  c('moderate-Treat', 'mild-Treat'),
  c('severe-Baseline', 'moderate-Baseline'),
  c('severe-Treat', 'moderate-Treat')
)
patient_gp_comp <- list(
  c('good-Baseline', 'good-Treat'),
  c('pEKI-Baseline', 'pEKI-Treat'),
  c('pEKD-Baseline', 'pEKD-Treat'),
  c('good-Baseline', 'pEKD-Baseline'),
  c('good-Treat', 'pEKD-Treat'),
  c('pEKI-Baseline', 'pEKD-Baseline'),
  c('pEKI-Treat', 'pEKD-Treat'),
  c('good-Baseline', 'pEKI-Baseline'),
  c('good-Treat', 'pEKI-Treat')
)

gp_comp_map <- list(
  `group` = prog_group_comp,
  `mandard_group` = mandard_group_comp, 
  `response` = response_comp,
  `response_degree` = response_degree_comp,
  `patient_gp` = patient_gp_comp,
  `patient_gp_v2` = patient_gp_comp
)
gp_comp_map_pre <- lapply(X = gp_comp_map, FUN = function(subls) {
   subls[
     unlist(lapply(X = subls, FUN = function(comp){all(str_ends(string = comp, pattern = 'Baseline'))}))
   ]
})
gp_comp_diff_map <- list(
  `group` = list(c('R', 'NR')),
  `mandard_group` = list(c('good', 'poor')),
  `response` = list(c('Responsive', 'Mild')),
  `response_degree` = list(
    c('severe', 'mild'),
    c('moderate', 'mild'),
    c('severe', 'moderate')
  ),
  `patient_gp` = list(
    c('good', 'pEKI'),
    c('good', 'pEKD'),
    c('pEKI', 'pEKD')
  ),
  `patient_gp_v2` = list(
    c('good', 'pEKI'),
    c('good', 'pEKD'),
    c('pEKI', 'pEKD')
  )
)
gp_lvls <- list(
  `mandard_group` = c('good', 'poor'),
  `response` = c('Responsive', 'Mild'),
  `response_degree` = c('severe', 'moderate', 'mild'),
  `group` = c('R', 'NR'),
  `patient_gp` = c('good', 'pEKD', 'pEKI', 'poor'),
  `patient_gp_v2` = c('good', 'pEKD', 'pEKI', 'poor')
)


clin_cols <- c('patient', 'treatment', 'prognosis', 'group', 
               'treatment_group', 'mandard_score', 'mandard_group', 
               'tumor_cell_remain')
clin_cols_full <- c(
  'patient', 'treatment', 'prognosis', 'group', 
  'treatment_group', 'mandard_score', 'mandard_group', 
  'tumor_cell_remain', 'stage',  'response_degree',
  'mandard_group_correct', 'treatment_sum', 'response'
)
wes_clin_cols <- c(
  'patient_id', 'prognosis', 'group', 'treatment_group', 
  'mandard_score', 'mandard_group', 'tumor_cell_remain', 
  'response_degree', 'stage'
)

onco_colors <- list(
    group = c(R='#EE0000FF', NR='#3B4992FF'),
    response_degree = c(severe='#D43F3AFF', moderate='#EEA236FF', mild='#5CB85CFF'),
    mandard_group = c(good='#EE0000FF', poor='#3B4992FF'),
    ek_group = c(R='#D43F3AFF', NR_highEK='#5CB85CFF', NR_lowEK='#46B8DAFF'),
    sample_type = c(Baseline='#1F77B4FF', Treat='#2CA02CFF'),
    treatment_sum = c(chemo='#749B58FF', `chemo+immune`='#F0E685FF')
)

scanpy_pal <- readRDS('/nfsshare2/nfshome/ming/projects/ESCC/assets/sc_palettes.rds')

# transform coordinate using logNp
lognp_trans <- function(base = 10, n = 1) {
  scales::trans_new(name='lognp', 
                    transform = function(x){log(base = base, x + n)},
                    inverse = function(x){base^x - n})
}

# set jupyter notebook figure option --------------------------------------------------------------
jp_opt <- function(res = 150, wd=5, hg=4) {
  options(repr.plot.res = res, repr.plot.width = wd, repr.plot.height = hg)
}

# add clinical info -----------------------------------------------------------------------------
add_clin_info <- function(df, columns=clin_cols_full, ftsv=NULL, fxlsx=NULL, merge_by='patient') {
  if (! is.null(fxlsx)) {
    clin <- readxl::read_xlsx(fxlsx, sheet='clinical_info')
  } else if (! is.null(ftsv)) {
    clin <- read_tsv(ftsv, show_col_types = F)
  }
  if (is.null(columns)) {
    columns <- colnames(clin)
  }
  uniq_cols <- setdiff(columns, colnames(df))
  loginfo('these clinial info will be added: %s', paste(uniq_cols, collapse=', '))
  # add merge_by columns
  uniq_cols <- unique(c(merge_by, uniq_cols))
  df <- clin %>%
    dplyr::select(all_of(uniq_cols)) %>%
    merge(df, by=merge_by, all.y = T)
  return(df)
}

# sankey plot ---------------------------------------------------------------------------------------
sankey_plot <- function(df, fill, outfile, mapping, fwidth = 6, fheight = 7) {
  p <- ggplot(data = df, mapping = mapping) + 
    ggalluvial::geom_alluvium(aes(fill = .data[[fill]])) +
    ggalluvial::geom_stratum(width = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_hue(l = 45) +
    labs(fill = '')
  if (is.null(outfile)) {
    return(p)
  } else {
    ggsave(filename = outfile, plot = p, width = fwidth, height = fheight)
  }
}

# Heatmap -------------------------------------------------------------------------------------------
treat_heatmap <- function(mat,
                          name,
                          col_title='',
                          return_ht = FALSE,
                          outfile = NULL,
                          fig_size = c(9, 7),
                          ...) {
  # treatment compare heatmap
  # Args:
  #   mat: heatmap matrix.
  #   name: matrix legend name.
  #   col_title: column title.
  #   return_ht: return the heatmap object.
  #   outfile: output file, if NULL, return_ht will set to TRUE.
  #   fig_size: output heatmap size, c(width, height)
  # Returns:
  #   complexheatmap object
  hp <- Heatmap(matrix = mat,
                name = name,
                row_names_side = 'left', 
                show_column_names = F,
                column_title = col_title,
                ...) 
  if (is.null(outfile)) {
    return(hp)
  } else if (tools::file_ext(outfile) == 'pdf') {
    pdf(file = outfile, width = fig_size[1], height = fig_size[2])
    ComplexHeatmap::draw(hp)
    dev.off()
  } else {
    png(file = outfile, width = fig_size[1], height = fig_size[2], units = 'in', res = 300)
    ComplexHeatmap::draw(hp)
    dev.off()
  }
}

# DE --------------------------------------------------------------------------------------------------
init_logging <- function(logfile = NULL, force = TRUE) {
  logReset()
  if (is.null(logfile)) {
    addHandler(handler = writeToConsole, logger = '')
  } else if (force) {
    if (file.exists(logfile)) {
      unlink(logfile)
    }
    addHandler(handler = writeToFile, logger = '', file = logfile)
  } else {
    addHandler(handler = writeToFile, logger = '', file = logfile)
  }
}


do_deseq2 <- function(exprs, treat, ref, cal_pct = FALSE, min_samp_per_gp = 3, parallel = TRUE) {
  # Do DESeq2 to find DEGs 
  # Args:
  #   exprs: raw counts matrix, gene (row) x sample (column)
  #   treat: samples used as treat group
  #   ref: samples used as ref group
  #   min_samp_per_gp: min samples each group should have.
  #   parallel: do DESeq parallel.
  #   cal_pct: calculate expression percent of each gene in treat/ref group.
  # Returns:
  #   a data frame of DEGs or NULL
  require(logging)
  require(DESeq2)
  require(tidyverse)
  # process samples
  t_len <- length(treat)
  treat <- intersect(colnames(exprs), treat)
  r_len <- length(ref)
  ref <- intersect(colnames(exprs), ref)
  loginfo('%g/%g valid samples in treat, %g/%g valid samples in ref', length(treat), t_len, length(ref), r_len)
  if (min(length(treat), length(ref)) < min_samp_per_gp) {
    loginfo('samples not enough for DE (at least %g samples per group), skip', min_samp_per_gp)
    return(NULL)
  }
  
  # exclude genes not expressed in any sample
  n_gene_raw <- nrow(exprs)
  samp_cnt_per_gene <- rowSums(exprs[, c(treat, ref)] > 0)
  used_genes <- names(samp_cnt_per_gene)[samp_cnt_per_gene > 0]
  loginfo('%g/%g genes used for analysis (other genes are expressed in neither `treat` nor `ref`)',
          length(used_genes), n_gene_raw)
  
  # DESeq2: prepare exprs matrix and coldata
  coldata <- data.frame(sample = c(treat, ref), 
                        group = c(rep('treat', length(treat)),
                                  rep('ref', length(ref)))) %>% 
    column_to_rownames('sample')
  coldata$group <- relevel(as.factor(coldata$group), ref = 'ref')  # set treat vs ctrl (ref)
  subexprs <- exprs[used_genes, rownames(coldata)] # ensure sample name in same orde
  
  # DESeq2: run
  dds <- DESeqDataSetFromMatrix(countData = subexprs,
                                colData = coldata,
                                design= ~ group)
  dds <- DESeq(dds, parallel = parallel)
  # escape adj.p=NA: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  res <- results(dds, independentFiltering = F)
  loginfo('DESeq2 finish!')
  
  # obtain results
  df_res <- as.data.frame(res) %>%
    rownames_to_column('symbol') 
  
  # calcualte gene epxression percent
  if (cal_pct) {
    df_res <- data.frame(
      pct_treat = rowSums(exprs[, treat] > 0) / length(treat),
      pct_ref = rowSums(exprs[, ref] > 0) / length(ref)
    ) %>%
    rownames_to_column('symbol') %>%
    merge(df_res, by = 'symbol')
  }
  
  df_res <- df_res %>%
    arrange(-log2FoldChange)

  return(df_res)
}

volcano_plot <- function(df,
                         gene,
                         pval,
                         log2fc,
                         color_by,
                         p_cut = NULL,
                         log2fc_cut = NULL,
                         labels = NULL,
                         pt_size = 2,
                         col_map = c(`up-regulated` = '#EE0000', 
                                     `down-regulated` = '#3B4992',
                                     `unchanged` = 'gray50',
                                     `unexpressed` = '#79ae98'),
                         ...) {
  # valcano plot for DEGs
  # Args:
  #   df: DEGs results, data frame
  #   gene: column in df indicates gene
  #   pval: column in df indicates signif
  #   log2fc: column in df indicates log2-transformed fold changes
  #   color_by: DEG types, used to color point.
  #   p_cut: cutoff of pval to determine DEGs, used for adding cutoff line
  #   log2fc_fc: cutoff of log2fc to determine DEGs, used for adding cutoff line
  #   labels: genes to label the point.
  #   pt_size: ponit size
  # Returns:
  #   a ggplot object
  require(tidyverse)
  require(ggpubr)
  require(ggrepel)
  require(logging)
  df <- df %>% 
    mutate(log_pval = -log10(.data[[pval]]))
  
  # main plot
  p <- ggscatter(data = df, x = log2fc, y = 'log_pval', color = color_by,
                 alpha = 0.6, size = pt_size) +
    scale_color_manual(values = col_map) +
    labs(color = '', y = str_glue('-log10({pval})'), x = 'log2 fold change')
  # add threshold line
  if (!is.null(p_cut)) {
    p <- p + geom_hline(yintercept = -log10(p_cut), linetype = 'dashed')
  }
  if (!is.null(log2fc_cut)) {
    p <- p + geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = 'dashed')
  }
  
  # add label
  if (!is.null(labels)) {
    subdf <- df %>% 
      filter(.data[[gene]] %in% labels)
    if (nrow(subdf) > 0) {
      p <- p +
        ggrepel::geom_text_repel(data = subdf, aes(label = .data[[gene]]), min.segment.length = 0.001)
    }
    lost_genes <- setdiff(labels, subdf[[gene]])
    if (length(lost_genes) > 0) {
      loginfo('%g genes not found: %s', length(lost_genes), paste(lost_genes, collapse = ', '))
    }
  }
  return(p)
}

# scatter with reg line -------------------------------------------------------------------------------
list.get <- function(x, name, default = NA) {
  # get list element and return a default
  return(ifelse(test = (name %in% names(x)), yes = x[[name]], no = default))
}

scatter_with_reg <- function(df,
                             x,
                             y,
                             color = 'black',
                             add_reg = TRUE,
                             add_cor = TRUE,
                             reg_pos = c(0.05, 0.95),
                             reg_hjust = 0,
                             cor_method = 'pearson',
                             cor_pos = c(0.05, 0.85),
                             cor_hjust = 0,
                             ...) {
  
  p <- ggscatter(data = df, x = x, y = y, color=color, 
                 size = list.get(x = list(...), name = 'size', default = 3),
                 alpha = list.get(x = list(...), name = 'alpha', default = 0.6),
                 palette = list.get(x = list(...), name = 'palette', default = 'nejm'),
                 add = 'reg.line')
  # add reg line labels
  if (color %in% colnames(df) & add_reg) {
    p <- p +
      stat_regline_equation(
        aes(
          color = .data[[color]],
          label = paste(after_stat(adj.rr.label), '*","~', after_stat(eq.label))
        ),
        label.x.npc = reg_pos[1],
        label.y.npc = reg_pos[2],
        hjust = reg_hjust
      )
  } else if (add_reg) {
    p <- p +
      stat_regline_equation(
        aes(label = paste(after_stat(adj.rr.label), '*","~', after_stat(eq.label))),
        color = color,
        label.x.npc = reg_pos[1],
        label.y.npc = reg_pos[2],
        hjust = reg_hjust
      )
  }
  # add correlation labels
  if (color %in% colnames(df) & add_cor) {
    p <- p +
      stat_cor(
        aes(
          color = .data[[color]],
          label = paste(stringr::str_to_title(cor_method), '==', after_stat(r))
        ),
        method = cor_method, 
        label.x.npc = cor_pos[1],
        label.y.npc = cor_pos[2],
        hjust = cor_hjust
      )
  } else if (add_cor) {
    p <- p +
      stat_cor(
        aes(label = paste(stringr::str_to_title(cor_method), '==', after_stat(r))),
        color = color,
        method = cor_method, 
        label.x.npc = cor_pos[1],
        label.y.npc = cor_pos[2],
        hjust = cor_hjust
      )
  }
  return(p)
}

# scatter with fit -------------------------------------------------------------------------------------
fit.lm <- function(x, y, df) {
  fit_summary <- summary.lm(lm(data = df, formula = str_glue('{y} ~ {x}')))
  return(list(slope = fit_summary$coefficients[x, 'Estimate'],
              intercept = fit_summary$coefficients['(Intercept)', 'Estimate'],
              r2 = fit_summary$r.squared,
              p = pf(fit_summary$fstatistic[1], 
                     fit_summary$fstatistic[2],
                     fit_summary$fstatistic[3], 
                     lower.tail = F)))
}

scatter_with_fit <- function(df, x, y, n_signif = 2, pt_size = 4, add_corr = TRUE, corr_method = 'pearson', ...) {
  vals <- fit.lm(x = x, y = y, df = df)
  lbl_info <- substitute(expr = Y==slope %.% X + intercept*','~italic(R^2)==r2*','~italic(P)==p,
                         env = lapply(vals, function(x){ifelse(is.numeric(x), signif(x, n_signif), x)}))
  p <- ggplot(df, aes_string(x = x, y = y, ...)) +
    geom_point(size = pt_size, alpha = 0.6) +
    geom_abline(slope = vals$slope, intercept = vals$intercept, color = 'red') +
    annotate(geom = 'text', 
             x = max(df[[x]]) * 0.5 + min(df[[x]]) * 0.5,
             y = max(df[[y]]) * 0.75 + min(df[[y]]) * 0.25,
             color = 'red', 
             label = lbl_info)
  if (add_corr) {
    cor_coef <- cor(x = df[[x]], df[[y]], method = corr_method)
    p <- p +
      annotate(geom = 'text',
               x = max(df[[x]]) * 0.2 + min(df[[x]]) * 0.8,
               y = max(df[[y]]) * 0.65 + min(df[[y]]) * 0.35,
               color = 'red',
               label = str_glue('Pearson: {signif(cor_coef, n_signif)}'))
  }
  return(p)
}

scatter_with_fit2 <- function(df, x, y, color_by, n_signif = 2, pt_size = 4, corr_method = 'pearson') {
  # Fit lm for each color_by
  method <- str_to_title(corr_method)
  
  # set color for each color_by group
  gps <- sort(unique(df[[color_by]]))
  col_pal <- scales::hue_pal()(length(gps))
  names(col_pal) <- gps
  
  # base plot
  p <- df %>% 
    ggplot(aes_string(x, y, color = color_by)) +
    geom_point(size = pt_size, alpha = 0.6) +
    scale_color_manual(values = col_pal)
  
  # add lines & obtain lm line info
  anno_ls <- c()
  for (gp in gps) {
    subdf <- df %>% 
      filter(.data[[color_by]] == gp)
    vals <- fit.lm(x = x, y = y, df = subdf) 
    vals <- lapply(vals, function(x){ifelse(is.numeric(x), signif(x, n_signif), x)})
    # add line
    p <- p + 
      geom_abline(slope = vals$slope, intercept = vals$intercept, color = col_pal[gp])
    # annotation info
    sign_str <- ifelse(vals$intercept < 0, '-', '+')
    corr <- signif(cor(x = subdf[[x]], subdf[[y]], method = corr_method), n_signif)
    anno_ls <- c(anno_ls, str_glue('Y=={vals$slope}%.%X{sign_str}{abs(vals$intercept)}*","~italic(R^2)=={vals$r2}*","~{method}=={corr}'))
  }
  
  # add annotation
  lbl_x <- max(df[[x]]) * 0.1 + min(df[[x]]) * 0.9
  lbl_y_min <- max(df[[y]]) * 0.75 + min(df[[y]]) * 0.25
  lbl_y_max <- max(df[[y]]) * 0.85 + min(df[[y]]) * 0.15
  p <- p + 
    annotate(geom = 'text', 
             label = anno_ls, 
             x = lbl_x, 
             y = seq(lbl_y_min, lbl_y_max, length = length(gps)), 
             parse = T, hjust = 0, color = col_pal[gps])
  return(p)
}


# TCR sharing network plot ----------------------------------------------------------------------------------------
exhau_cyto_netplot2 <- function(df_node, df_edge, pt_fill, edge_color, edge_weight='weight') {
  # Transition Network plot show TCR sharing and cell type exhuastion/cytotocity score.
  # Args:
  #   df_node: data frame of network nodes, include columns: node, exhau, cyto, ...
  #   df_edge: data frame of network edges, include columns: start_exhau, start_cyto, end_exhau, end_cyto, edge_type, ...
  #   pt_fill: color network nodes by this column in df_node
  #   edge_color: color network edge by this column in df_edge
  #   edge_weight: column in df_edge indicates edge wight (usually the TCR sharing score)
  # Returns:
  #   a ggplot object.
  require(ggrepel)
  p <- ggplot(df_node, aes(exhau, cyto))
  # add edge
  edges <- df_edge %>%
    pull(edge_type) %>%
    unique()
  for (e in edges) {
    tmp <- df_edge %>% 
      filter(edge_type == e) %>% 
      distinct()
    if (nrow(tmp) > 1) {
      stop(sprintf('edge: %s more one records', e))
    }
    p <- p + 
      geom_segment(data = tmp, 
                   show.legend = T,
                   aes(x = `start_exhau`, y = `start_cyto`, xend = `end_exhau`, yend = `end_cyto`, color = .data[[edge_color]]), 
                   linewidth = tmp[[edge_weight]] * 40,
                   alpha = 0.8)
  }
  # add node
  p <- p + 
    geom_point(aes(fill = .data[[pt_fill]]), size = 8, alpha = 0.8, shape = 21) +
    geom_text_repel(aes(label = node), show.legend = F, size = 6) + 
    labs(x = 'Exhaustion Score', y = 'Cytotoxicity Score') 
  return(p)
}

cal_edge_weight_diff <- function(a, b) {
  # a - b 
  cols_by <- a %>%
    select(-weight) %>%
    colnames()
  df <- merge(a, b, by = cols_by, all = T) %>%
    mutate(weight.x = replace_na(weight.x, 0),
           weight.y = replace_na(weight.y, 0)) %>%
    mutate(weight.diff = weight.x - weight.y) %>%
    mutate(weight = abs(weight.diff)) %>%
    filter(weight != 0)
  loginfo('a(%g) - b(%g) => %g diff edges.', nrow(a), nrow(b), nrow(df))
  return(df)
}

exhau_cyto_netplot <- function(df_node, df_edge, color_by, edge_weight='weight', text_size = 6, pt_size = 8) {
  # Transition Network plot show TCR sharing and cell type exhuastion/cytotocity score.
  # Args:
  #   df_node: data frame of network nodes, include columns: node, exhau, cyto, ...
  #   df_edge: data frame of network edges, include columns: start_exhau, start_cyto, end_exhau, end_cyto, edge_type, ...
  #   color_by: color network nodes by this column in df_node
  #   edge_weight: column in df_edge indicates edge wight (usually the TCR sharing score)
  # Returns:
  #   a ggplot object.
  require(ggrepel)
  p <- ggplot(df_node, aes(exhau, cyto))
  # add edge
  edges <- df_edge %>%
    pull(edge_type) %>%
    unique()
  for (e in edges) {
    tmp <- df_edge %>% 
      filter(edge_type == e) %>% 
      distinct()
    if (nrow(tmp) > 1) {
      stop(sprintf('edge: %s more one records', e))
    }
    p <- p + 
      geom_segment(data = tmp, 
                   aes(x = `start_exhau`, y = `start_cyto`, xend = `end_exhau`, yend = `end_cyto`), 
                   linewidth = tmp[[edge_weight]] * 40,
                   alpha = 0.8,
                   color = 'gray')
  }
  # add node
  p <- p + 
    geom_point(aes(color = .data[[color_by]]), size = pt_size, alpha = 0.8) +
    geom_text_repel(aes(label = node), show.legend = F, size = text_size) + 
    labs( x = 'Exhaustion Score', y = 'Cytotoxicity Score')
  return(p)
}


# heatmap to show batch effect ------------------------------------------------------------------------------------

batch_heatmap <- function(data, sample, cluster, outfile, 
                          width, height, show=FALSE, digit=1){
  # Heatmap to check batch effector.
  #
  # Args:
  #
  #   data: data frame.
  #   sample: one column in data indicating sample info.
  #   cluster: one column in data indicating cluster info.
  #   outfile: graph output file path.
  #   width: figure width.
  #   height: figure height.
  #   show: show figure.
  #   digit: digit of percentage.
  require(tidyverse)
  require(patchwork)
  require(ggsci)
  df <- as.matrix(table(data[[sample]], data[[cluster]]))
  # percent of clusters in each sample
  dp <- reshape2::melt(round(df/rowSums(df) * 100, digit))
  dp$Var2 <- factor(dp$Var2, levels = sort(unique(dp$Var2)))
  dp$Var1 <- factor(dp$Var1, levels = sort(unique(dp$Var1)))
  p1 <- ggplot(data = dp, aes(x = Var2, 
                              y = Var1, 
                              fill = value,
                              label = value)) + 
    geom_tile() + 
    geom_text() + 
    labs(x = cluster, y = sample, fill = 'percent', 
         title = sprintf('percent of %s in %s', cluster, sample))
  # percent of samples in each cluster
  dp <- reshape2::melt(round(t(df)/colSums(df) * 100, digit))
  dp$Var1 <- factor(dp$Var1, levels = sort(unique(dp$Var1)))
  dp$Var2 <- factor(dp$Var2, levels = sort(unique(dp$Var2)))
  p2 <- ggplot(data = dp, aes(x = Var1, 
                              y = Var2, 
                              fill = value, 
                              label = value)) + 
    geom_tile() + 
    geom_text() + 
    labs(x = cluster, y = sample, fill = 'percent', 
         title = sprintf('percent of %s in %s', sample, cluster)) 
  p <- (p1 | p2) &
    scale_x_discrete(position = 'top') &
    scale_fill_gsea() &
    theme(plot.title = element_text(hjust = 0),
          axis.text.x = element_text(angle = 60, hjust = 0.1, vjust = 0.1))
  ggsave(filename = outfile, width = width, height = height, plot = p)
  if (show) {
    print(p)
  }
}


# functions for cell compostion -----------------------------------------------------------------------------------

cal_cell_comp <- function(df_info, samp_cnt, clinical_cols=c('patient', 'prognosis', 'treatment', 'group')) {
  # Calculate cell composition by patient & sample_type.
  # Args:
  #   df_info: dataframe, merged obs & clinical info of specific cell type.
  #   samp_cnt: dataframe, cell counts of each sample, at least 2 columns: sample, n_cell_samp
  # Returns:
  #   a dataframe with calculated cell compositions.
  
  # count frequency of each cell type
  df <- df_info %>% 
    count(patient, sample_type, cell_type, name='freq')
  # add lost cell types (0 cells in the cell type of one sample)
  pat_samp <- df %>%
    select(patient, sample_type) %>%
    distinct()
  all_ctypes <- df %>%
    select(cell_type) %>%
    distinct() %>%
    pull(cell_type)
  for (i in 1:nrow(pat_samp)) {
    id_pat <- pat_samp[['patient']][i]
    id_stype <- pat_samp[['sample_type']][i]
    exist_ctypes <- df %>%
      filter(patient == id_pat & sample_type == id_stype) %>%
      distinct() %>%
      pull(cell_type)
    lost_ctypes <- setdiff(all_ctypes, exist_ctypes)
    for (ctype in lost_ctypes) {
      df <- df %>%
        add_row(patient = id_pat,
                sample_type = id_stype,
                cell_type = ctype, 
                freq = 0)
    }
  }
  # fraction of cell types: devided by current total cell number of a sample (by major cell type)
  df <- df %>%
    group_by(patient, sample_type) %>%
    mutate(pct = 100 * freq / sum(freq))
  # fraction of cell types: devided by total cells of a sample
  df <- ungroup(df) %>% 
    unite(col = 'sample', patient, sample_type, remove=FALSE, sep='-') %>% 
    merge(samp_cnt, by = 'sample') %>% 
    mutate(pct_by_total_cell = 100 * freq / n_cell_samp)
  # add clinical info (one patient one recoerd)
  df <- df_info %>%
    select(all_of(clinical_cols)) %>%
    distinct() %>%
    merge(df, by = 'patient')
  return(df)
}

cell_comp_boxplot <- function(df, 
                              x=c('group', 'sample_type'),
                              y='pct', 
                              xorder=c('R-Baseline', 'R-Treat', 
                                       'NR-Baseline', 'NR-Treat', 
                                       'Out-Baseline', 'Out-Treat',
                                       'NA-Baseline', 'NA-Treat'),
                              pt_fill='group', 
                              pair_by='patient',
                              fill_order=c('R', 'NR', 'Out', 'NA'), 
                              facet_by='cell_type', 
                              wrap_free='free',
                              ncol=4,
                              xangle = 45,
                              xtitle = '',
                              ytitle = 'percentage',
                              seed_i = 0,
                              ...) {
  # Boxplot of Cell Composition.
  # Args:
  #   df: dataframe indicating cell composition info.
  #   xorder: x lab orders.
  #   fill_order: fill order.
  #   x: one or more column names in dataframe, defining the x lab, multi-columns unite by '-'.
  #   y: column name of df defining cell compostion.
  #   pt_fill: column of df to fill point color by (scale_fill).
  #   pair_by: column of df indicating paired points which will be linked by line in final figure.
  #   facet_by: column of df to facet by. 
  #   wrap_free: facet scales to wrap free.
  #   ncol: number of columns when facet the plot.
  # Returns:
  #   a ggplot object
  require(ggplot2)
  require(ggsci)
  # set x, legend order
  df <- unite(data = df, col = 'plot_x', all_of(x), sep = '-', remove = FALSE)
  if (is.null(xorder)) {
    xorder <- sort(unique(df$plot_x))
  }
  df$plot_x <- factor(x = df$plot_x,
                      levels = intersect(xorder, df$plot_x))
  if ((!is.null(fill_order)) & (!is.null(pt_fill))) {
    df[[pt_fill]] <- factor(x = df[[pt_fill]],
                            levels = intersect(fill_order, df[[pt_fill]]))
  }
  
  # add jitter
  set.seed(seed_i)
  psudo_count <- runif(nrow(df), -0.3, 0.3)
  df <- df %>% 
    mutate(jitter_x = as.numeric(plot_x) + psudo_count)
  # plot
  p <- ggplot(data = df, aes(x = plot_x, y = .data[[y]])) + 
    geom_boxplot(outlier.shape = NA, size = 1, fatten = 1.5, 
                 width = 0.6, position = position_dodge(1))
  if (!is.null(pt_fill)) {
    p <- p + 
      geom_point(aes(x = jitter_x, fill = .data[[pt_fill]]), shape = 21, alpha = 0.7,
                 size = list.get(x = list(...), name = 'size', default = 3)) + 
      scale_fill_d3(na.value = 'grey50')
  } else {
    p <- p + 
      geom_point(aes(x = jitter_x), fill = 'grey50', shape = 21, alpha = 0.7,
                 size = list.get(x = list(...), name = 'size', default = 3))
  }
  # add pair line
  if (!is.null(pair_by)) {
    p <- p + 
      geom_line(mapping = aes(x = jitter_x, 
                              y = .data[[y]], 
                              group = .data[[pair_by]]),
                color = 'gray50')
  }
  # facet by
  if (!is.null(facet_by)) {
    p <- p + 
      facet_wrap(~ .data[[facet_by]], scales = wrap_free, ncol = ncol)
  }
  # plot style
  p <- p +
    labs(x = xtitle, y = ytitle) + 
    theme(axis.text.x = element_text(angle = xangle, hjust = 0.95))
  return(p)
}

cal_ratio2malig <- function(malig, cell_comp, pseudo = 0, adds = 0.05) {
  # Calculate cell type percent to malignant cells + cell type cells.
  # Args:
  #   malig: data frame of malignant cell stat, at least two column: sample, n_malig
  #   cell_comp: data frame of cell compositions.
  #   pseudo: pseudo cells to assign as malignant cells if the malignant cells of a sample is 0. 
  #           (pseudo < 1 means fraction) 
  #   adds: adds cells to add to malignant cells of a sample (adds < 0 means fraction).
  # Returns:
  #    updated cell composition dataframe, add 2 new columns: n_malig, ratio2malig
  res <- cell_comp %>%
    mutate(sample = paste(patient, sample_type, sep = '-')) %>%
    merge(malig, by = 'sample', all.x = T)
  if (pseudo < 1 && pseudo > 0) {
    res <- res %>%
      mutate(n_malig = if_else(condition = is.na(n_malig),
                               true = as.double(pseudo * n_cell_samp),
                               false = as.double(n_malig)))
  } else if (pseudo >= 0) {
    res <- res %>%
      mutate(n_malig = if_else(condition = is.na(n_malig),
                               true = as.double(pseudo),
                               false = as.double(n_malig)))
  }
  if (adds < 1 && adds > 0) {
    res <- res %>% 
      mutate(n_malig = adds * n_cell_samp + n_malig)
  } else if (adds >= 0) {
    res <- res %>%
      mutate(n_malig = adds + n_malig)
  }
  res <- res %>%
    mutate(ratio2malig = round(freq / (n_malig + freq), 2))
  return(res)
}


# degs filtering ------------------------------------------------------------------------------------------------------
create_dir <- function(dir) {
  # Check and Create dirs
  if (!dir.exists(dir)) {
    dir.create(path = dir, recursive = T)
  }
}


flt_gene <- function(infile, outfile, 
                     max_adj_p = 0.01, 
                     min_porp = 0.25, 
                     min_log2fc = 1, 
                     max_genes = 500) {
  # filtr to obtain final DEGs for scanpy results.
  # Args:
  #   infile: csv file, at least these columns: symbol, log2fc, pval, pval_adj, prop1, prop2.
  #   outfile: tsv file of filtered results, contaions a new column: deg_type.
  #   max_adj_p: max value of adjust p value (pval_adj).
  #   min_prop: min proporion the gene expressed in either group (prop1, prop2).
  #   min_log2fc: min log2-fold changes (log2fc).
  #   max_genes: max number of up/down-regulate genes, sort by log2fc.
  # Returns:
  #   No returns.
  require(tidyverse)
  writeLines(sprintf('thresholds: max_adj_p=%g, min_porp=%g, min_log2fc=%g, max_genes=%g',
                     max_adj_p, min_porp, min_log2fc, max_genes))
  
  # filter
  df <- read_csv(infile, show_col_types = F) %>% 
    filter(abs(log2fc) >= min_log2fc & pmax(prop1, prop2) >= min_porp & pval_adj <= max_adj_p) %>% 
    mutate(deg_type = case_when(
      log2fc > 0 ~ 'up',
      log2fc < 0 ~ 'down',
      log2fc == 0 ~ 'unchange'
    )) %>% 
    arrange(-log2fc)
  ct <- as.list(table(df$deg_type))
  writeLines(sprintf("before selecting: up: %d, down: %d", ct$up, ct$down))
  
  # only top max_degs
  df <- df %>% 
    group_by(deg_type) %>% 
    top_n(n = max_genes, wt = abs(log2fc))
  ct <- as.list(table(df$deg_type))
  writeLines(sprintf("after selecting: up: %d, down: %d", ct$up, ct$down))
  
  write_tsv(x = df, file = outfile)
}

flt_gene_seu <- function(infile, outfile, 
                         max_adj_p = 0.01, 
                         min_porp = 0.25, 
                         min_log2fc = 1, 
                         max_genes = 500) {
  # filtr to obtain final DEGs for seurat results
  # Args:
  #   infile: csv file, at least these columns: symbol, avg_log2FC, p_val, p_val_adj, pct.1, pct.2
  #   outfile: tsv file of filtered results, contaions a new column: deg_type.
  #   max_adj_p: max value of adjust p value (p_val_adj).
  #   min_prop: min proporion the gene expressed in either group (pct.1, pct.2).
  #   min_log2fc: min log2-fold changes (avg_log2FC).
  #   max_genes: max number of up/down-regulate genes, sort by avg_log2FC
  # Returns:
  #   No returns.
  require(tidyverse)
  writeLines(sprintf('thresholds: max_adj_p=%g, min_porp=%g, min_log2fc=%g, max_genes=%g',
                     max_adj_p, min_porp, min_log2fc, max_genes))
  
  # filter
  df <- read_csv(infile, show_col_types = F) %>% 
    filter(abs(avg_log2FC) >= min_log2fc & pmax(pct.1, pct.2) >= min_porp & p_val_adj <= max_adj_p) %>% 
    mutate(deg_type = case_when(
      avg_log2FC > 0 ~ 'up',
      avg_log2FC < 0 ~ 'down',
      avg_log2FC == 0 ~ 'unchange'
    )) %>% 
    arrange(-avg_log2FC)
  ct <- as.list(table(df$deg_type))
  writeLines(sprintf("before selecting: up: %d, down: %d", ct$up, ct$down))
  
  # only top max_degs
  df <- df %>% 
    group_by(deg_type) %>% 
    top_n(n = max_genes, wt = abs(avg_log2FC))
  ct <- as.list(table(df$deg_type))
  writeLines(sprintf("after selecting: up: %d, down: %d", ct$up, ct$down))
  
  write_tsv(x = df, file = outfile)
}

# data transform --------------------------------------------------------------------------------------------------

create_seurat <- function(count_dir,
                          fmeta = NULL,
                          fhvg = NULL,
                          fpca = NULL,
                          fumap = NULL,
                          pca_key = 'PC_',
                          umap_key = 'UMAP_') {
  # Create Seurat object manually.
  # Args:
  #   count_dir: 10x-like mtx file directory.
  #   fmeta: meta data file (with header, first column is cell).
  #   fhvg: hvg file (no header, one row one gene).
  #   fpca: pca coordinates file (with header, first column should indicate cell).
  #   fumap: UMAP coordinates file (with header, first column should indicate cell).
  #   pca_key: prefix of PC names in fpca, e.g. "PC_".
  #   umap_key: prefix of UMAP names in fumap, e.g. "UMAP_".
  # Return:
  #   Seurat Object.
  require(tidyverse)
  require(Seurat)
  require(logging)
  # create basical seurat object
  loginfo('creating seurat object ...')
  scrna <- Read10X(count_dir) %>% 
    CreateSeuratObject(min.cells = 0, min.features = 0)
  # add meta data
  loginfo('add meta data ...')
  if (! is.null(fmeta)) {
    meta <- read_csv(file = fmeta, show_col_types = F) %>%
      as.data.frame()
    rownames(meta) <- meta[[1]]
    meta <- meta[, 2:length(meta)]
    tmp <- scrna@meta.data %>% 
      merge(meta, by = 'row.names', all.x = T) %>%
      column_to_rownames('Row.names')
    scrna@meta.data <- tmp[colnames(scrna), ]
  }
  # add variable features
  loginfo('add variable features ...')
  if (! is.null(fhvg)) {
    VariableFeatures(scrna) <- read_csv(fhvg, col_names = F, show_col_types = F)$X1
  }
  # add pca coordinates
  loginfo('add PCA coordinates ...')
  if (! is.null(fpca)) {
    coords <- read_csv(file = fpca, show_col_types = F) %>%
      as.data.frame()
    rownames(coords) <- coords[[1]]
    coords <- coords[, 2:length(coords)] %>%
      as.matrix()
    scrna@reductions$pca <- CreateDimReducObject(embeddings = coords, assay = 'RNA', key = pca_key)
  }
  # add umap coordinates
  loginfo('add UMAP coordinates ...')
  if (! is.null(fumap)) {
    coords <- read_csv(file = fumap, show_col_types = F) %>%
      as.data.frame()
    rownames(coords) <- coords[[1]]
    coords <- coords[, 2:length(coords)] %>%
      as.matrix()
    scrna@reductions$umap <- CreateDimReducObject(embeddings = coords, assay = 'RNA', key = umap_key)
  }
  return(scrna)
}

# monocle2 --------------------------------------------------------------------------------------------------------

anova_test <- function(obj, by, genes = NULL, ncores = 30) {
  # ANOVA by given group in metadata for seurat object.
  # Args:
  #   obj: seurat object.
  #   by: column in meta.data to group cells by.
  #   genes: genes to calculate, if NULL use all genes.
  #   ncores: n cores to use parallel.
  # Return:
  #   a data frame with F statistics, p, adjust p, and average expression of each gene.
  if (is.null(genes)) {
    used_genes <- rownames(obj)
  } else {
    used_genes <- intersect(rownames(obj), genes)
    logging::loginfo('%g/%g genes in the object', length(used_genes), length(genes))
  }
  
  # ANOVA
  anova_ls <- parallel::mclapply(used_genes, mc.cores = ncores, FUN = function(gene) {
    df <- obj@assays$RNA@data[gene, ] %>%
      as.data.frame() %>%
      rename(g_exprs = '.') 
    df[rownames(obj@meta.data), 'cell_group'] <- obj@meta.data[[by]]
    tmp <- anova(aov(g_exprs ~ cell_group, data = df))['cell_group', ] %>%
      mutate(symbol = gene) %>% 
      select(symbol, Fval = 'F value', p = 'Pr(>F)')
    rownames(tmp) <- NULL
    # calculate average expression
    avg_expr <- AverageExpression(obj, group.by = by, features = gene, slot = 'data')[['RNA']]
    tmp <- cbind.data.frame(tmp, avg_expr)
    return(tmp)
  })
  res <- do.call(what = rbind.data.frame,
                 args = anova_ls)
  res[['p.adj']] <- p.adjust(res$p, method = 'BH')
  return(res)
}
