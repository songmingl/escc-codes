# R functions can frequenctly used

# Create an output directory if needed. The historical notebooks call this
# helper but it was not bundled with the original helper file.
create_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

# Historical notebooks populate these list objects panel by panel. Initialize
# them here so the notebooks also run in a clean R session.
if (!exists("gp_lvls", inherits = FALSE)) gp_lvls <- list()
if (!exists("gp_comp_map", inherits = FALSE)) gp_comp_map <- list()
if (!exists("gp_comp_diff_map", inherits = FALSE)) gp_comp_diff_map <- list()
if (!exists("gp_comp_map_pre", inherits = FALSE)) gp_comp_map_pre <- list()

BatchHeatmap <- function(data, sample, cluster, outfile, width, height, 
                          show=FALSE, digit = 1) {
    # Heatmap to check batch effector.
    #
    # Args:
    #
    #   data: data frame.
    #   sample: one column in data indicating sample info.
    #   cluster: one column in data indicating cluster info.
    #   outfile: graph output file path.
    #   width: figure width.
    #   heigth: figure height.
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

# cell composition calculating ----------------------------------------

MergeClinialObs <- function(fclin, fobs, by = 'sample') {
    # Merge clinical and obs information.
    # Args:
    #   fclin: clinical info file.
    #   fobs: scanpy obs file.
    #   by: column to merge by.
    # Returns:
    #   merge data frame
    require(tidyverse)
    # clinical info
    clinic <- read_tsv(fclin, show_col_types = F)
    # cell cluster info
    info <- read_csv(fobs, show_col_types = F) %>% 
        rename(barcode = `...1`) %>%
        merge(clinic, by = by)
    return(info)
}

MergeClinialObs2 <- function(df_clin, fobs, by) {
    # Merge clinical and obs information.
    # Args:
    #   df_clin: clinical info data frame.
    #   fobs: scanpy obs file.
    #   by: column to merge by.
    # Returns:
    #   merge data frame
    require(tidyverse)
    # cell cluster info
    info <- read_csv(fobs, show_col_types = F) %>% 
        rename(barcode = `...1`) %>%
        separate(col = sample, into = c('patient', 'sample_type'), sep = '-') %>%
        mutate(sample = paste(patient, sample_type, sep = '-')) %>%
        merge(df_clin, by = by)
    return(info)
}

CellCompBoxplot <- function(df, cond, frac, ctype = 'cell_type', patient = 'patient', ncol=4) {
    # Plot Cell Composition using Boxplot, for Paired Samples.
    # Args:
    #   df: data frame to plot.
    #   cond: column of df indicating conditions, e.g. sample_type (baseline vs treat).
    #   frac: column of df indicating cell fraction.
    #   ctype: column of df indicating cell type.
    #   patient: column of df indicating patient id.
    # Returns:
    #   a ggplot object.
    require(ggplot2)
    require(ggpubr)
    p <- ggplot(data = df, aes(x = .data[[cond]], 
                               y = .data[[frac]], 
                               fill = .data[[cond]])) +
        geom_boxplot(outlier.size = 0) + 
        geom_point(size = 2, alpha = 0.5) + # add point
        geom_line(aes(group = .data[[patient]]), color = 'gray', alpha = 0.8) + # add line between paired samples
        facet_wrap(~ .data[[ctype]], ncol = ncol, scales = 'free') + # facet by cell type
        stat_compare_means(method = 'wilcox.test', paired = T) + 
        theme(legend.position = 'none')
    return(p)
}

CellCompBoxplot2 <- function(df, cond, frac, ctype = 'cell_type', pair_by = NULL, ncol=4) {
    # Plot Cell Composition using Boxplot, for Both Paired & Unpaired Samples.
    # Args:
    #   df: data frame to plot.
    #   cond: column of df indicating conditions, e.g. sample_type (baseline vs treat).
    #   frac: column of df indicating cell fraction.
    #   ctype: column of df indicating cell type.
    #   pair_by: column of df indicating paired info, e.g. patient id.
    # Returns:
    #   a ggplot object.
    require(ggplot2)
    require(ggpubr)
    p <- ggplot(data = df, aes(x = .data[[cond]], 
                               y = .data[[frac]], 
                               fill = .data[[cond]])) +
        geom_boxplot(outlier.size = 0) +  
        geom_point(size = 2, alpha = 0.5) # add point
    if (is.null(pair_by)) {
        p <- p + 
            facet_wrap(~ .data[[ctype]], ncol = ncol, scales = 'free') + # facet by cell type
            stat_compare_means(method = 'wilcox.test')
    } else {
        p <- p +
            geom_line(aes(group = .data[[pair_by]]), color = 'gray', alpha = 0.8) + # add line between paired samples
            facet_wrap(~ .data[[ctype]], ncol = ncol, scales = 'free') + # facet by cell type
            stat_compare_means(method = 'wilcox.test', paired = T) 
    }
    p <- p + theme(legend.position = 'none')
    return(p)
}


# Reproducibility helpers used by the bundled Figure notebooks -----------------
jp_opt <- function(res = 150, wd = 5, hg = 4) {
    options(repr.plot.res = res, repr.plot.width = wd, repr.plot.height = hg)
}

add_clin_info <- function(df, columns = NULL, ftsv = NULL, fxlsx = NULL, merge_by = 'patient') {
    if (!is.null(fxlsx)) {
        clin <- readxl::read_xlsx(fxlsx, sheet = 'clinical_info')
    } else if (!is.null(ftsv)) {
        clin <- readr::read_tsv(ftsv, show_col_types = FALSE)
    } else {
        stop('Provide ftsv or fxlsx for clinical information.')
    }
    if (is.null(columns)) columns <- colnames(clin)
    cols <- unique(c(merge_by, setdiff(columns, colnames(df))))
    clin <- dplyr::select(clin, dplyr::all_of(cols))
    merge(clin, df, by = merge_by, all.y = TRUE)
}

cell_comp_boxplot <- function(df, x = c('group', 'sample_type'), y = 'pct',
                              xorder = NULL, pt_fill = 'group', pair_by = 'patient',
                              fill_order = NULL, facet_by = 'cell_type', wrap_free = 'free',
                              ncol = 4, xangle = 45, xtitle = '', ytitle = 'percentage',
                              seed_i = 0, ...) {
    df <- tidyr::unite(df, col = 'plot_x', dplyr::all_of(x), sep = '-', remove = FALSE)
    if (is.null(xorder)) xorder <- sort(unique(df$plot_x))
    df$plot_x <- factor(df$plot_x, levels = intersect(xorder, unique(df$plot_x)))
    if (!is.null(pt_fill) && !is.null(fill_order)) {
        df[[pt_fill]] <- factor(df[[pt_fill]], levels = intersect(fill_order, unique(df[[pt_fill]])))
    }
    set.seed(seed_i)
    df$jitter_x <- as.numeric(df$plot_x) + stats::runif(nrow(df), -0.3, 0.3)
    dots <- list(...); point_size <- if (!is.null(dots$size)) dots$size else 3
    p <- ggplot2::ggplot(df, ggplot2::aes(x = plot_x, y = .data[[y]])) +
        ggplot2::geom_boxplot(outlier.shape = NA, size = 1, fatten = 1.5, width = 0.6)
    if (!is.null(pt_fill)) {
        p <- p + ggplot2::geom_point(ggplot2::aes(x = jitter_x, fill = .data[[pt_fill]]), shape = 21, alpha = 0.7, size = point_size) +
            ggsci::scale_fill_d3(na.value = 'grey50')
    } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(x = jitter_x), fill = 'grey50', shape = 21, alpha = 0.7, size = point_size)
    }
    if (!is.null(pair_by)) {
        p <- p + ggplot2::geom_line(ggplot2::aes(x = jitter_x, y = .data[[y]], group = .data[[pair_by]]), color = 'gray50')
    }
    if (!is.null(facet_by)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste('~', facet_by)), scales = wrap_free, ncol = ncol)
    p + ggplot2::labs(x = xtitle, y = ytitle) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xangle, hjust = 0.95))
}

# Equivalent network plot for Figure 3b. The original helper was not bundled.
exhau_cyto_netplot <- function(df_node, df_edge, color_by = 'node_type', text_size = 4.5, pt_size = 6) {
    edge_cols <- intersect(c('from', 'to', 'weight'), names(df_edge))
    edge_data <- dplyr::select(df_edge, dplyr::all_of(edge_cols))
    node_name <- if ('node' %in% names(df_node)) 'node' else names(df_node)[1]
    node_data <- dplyr::rename(df_node, name = dplyr::all_of(node_name))
    required <- unique(c(edge_data$from, edge_data$to))
    missing <- setdiff(required, node_data$name)
    if (length(missing) > 0) {
        node_data <- dplyr::bind_rows(node_data, data.frame(name = missing, node_type = 'others'))
    }
    graph <- igraph::graph_from_data_frame(edge_data, vertices = node_data, directed = FALSE)
    ggraph::ggraph(graph, layout = 'fr') +
        ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), alpha = 0.65, colour = 'grey50') +
        ggraph::geom_node_point(ggplot2::aes_string(color = color_by), size = pt_size) +
        ggraph::geom_node_text(ggplot2::aes(label = name), size = text_size, repel = TRUE) +
        ggraph::scale_edge_width(range = c(0.3, 2.5), guide = 'none') +
        ggplot2::theme_void()
}


# Figure 2/4 helpers retained from the original project helper library.
lognp_trans <- function(base = 10, n = 1) {
    scales::trans_new(name = 'lognp', transform = function(x) log(x + n, base = base), inverse = function(x) base^x - n)
}

treat_heatmap <- function(mat, name, col_title = '', return_ht = FALSE, outfile = NULL, fig_size = c(9, 7), ...) {
    hp <- ComplexHeatmap::Heatmap(matrix = mat, name = name, row_names_side = 'left', show_column_names = FALSE, column_title = col_title, ...)
    if (is.null(outfile) || return_ht) return(hp)
    if (tools::file_ext(outfile) == 'pdf') {
        grDevices::pdf(outfile, width = fig_size[1], height = fig_size[2]); ComplexHeatmap::draw(hp); grDevices::dev.off()
    } else {
        grDevices::png(outfile, width = fig_size[1], height = fig_size[2], units = 'in', res = 300); ComplexHeatmap::draw(hp); grDevices::dev.off()
    }
    invisible(hp)
}


cal_cell_comp <- function(df_info, samp_cnt, clinical_cols = c('patient', 'prognosis', 'treatment', 'group')) {
    df <- dplyr::count(df_info, patient, sample_type, cell_type, name = 'freq')
    pat_samp <- dplyr::distinct(dplyr::select(df, patient, sample_type))
    all_ctypes <- dplyr::pull(dplyr::distinct(dplyr::select(df, cell_type)), cell_type)
    for (i in seq_len(nrow(pat_samp))) {
        id_pat <- pat_samp$patient[i]; id_stype <- pat_samp$sample_type[i]
        existing <- dplyr::pull(dplyr::filter(df, patient == id_pat, sample_type == id_stype), cell_type)
        for (ctype in setdiff(all_ctypes, existing)) {
            df <- dplyr::bind_rows(df, data.frame(patient = id_pat, sample_type = id_stype, cell_type = ctype, freq = 0))
        }
    }
    df <- df %>% dplyr::group_by(patient, sample_type) %>% dplyr::mutate(pct = 100 * freq / sum(freq)) %>% dplyr::ungroup() %>%
        tidyr::unite(col = 'sample', patient, sample_type, remove = FALSE, sep = '-') %>%
        merge(samp_cnt, by = 'sample') %>% dplyr::mutate(pct_by_total_cell = 100 * freq / n_cell_samp)
    clinical_cols <- intersect(clinical_cols, names(df_info))
    clinical <- df_info %>% dplyr::select(dplyr::all_of(clinical_cols)) %>% dplyr::distinct()
    merge(clinical, df, by = 'patient')
}


# Notebook 交互式统计核查 ------------------------------------------------------

format_panel_p <- function(x, digits = 2) {
    ifelse(
        x < 0.001,
        formatC(x, format = 'e', digits = digits),
        formatC(x, format = 'f', digits = 3)
    )
}

panel_wilcox_stats <- function(data, x, y, comparisons, facet_by = NULL,
                               use_paired = FALSE, pair_by = NULL,
                               paired_comparisons = NULL,
                               apply_correction = FALSE,
                               correction_method = 'BH',
                               alternative = 'two.sided') {
    stopifnot(is.data.frame(data), length(x) >= 1, length(y) == 1)
    required <- unique(c(x, y, facet_by, pair_by))
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0) {
        stop('Missing columns: ', paste(missing_cols, collapse = ', '))
    }

    df <- data
    df$.panel_group <- if (length(x) == 1) {
        as.character(df[[x]])
    } else {
        do.call(paste, c(lapply(df[x], as.character), sep = '-'))
    }
    df$.panel_y <- df[[y]]
    if (is.null(facet_by)) {
        df$.panel_facet_key <- '__all__'
    } else {
        df$.panel_facet_key <- do.call(
            paste,
            c(lapply(df[facet_by], as.character), sep = '\r')
        )
    }

    comparison_key <- function(z) paste(z, collapse = '\r')
    paired_keys <- if (is.null(paired_comparisons)) {
        vapply(comparisons, comparison_key, character(1))
    } else {
        vapply(paired_comparisons, comparison_key, character(1))
    }

    results <- list()
    result_i <- 0L
    for (facet_key in unique(df$.panel_facet_key)) {
        facet_data <- df[df$.panel_facet_key == facet_key, , drop = FALSE]
        facet_values <- if (is.null(facet_by)) NULL else facet_data[1, facet_by, drop = FALSE]

        for (comparison_i in seq_along(comparisons)) {
            comparison <- comparisons[[comparison_i]]
            if (length(comparison) != 2) stop('Each comparison must contain two groups.')
            comparison_data <- facet_data[
                facet_data$.panel_group %in% comparison & !is.na(facet_data$.panel_y),
                , drop = FALSE
            ]
            paired_this_comparison <- isTRUE(use_paired) &&
                comparison_key(comparison) %in% paired_keys

            if (paired_this_comparison) {
                if (is.null(pair_by)) {
                    stop('pair_by is required when use_paired = TRUE.')
                }
                comparison_data$.pair_id <- if (length(pair_by) == 1) {
                    as.character(comparison_data[[pair_by]])
                } else {
                    do.call(paste, c(lapply(comparison_data[pair_by], as.character), sep = '\r'))
                }
                duplicate_pairs <- duplicated(comparison_data[c('.pair_id', '.panel_group')])
                if (any(duplicate_pairs)) {
                    stop('Paired Wilcoxon requires one value per pair and group. Duplicate pair/group rows found.')
                }
                wide <- tidyr::pivot_wider(
                    comparison_data[c('.pair_id', '.panel_group', '.panel_y')],
                    names_from = '.panel_group', values_from = '.panel_y'
                )
                wide <- tidyr::drop_na(wide, dplyr::all_of(comparison))
                x1 <- wide[[comparison[1]]]
                x2 <- wide[[comparison[2]]]
                n1 <- n2 <- nrow(wide)
                n_pairs <- nrow(wide)
            } else {
                x1 <- comparison_data$.panel_y[comparison_data$.panel_group == comparison[1]]
                x2 <- comparison_data$.panel_y[comparison_data$.panel_group == comparison[2]]
                n1 <- length(x1)
                n2 <- length(x2)
                n_pairs <- NA_integer_
            }
            if (length(x1) == 0 || length(x2) == 0) next

            test <- stats::wilcox.test(
                x1, x2, paired = paired_this_comparison,
                alternative = alternative
            )
            result_i <- result_i + 1L
            row <- data.frame(
                group1 = comparison[1],
                group2 = comparison[2],
                comparison_order = comparison_i,
                paired = paired_this_comparison,
                n1 = n1,
                n2 = n2,
                n_pairs = n_pairs,
                statistic_W = unname(test$statistic),
                p_raw = test$p.value,
                stringsAsFactors = FALSE
            )
            if (!is.null(facet_by)) row <- cbind(facet_values, row)
            results[[result_i]] <- row
        }
    }
    if (length(results) == 0) stop('No valid Wilcoxon comparisons were produced.')

    stat_df <- dplyr::bind_rows(results)
    stat_df$p_adj <- if (isTRUE(apply_correction)) {
        stats::p.adjust(stat_df$p_raw, method = correction_method)
    } else {
        stat_df$p_raw
    }
    stat_df$p_display <- stat_df$p_adj
    stat_df$label <- format_panel_p(stat_df$p_display)

    position_data <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c('.panel_facet_key')))) |>
        dplyr::summarise(
            .panel_ymin = min(.panel_y, na.rm = TRUE),
            .panel_ymax = max(.panel_y, na.rm = TRUE),
            .groups = 'drop'
        )
    stat_df$.panel_facet_key <- if (is.null(facet_by)) {
        '__all__'
    } else {
        do.call(paste, c(lapply(stat_df[facet_by], as.character), sep = '\r'))
    }
    stat_df <- dplyr::left_join(stat_df, position_data, by = '.panel_facet_key') |>
        dplyr::group_by(.panel_facet_key) |>
        dplyr::arrange(comparison_order, .by_group = TRUE) |>
        dplyr::mutate(
            .panel_span = dplyr::if_else(
                is.finite(.panel_ymax - .panel_ymin) & (.panel_ymax - .panel_ymin) > 0,
                .panel_ymax - .panel_ymin,
                pmax(abs(.panel_ymax), 1)
            ),
            y.position = .panel_ymax + (0.08 + 0.12 * (dplyr::row_number() - 1)) * .panel_span
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.panel_facet_key, -.panel_ymin, -.panel_ymax, -.panel_span)

    stat_df
}

add_panel_wilcox <- function(plot, stat_df, hide_ns = FALSE, tip_length = 0.01) {
    plot + ggpubr::stat_pvalue_manual(
        stat_df,
        xmin = 'group1', xmax = 'group2',
        y.position = 'y.position', label = 'label',
        hide.ns = hide_ns, tip.length = tip_length,
        inherit.aes = FALSE
    )
}

configure_panel_wilcox <- function(plot, ...) {
    stat_df <- panel_wilcox_stats(data = plot$data, ...)
    configured_plot <- add_panel_wilcox(plot, stat_df)
    attr(configured_plot, 'panel_wilcox_stats') <- stat_df
    configured_plot
}

get_panel_wilcox_stats <- function(plot) {
    attr(plot, 'panel_wilcox_stats')
}
