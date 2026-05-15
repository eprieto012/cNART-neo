# --- Libraries ---
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(pdfCluster)
library(patchwork)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(scales)
library(scDblFinder)
library(Seurat)
library(stringr)
library(tibble)
library(tidyr)
library(viridis)
# --- Shared plotting style ---
theme_ic <- theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 0.5,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "plain"),
    plot.title = element_text(hjust = 0),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )
# --- Functions ---
# 1
make_pdcd1_controls <- function(md, comp, title, n_pick = NULL, direction = c("high","low")) {
  direction <- match.arg(direction)
  sum_df <- md %>%
    filter(tag.compartment == comp, is.na(pool)) %>%
    group_by(chains_combi, cdr3_combi, v_gene_combi, j_gene_combi, cdr3_vj) %>%
    summarise(
      n_cells      = n(),
      mean_PDCD1   = mean(PDCD1_mag, na.rm = TRUE),
      median_PDCD1 = median(PDCD1_mag, na.rm = TRUE),
      max_PDCD1    = max(PDCD1_mag, na.rm = TRUE),
      .groups      = "drop"
    )
  pick_fun <- if (direction == "high") dplyr::slice_max else dplyr::slice_min
  pool_lab <- if (direction == "high") "High PDCD1 controls" else "Low PDCD1 controls"
  picked <- pick_fun(sum_df, max_PDCD1, n = n_pick, with_ties = FALSE) %>%
    mutate(.picked = TRUE)
  sum_df2 <- sum_df %>%
    left_join(picked %>% select(chains_combi, cdr3_combi, v_gene_combi, j_gene_combi, cdr3_vj, .picked),
              by = c("chains_combi","cdr3_combi","v_gene_combi","j_gene_combi","cdr3_vj")) %>%
    mutate(.picked = ifelse(is.na(.picked), FALSE, .picked))
  p <- ggplot(sum_df2, aes(x = max_PDCD1, y = n_cells)) +
    geom_point(size = 2.5, alpha = 0.25) +
    geom_point(data = subset(sum_df2, .picked), shape = 21, size = 5, stroke = 0.8, alpha = 0.75, fill = "red3") +
    labs(
      x = "Max PDCD1 expression (per clonotype)",
      y = "Number of cells (clonotype size)"
    ) +
    theme_classic() +
    ggtitle(paste0(title)) +
    theme(aspect.ratio = 1)
  tbl <- picked %>%
    mutate(
      pool   = pool_lab,
      source = comp,
      top100 = "N",
      top70  = "N",
      n_cells = as.character(n_cells)
    ) %>%
    select(pool, chains_combi, cdr3_combi, v_gene_combi, j_gene_combi,
           cdr3_vj, n_cells, source, top100, top70)
  list(plot = p, table = tbl)
}

# 2
validation_collapse_col <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste(sort(as.character(x)), collapse = " | ")
}

# 3
make_trb_pheno <- function(chains, cdr3, v, j, sep = ";") {
  ch <- str_split(chains %||% "", fixed(sep))[[1]] %>% str_trim()
  cd <- str_split(cdr3   %||% "", fixed(sep))[[1]] %>% str_trim()
  vg <- str_split(v      %||% "", fixed(sep))[[1]] %>% str_trim()
  jg <- str_split(j      %||% "", fixed(sep))[[1]] %>% str_trim()
  n <- max(length(ch), length(cd), length(vg), length(jg))
  if (n == 0) return(NA_character_)
  length(ch) <- n; length(cd) <- n; length(vg) <- n; length(jg) <- n
  idx_trb <- which(ch == "TRB")
  if (length(idx_trb) == 0) return(NA_character_)
  out <- paste0(cd[idx_trb], "_", vg[idx_trb], "_", jg[idx_trb])
  out <- out[!is.na(out) & out != "NA_NA_NA" & out != "__"]
  if (length(out) == 0) NA_character_ else paste(out, collapse = ";")
}

# 4
make_qual_palette <- function(n, palette = "Set3") {
  max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  base <- RColorBrewer::brewer.pal(max_n, palette)
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

# 5
barplot_one_single <- function(d, tcr_colors, cols_pool, cols_reactive, plot_title, stats_subtitle, x_intercept = NULL, y_intercept = NULL) {
  d_poolline <- d %>%
    dplyr::mutate(
      x_start = cell_rank - 0.5,
      x_end = cell_rank + 0.5,
      ymin = -0.065,
      ymax = -0.02
    )
  d_reactiveline <- d %>%
    dplyr::mutate(
      x_start = cell_rank - 0.5,
      x_end = cell_rank + 0.5,
      ymin = -0.13,
      ymax = -0.085
    )
  d_entpd1line <- d %>%
    dplyr::mutate(
      x_start = cell_rank - 0.5,
      x_end = cell_rank + 0.5,
      ymin = -0.195,
      ymax = -0.15
    )
  p <- ggplot() +
    geom_col(
      data = d,
      aes(x = cell_rank, y = PDCD1_mag, fill = tcr_id_plot),
      width = 1
    ) +
    scale_fill_manual(
      values = tcr_colors,
      breaks = tcr_breaks,
      drop = FALSE,
      name = "TCR.ID",
      guide = guide_legend(
        ncol = tcr_ncol,
        byrow = TRUE,
        order = 1,
        override.aes = list(alpha = 1)
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_rect(
      data = d_poolline,
      aes(xmin = x_start, xmax = x_end, ymin = ymin, ymax = ymax, fill = pool_plot),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = cols_pool,
      breaks = pool_breaks,
      drop = FALSE,
      name = "Pool",
      guide = guide_legend(
        ncol = pool_ncol,
        byrow = TRUE,
        order = 2,
        override.aes = list(alpha = 1)
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_rect(
      data = d_reactiveline,
      aes(xmin = x_start, xmax = x_end, ymin = ymin, ymax = ymax, fill = reactive_plot),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = cols_reactive,
      breaks = reactive_breaks,
      drop = FALSE,
      name = "Reactive",
      guide = guide_legend(
        ncol = reactive_ncol,
        byrow = TRUE,
        order = 3
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_rect(
      data = d_entpd1line,
      aes(xmin = x_start, xmax = x_end, ymin = ymin, ymax = ymax, fill = entpd1_plot),
      inherit.aes = FALSE
    ) +
    scale_fill_gradientn(
      name = "ENTPD1_mag",
      colours = c("grey95", "#C6B5E9", "#7A4CC2", "#2B004F"),
      limits = c(0, 0.5),
      oob = scales::squish
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(-0.205, 1),
      expand = c(0, 0),
      oob = scales::squish
    ) +
    labs(
      title = plot_title,
      subtitle = stats_subtitle,
      x = NULL,
      y = "PDCD1_mag (MAGIC)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "plain"),
      plot.subtitle = element_text(size = 10),
      aspect.ratio = 0.2,
      plot.margin = margin(t = 5.5, r = 5.5, b = 60, l = 5.5),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    coord_cartesian(clip = "off") +
    annotate("text", x = 0, y = -0.0425, label = "Pool", angle = 0, size = 3, hjust = 1, vjust = 0.5) +
    annotate("text", x = 0, y = -0.1075, label = "Reactive", angle = 0, size = 3, hjust = 1, vjust = 0.5) +
    annotate("text", x = 0, y = -0.1725, label = "ENTPD1", angle = 0, size = 3, hjust = 1, vjust = 0.5)
  if (!is.null(x_intercept)) {
    p <- p +
      geom_vline(
        xintercept = x_intercept,
        linetype = "dashed",
        linewidth = 0.25,
        colour = "black"
      )
  }
  if (!is.null(y_intercept)) {
    p <- p +
      geom_hline(
        yintercept = y_intercept,
        linetype = "dashed",
        linewidth = 0.25,
        colour = "black"
      )
  }
  p
}

# 6
make_barplot_single <- function(df, comp, threshold, top_n = NULL, tcr_colors, cols_pool, cols_reactive, pool_plot_col = "pool_plot") {
  d0 <- df %>%
    dplyr::filter(tag.compartment == comp) %>%
    dplyr::filter(!is.na(PDCD1_mag)) %>%
    dplyr::arrange(dplyr::desc(PDCD1_mag))
  if (!is.null(top_n)) {
    d0 <- dplyr::slice_head(d0, n = top_n)
  }
  d0 <- d0 %>%
    dplyr::mutate(cell_rank = dplyr::row_number())
  category_levels <- c(setdiff(names(cols_pool), "NA"), "NA")
  .mk_stats <- function(d) {
    cnt_pool <- table(factor(d[[pool_plot_col]], levels = category_levels))
    n_tcr <- dplyr::n_distinct(d$tcr_id_plot[d$tcr_id_plot != "NA"])
    list(cnt_pool = cnt_pool, n_tcr = n_tcr)
  }
  .fmt_stats <- function(s, label) {
    pool_txt <- paste0(
      names(s$cnt_pool),
      " ", as.integer(s$cnt_pool), "/", sum(s$cnt_pool),
      collapse = "; "
    )
    paste0(label, " | ", pool_txt)
  }
  d_hi <- d0 %>% dplyr::filter(PDCD1_mag >= threshold)
  d_lo <- d0 %>% dplyr::filter(PDCD1_mag < threshold)
  s_hi <- .mk_stats(d_hi)
  s_lo <- .mk_stats(d_lo)
  stats_subtitle <- paste(
    .fmt_stats(s_hi, paste0("PDCD1_mag ≥ ", signif(threshold, 4))),
    .fmt_stats(s_lo, paste0("PDCD1_mag < ", signif(threshold, 4))),
    sep = "\n"
  )
  n_above <- nrow(d_hi)
  x_intercept <- if (n_above > 0 && n_above < nrow(d0)) n_above + 0.5 else NULL
  y_intercept <- threshold
  
  barplot_one_single(
    d0,
    tcr_colors = tcr_colors,
    cols_pool = cols_pool,
    cols_reactive = cols_reactive,
    plot_title = comp,
    stats_subtitle = stats_subtitle,
    x_intercept = x_intercept,
    y_intercept = y_intercept
  )
}

# 7
plot_cluster_percentages <- function(object, 
                                     clusters = NULL, 
                                     category = NULL, 
                                     title = NULL, 
                                     colors = NULL, 
                                     horizontal = FALSE, 
                                     limits = NULL, 
                                     exclude = NULL, 
                                     remove = TRUE,
                                     print_percentages = TRUE,
                                     digits = 2) {
  data <- FetchData(object, vars = c(clusters, category))
  if (remove) {
    data <- data %>%
      filter(
        !is.na(.data[[clusters]]),
        !is.na(.data[[category]]),
        if (!is.null(exclude)) !.data[[category]] %in% exclude else TRUE
      )
  }
  if (!is.null(colors)) {
    valid_levels <- intersect(names(colors), unique(data[[category]]))
    data[[category]] <- factor(data[[category]], levels = valid_levels)
    colors <- colors[valid_levels]
  } else {
    data[[category]] <- factor(data[[category]])
  }
  counts <- data %>%
    group_by(across(all_of(c(clusters, category)))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(clusters))) %>%
    mutate(
      total_group = sum(n),
      percentage = n / total_group * 100,
      percentage_round = round(percentage, digits)
    ) %>%
    ungroup()
  if (print_percentages) {
    print(counts)
    cat("\nPiechart percentages for:", title, "\n")
    pie_percentages <- counts %>%
      arrange(.data[[clusters]], .data[[category]]) %>%
      select(
        all_of(clusters),
        all_of(category),
        n,
        total_group,
        percentage_round
      )
    print(pie_percentages)
  }
  p <- ggplot(counts, aes_string(x = clusters, y = "percentage", fill = category)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(x = clusters, y = "Percentage", fill = category, title = title) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank()
    )
  if (!is.null(limits)) {
    p <- p + ylim(limits)
  }
  if (horizontal) {
    p <- p + coord_flip()
  }
  q <- ggplot(counts, aes(x = "", y = percentage, fill = .data[[category]])) +
    geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors) +
    labs(title = title) +
    theme_void() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10)
    )
  if (!is.null(clusters) && clusters %in% colnames(counts)) {
    q <- q + facet_wrap(as.formula(paste("~", clusters)), ncol = 2)
  }
  wrap_plots(p, q, widths = c(2, 1))
}

# 8
select_representative_indices <- function(meta, n_total = 1500, tag_col = "tag") {
  stopifnot(tag_col %in% colnames(meta))
  tag <- meta[[tag_col]]
  idx_all <- seq_len(nrow(meta))
  keep <- !is.na(tag)
  tag <- tag[keep]
  idx_all <- idx_all[keep]
  if (length(idx_all) <= n_total) {
    return(idx_all)
  }
  tag_table <- table(tag)
  tag_levels <- names(tag_table)
  if (length(tag_levels) == 1) {
    return(idx_all[seq_len(min(n_total, length(idx_all)))])
  }
  prop <- as.numeric(tag_table) / sum(tag_table)
  raw_alloc <- prop * n_total
  alloc <- floor(raw_alloc)
  remainder <- n_total - sum(alloc)
  if (remainder > 0) {
    frac <- raw_alloc - alloc
    ord <- order(frac, decreasing = TRUE)
    alloc[ord[seq_len(remainder)]] <- alloc[ord[seq_len(remainder)]] + 1
  }
  alloc <- pmin(alloc, as.numeric(tag_table))
  names(alloc) <- tag_levels
  idx_by_tag <- split(idx_all, tag)
  idx_selected <- unlist(lapply(tag_levels, function(g) {
    idx_by_tag[[g]][seq_len(alloc[g])]
  }), use.names = FALSE)
  idx_selected <- sort(idx_selected)
  idx_selected
}

# 9
annotate_plot_til_pools <- function(obj, pool_col, out_col, plot_name) {
  til_pools <- grep(
    til_pattern,
    unique(na.omit(obj@meta.data[[pool_col]])),
    value = TRUE
  )
  obj@meta.data[[out_col]] <- dplyr::case_when(
    obj@meta.data[[pool_col]] %in% til_pools ~ "TIL_pool_non_tested",
    TRUE ~ "No"
  )
  obj@meta.data[[out_col]] <- dplyr::case_when(
    obj$Reactive == "Yes" ~ "TIL_pool_reactive",
    TRUE ~ obj@meta.data[[out_col]]
  )
  plot_umap <- DimPlot(
    obj,
    reduction = "umap",
    group.by = out_col,
    cols = colors,
    pt.size = ifelse(obj@meta.data[[out_col]] == "No", 0.25, 1),
    alpha = ifelse(obj@meta.data[[out_col]] == "No", 0.25, 1),
    label.size = 5
  ) +
    theme(aspect.ratio = 1)
  md <- obj@meta.data
  cells.keep.comp <- rownames(md)[!is.na(md$tag.compartment)]
  f1 <- plot_cluster_percentages(
    subset(obj, cells = cells.keep.comp),
    clusters = "tag.compartment",
    category = out_col,
    title = paste0("tag_", out_col),
    colors = colors,
    horizontal = TRUE,
    exclude = "NA",
    limits = NULL,
    remove = FALSE
  )
  cells.keep.state <- rownames(md)[!is.na(md$cellstate_manual)]
  f2 <- plot_cluster_percentages(
    subset(obj, cells = cells.keep.state),
    clusters = "cellstate_manual",
    category = out_col,
    title = paste0("cellstate_", out_col),
    colors = colors,
    horizontal = TRUE,
    exclude = "NA",
    limits = NULL,
    remove = FALSE
  )
  tcr_states <- obj@meta.data %>%
    as.data.frame() %>%
    dplyr::filter(
      .data[[out_col]] %in% c("TIL_pool_non_tested", "TIL_pool_reactive"),
      !is.na(cdr3_vj),
      cdr3_vj != ""
    ) %>%
    dplyr::count(cdr3_vj, .data[[out_col]], cellstate_manual, name = "n_cells") %>%
    dplyr::group_by(cdr3_vj, .data[[out_col]]) %>%
    dplyr::summarise(
      n_cellstates = dplyr::n(),
      cellstates_found = paste0(
        cellstate_manual, " (n=", n_cells, ")",
        collapse = " | "
      ),
      .groups = "drop"
    ) %>%
    dplyr::rename(source = .data[[out_col]]) %>%
    dplyr::arrange(cdr3_vj, source)
  plots_list[[plot_name]] <<- wrap_plots(f1, f2, plot_umap, nrow = 2)
  list(
    obj = obj,
    plot = plots_list[[plot_name]],
    tcr_states = tcr_states
  )
}