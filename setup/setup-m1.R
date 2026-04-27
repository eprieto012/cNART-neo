# --- Python ---
library(reticulate)
use_virtualenv("/Users/endikaprieto/venv", required = TRUE)
py_config()
# --- Libraries ---
library(data.table)
library(djvdj)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggridges)
library(keras3)
library(patchwork)
library(pdfCluster)
library(purrr)
library(RColorBrewer)
library(scGate)
reticulate::import("magic")
library(Rmagic)
library(readr)
library(scDblFinder)
library(Seurat)
library(SoupX)
library(stringr)
library(tensorflow)
library(tidyr)
library(Trex)
library(viridis)
# --- Set seed for R and python ---
set.seed(123)
py_run_string("
import numpy as np
import random
np.random.seed(123)
random.seed(123)
")
# --- Plotting style ---
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
plot_tcr_frequencies <- function(meta,
                                 cdr3_col, v_gene_col, j_gene_col, chains_col,
                                 samples = c("Sample1", "Sample2", "Sample3", "Sample4"),
                                 chain_filter = c("TRA;TRA;TRB", "TRA;TRB"),
                                 palette_name = "turbo",
                                 save_plot = NULL,
                                 save_table = NULL,
                                 save_shared_table = NULL,
                                 plot_width = 20,
                                 plot_height = 15,
                                 export_decimals_comma = TRUE,
                                 export_digits = 4,
                                 y_log10 = FALSE,
                                 y_pseudocount = 1e-3) {
  meta <- meta %>%
    dplyr::mutate(full_tcr = paste(.data[[cdr3_col]],
                                   .data[[v_gene_col]],
                                   .data[[j_gene_col]],
                                   sep = "_")) %>%
    dplyr::select(tag, full_tcr, chains = dplyr::all_of(chains_col)) %>%
    dplyr::filter(!is.na(full_tcr), full_tcr != "NA_NA_NA") %>%
    dplyr::filter(chains %in% chain_filter) %>%
    dplyr::filter(tag %in% samples)
  freq_table <- meta %>%
    dplyr::group_by(tag, full_tcr) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(tag) %>%
    dplyr::mutate(freq = count / sum(count) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(tag, dplyr::desc(freq))
  clones <- lapply(samples, function(s) unique(freq_table$full_tcr[freq_table$tag == s]))
  names(clones) <- samples
  pairwise_intersections <- combn(
    samples, 2, simplify = FALSE,
    FUN = function(x) length(intersect(clones[[x[1]]], clones[[x[2]]]))
  )
  freq_wide <- freq_table %>%
    dplyr::select(tag, full_tcr, freq) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = tag, values_from = freq)
  shared_table <- freq_wide %>%
    dplyr::mutate(
      n_present = apply(dplyr::across(dplyr::all_of(samples)), 1, function(x) sum(!is.na(x))),
      combo = apply(dplyr::across(dplyr::all_of(samples)), 1, function(x) {
        present <- samples[!is.na(x)]
        paste(present, collapse = " & ")
      })
    ) %>%
    dplyr::filter(n_present >= 2) %>%
    dplyr::select(combo, full_tcr, dplyr::all_of(samples)) %>%
    dplyr::arrange(nchar(combo), combo, full_tcr)
  shared_clones_use <- shared_table$full_tcr %>% unique()
  n_shared <- length(shared_clones_use)
  legend_title <- "Shared clonotypes (>=2 samples)"
  if (n_shared == 0) {
    freq_table_plot <- freq_table %>%
      dplyr::mutate(color_group = factor("Other", levels = "Other"))
    palette <- c(Other = "grey70")
  } else {
    cols <- if (requireNamespace("viridisLite", quietly = TRUE) &&
                palette_name %in% c("viridis", "magma", "inferno", "plasma", "cividis", "turbo")) {
      do.call(getFromNamespace(palette_name, "viridisLite"), list(n_shared))
    } else if (requireNamespace("viridisLite", quietly = TRUE)) {
      viridisLite::turbo(n_shared)
    } else {
      grDevices::hcl.colors(n_shared, palette = "Turbo")
    }
    names(cols) <- shared_clones_use
    freq_table_plot <- freq_table %>%
      dplyr::mutate(color_group = dplyr::if_else(full_tcr %in% shared_clones_use, full_tcr, "Other"))
    freq_table_plot$color_group <- factor(freq_table_plot$color_group,
                                          levels = c(shared_clones_use, "Other"))
    palette <- c(cols, Other = "grey70")
  }
  if (isTRUE(y_log10)) {
    freq_table_plot <- freq_table_plot %>%
      tidyr::complete(tag = samples, full_tcr, fill = list(count = 0, freq = 0)) %>%
      dplyr::mutate(
        freq = dplyr::if_else(freq <= 0, y_pseudocount, freq),
        color_group = dplyr::if_else(full_tcr %in% shared_clones_use, as.character(full_tcr), "Other"),
        color_group = factor(color_group, levels = c(shared_clones_use, "Other"))
      )
  }
  p <- ggplot2::ggplot(freq_table_plot,
                       ggplot2::aes(x = tag, y = freq, group = full_tcr, color = color_group)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::geom_line(alpha = 0.6) +
    ggplot2::scale_color_manual(values = palette, name = legend_title, drop = FALSE) +
    ggplot2::labs(y = "Frequency (%)", x = "Sample") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1.5)
  if (isTRUE(y_log10)) {
    p <- p +
      ggplot2::scale_y_log10(
        limits = c(1e-3, 100),
        breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100),
        labels = c("0.001", "0.01", "0.1", "1", "10", "100")
      )
  }
  if (!is.null(save_plot)) {
    ggplot2::ggsave(filename = save_plot, plot = p, width = plot_width, height = plot_height)
  }
  if (!is.null(save_table)) {
    readr::write_tsv(freq_table_plot, save_table, na = "")
  }
  if (!is.null(save_shared_table)) {
    to_export <- shared_table
    if (isTRUE(export_decimals_comma)) {
      to_export <- to_export %>%
        dplyr::mutate(
          dplyr::across(where(is.numeric),
                        ~ format(round(.x, export_digits),
                                 decimal.mark = ",",
                                 scientific = FALSE,
                                 trim = TRUE))
        )
    }
    readr::write_tsv(to_export, save_shared_table, na = "")
  }
  return(list(freq_table = freq_table_plot,
              shared_table = shared_table,
              plot = p,
              pairwise_intersections = pairwise_intersections,
              shared_clones_plot = shared_clones_use))
}
# 2
make_top_alluvial <- function(df, ref_col, ref_name, techniques, top = 10, save_path = NULL) {
  top10_per_sample <- df %>%
    dplyr::filter(!is.na(.data[[ref_col]]), .data[[ref_col]] != "NA_NA_NA") %>%
    group_by(sample, clonotype = .data[[ref_col]]) %>%
    summarise(freq_ref = n(), .groups = "drop") %>%
    arrange(sample, desc(freq_ref)) %>%
    group_by(sample) %>%
    slice_head(n = top) %>%
    ungroup()
  alluvial_df <- df %>%
    select(sample, all_of(techniques)) %>%
    pivot_longer(cols = all_of(techniques),
                 names_to = "technique",
                 values_to = "clonotype_value") %>%
    semi_join(top10_per_sample, by = c("sample", "clonotype_value" = "clonotype")) %>%
    group_by(sample, technique, clonotype_value) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    do({
      sample_name <- unique(.$sample)
      top10 <- top10_per_sample %>% filter(sample == sample_name) %>% pull(clonotype)
      tidyr::complete(., 
                      technique = techniques, 
                      clonotype_value = top10, 
                      fill = list(count = 0),
                      sample = sample_name)
    }) %>%
    ungroup() %>%
    group_by(sample, technique) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup() %>%
    mutate(technique = factor(technique, levels = techniques))
  plots <- lapply(unique(alluvial_df$sample), function(samp) {
    df_samp <- alluvial_df %>% filter(sample == samp)
    ggplot(df_samp,
           aes(x = technique, stratum = clonotype_value,
               alluvium = clonotype_value, y = freq, fill = clonotype_value)) +
      geom_flow(stat = "alluvium", lode.guidance = "forward") +
      geom_stratum(color = "black", size = 0.1) +
      scale_fill_viridis_d(option = "turbo") +
      theme_minimal() +
      theme(legend.position = "right", aspect.ratio = 1) +
      labs(title = samp,
           y = "relative frequency",
           x = "technique",
           fill = paste0(ref_name, " top 10 clonotypes"))
  })
  combined_plot <- wrap_plots(plots, ncol = 2) + plot_annotation(title = paste(ref_name, "top10"))
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 30, height = 15)
    message("Plot saved to: ", save_path)
  }
  return(combined_plot)
}
# 3
quietRibogenes <- function(genes) {
  intersect(genes, na.omit(scGate::genes.blacklist.default$Hs$Ribo))
}