suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (arg in args) {
    if (grepl("=", arg, fixed = TRUE)) {
      key <- sub("^--", "", sub("=.*$", "", arg))
      val <- sub("^--[^=]*=", "", arg)
      out[[key]] <- val
    }
  }
  out
}

args <- parse_args()
input_dir <- if (!is.null(args$input_dir)) args$input_dir else "input"
output_base <- if (!is.null(args$output_base)) args$output_base else "."

pcol <- if (!is.null(args$pcol)) args$pcol else "pvalue"
pthresh <- if (!is.null(args$pthresh)) as.numeric(args$pthresh) else 0.01
lfc_thresh <- if (!is.null(args$lfc)) as.numeric(args$lfc) else 1

if (!pcol %in% c("pvalue", "padj")) {
  stop("pcol must be 'pvalue' or 'padj'")
}

find_input <- function(pattern) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) != 1) {
    stop(sprintf("Expected 1 file for pattern %s, found %d", pattern, length(files)))
  }
  files[[1]]
}

path_em <- find_input("em\\.csv$")
path_ann <- find_input("annotations\\.csv$")
path_ss <- find_input("sample_sheet\\.csv$")

em <- read.delim(path_em, sep = "\t", header = TRUE, check.names = FALSE)
ann <- read.delim(path_ann, sep = "\t", header = TRUE, check.names = FALSE)
ss <- read.delim(path_ss, sep = "\t", header = TRUE, check.names = FALSE)

if (!"ID" %in% names(em)) {
  stop("Expected an 'ID' column in expression matrix")
}

if (!all(c("SAMPLE", "SAMPLE_GROUP") %in% names(ss))) {
  stop("sample_sheet must have SAMPLE and SAMPLE_GROUP columns")
}

if (!all(c("Gene ID", "Associated Gene Name") %in% names(ann))) {
  stop("annotations.csv must have 'Gene ID' and 'Associated Gene Name' columns")
}

# Prepare count matrix
rownames(em) <- em$ID
em$ID <- NULL

samples <- ss$SAMPLE
missing_samples <- setdiff(samples, colnames(em))
if (length(missing_samples) > 0) {
  stop(sprintf("Missing samples in em matrix: %s", paste(missing_samples, collapse = ", ")))
}

count_data <- em[, samples, drop = FALSE]
count_data <- round(as.matrix(count_data))

col_data <- data.frame(
  row.names = ss$SAMPLE,
  group = factor(ss$SAMPLE_GROUP, levels = c("gut", "duct", "node"))
)

# Ensure columns match rownames of col_data and apply sample order
col_data <- col_data[order(col_data$group, rownames(col_data)), , drop = FALSE]
count_data <- count_data[, rownames(col_data), drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ group
)

dds <- DESeq(dds)

res_with_symbols <- function(contrast) {
  res <- results(dds, contrast = contrast)
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)

  ann_idx <- match(res_df$GeneID, ann[["Gene ID"]])
  gene_symbol <- ann[["Associated Gene Name"]][ann_idx]
  empty <- is.na(gene_symbol) | gene_symbol == ""
  gene_symbol[empty] <- res_df$GeneID[empty]
  res_df$GeneSymbol <- gene_symbol
  res_df
}

contrasts <- list(
  gut_duct = c("group", "gut", "duct"),
  duct_node = c("group", "duct", "node"),
  node_gut = c("group", "node", "gut")
)

res_list <- lapply(contrasts, res_with_symbols)

output_step2a <- file.path(output_base, "output_step2a")
dir.create(output_step2a, showWarnings = FALSE, recursive = TRUE)

write_res_table <- function(res_df, label) {
  df <- res_df
  df <- df[is.finite(df$pvalue) & is.finite(df$padj) & is.finite(df$log2FoldChange), ]

  out <- data.frame(
    GeneID = df$GeneSymbol,
    log2FoldChange = df$log2FoldChange,
    pvalue = df$pvalue,
    padj = df$padj,
    stringsAsFactors = FALSE
  )

  out <- out[order(out$padj, out$pvalue, na.last = TRUE), ]

  out_path <- file.path(output_step2a, sprintf("%s.tsv", label))
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", out_path)
}

write_res_table(res_list$gut_duct, "de_gut_duct")
write_res_table(res_list$duct_node, "de_duct_node")
write_res_table(res_list$node_gut, "de_node_gut")

output_plots <- if (pcol == "padj") {
  file.path(output_base, "output_step2c")
} else {
  file.path(output_base, "output_step2b")
}

dir.create(output_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_plots, "volcano"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_plots, "ma"), showWarnings = FALSE, recursive = TRUE)

vst_data <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_data)
scaled_mat <- t(scale(t(vst_mat)))
scaled_mat[is.na(scaled_mat)] <- 0

# Map rownames to symbols for display
symbol_map <- res_list$gut_duct$GeneSymbol
names(symbol_map) <- res_list$gut_duct$GeneID
symbol_map <- symbol_map[rownames(scaled_mat)]
missing_symbol <- is.na(symbol_map) | symbol_map == ""
symbol_map[missing_symbol] <- rownames(scaled_mat)[missing_symbol]
rownames(scaled_mat) <- symbol_map

# PCA plot (scaled values)
gene_var <- apply(scaled_mat, 1, var)
pca_mat <- scaled_mat[gene_var > 0, , drop = FALSE]
pca <- prcomp(t(pca_mat), center = FALSE, scale. = FALSE)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = col_data$group
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 50) +
  theme_minimal() +
  labs(title = "PCA (scaled)", x = "PC1", y = "PC2")

ggsave(file.path(output_plots, "pca_samples.png"), pca_plot, width = 7, height = 5, dpi = 300)

# Heatmap of top DE genes (scaled values)
sig_mask <- lapply(res_list, function(x) x[[pcol]] < pthresh & abs(x$log2FoldChange) > lfc_thresh)

sig_genes <- unique(unlist(lapply(seq_along(res_list), function(i) {
  res_list[[i]]$GeneSymbol[sig_mask[[i]]]
})))

min_p_all <- Reduce(pmin, lapply(res_list, function(x) x[[pcol]]))
names(min_p_all) <- res_list$gut_duct$GeneSymbol
min_p_all <- min_p_all[is.finite(min_p_all)]

min_p <- min_p_all
if (length(sig_genes) > 0) {
  min_p <- min_p[names(min_p) %in% sig_genes]
}
min_p <- min_p[order(min_p, na.last = TRUE)]

top_n <- 50
top_genes <- names(min_p)[seq_len(min(top_n, length(min_p)))]

heatmap_mat <- scaled_mat[top_genes, , drop = FALSE]

annotation_col <- data.frame(Group = col_data$group)
rownames(annotation_col) <- rownames(col_data)

png(file.path(output_plots, "heatmap_top_degenes.png"), width = 1200, height = 1600, res = 150)
pheatmap(
  heatmap_mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 6,
  main = "Top DE genes (scaled)"
)
dev.off()

plot_volcano <- function(res_df, label, pcol, out_dir) {
  df <- res_df
  df <- df[is.finite(df[[pcol]]) & is.finite(df$log2FoldChange), ]
  df$Signif <- "non_sig"
  df$Signif[df[[pcol]] < pthresh & df$log2FoldChange > lfc_thresh] <- "up"
  df$Signif[df[[pcol]] < pthresh & df$log2FoldChange < -lfc_thresh] <- "down"
  df$Signif <- factor(df$Signif, levels = c("down", "non_sig", "up"))

  df$neglog10p <- -log10(df[[pcol]])
  top_up <- df[df$Signif == "up", ]
  top_up <- top_up[order(-top_up$log2FoldChange), ]
  top_up <- head(top_up, 5)

  top_down <- df[df$Signif == "down", ]
  top_down <- top_down[order(top_down$log2FoldChange), ]
  top_down <- head(top_down, 5)

  volcano <- ggplot(df, aes(x = log2FoldChange, y = neglog10p, color = Signif)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_manual(
      name = "Significant genes",
      values = c(down = "blue", non_sig = "black", up = "red"),
      labels = c(down = "Down", non_sig = "Non-sig", up = "Up")
    ) +
    geom_label_repel(
      data = rbind(top_up, top_down),
      aes(label = GeneSymbol),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE,
      label.size = 0.2,
      label.r = grid::unit(0.1, "lines")
    ) +
    theme_minimal() +
    labs(
      title = sprintf("Volcano: %s", label),
      x = "log2FoldChange",
      y = sprintf("-log10(%s)", pcol)
    )

  ggsave(
    file.path(out_dir, "volcano", sprintf("%s.png", label)),
    volcano,
    width = 7,
    height = 5,
    dpi = 300
  )
}

plot_ma <- function(res_df, label, pcol, out_dir) {
  df <- res_df
  df <- df[is.finite(df[[pcol]]) & is.finite(df$log2FoldChange) & is.finite(df$baseMean), ]
  df$Signif <- "non_sig"
  df$Signif[df[[pcol]] < pthresh & df$log2FoldChange > lfc_thresh] <- "up"
  df$Signif[df[[pcol]] < pthresh & df$log2FoldChange < -lfc_thresh] <- "down"
  df$Signif <- factor(df$Signif, levels = c("down", "non_sig", "up"))

  df$log10baseMean <- log10(df$baseMean + 1)
  top_up <- df[df$Signif == "up", ]
  top_up <- top_up[order(-top_up$log2FoldChange), ]
  top_up <- head(top_up, 5)

  top_down <- df[df$Signif == "down", ]
  top_down <- top_down[order(top_down$log2FoldChange), ]
  top_down <- head(top_down, 5)

  ma <- ggplot(df, aes(x = log10baseMean, y = log2FoldChange, color = Signif)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_manual(
      name = "Significant genes",
      values = c(down = "blue", non_sig = "black", up = "red"),
      labels = c(down = "Down", non_sig = "Non-sig", up = "Up")
    ) +
    geom_label_repel(
      data = rbind(top_up, top_down),
      aes(label = GeneSymbol),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE,
      label.size = 0.2,
      label.r = grid::unit(0.1, "lines")
    ) +
    theme_minimal() +
    labs(title = sprintf("MA: %s", label), x = "log10(baseMean + 1)", y = "log2FoldChange")

  ggsave(
    file.path(out_dir, "ma", sprintf("%s.png", label)),
    ma,
    width = 7,
    height = 5,
    dpi = 300
  )
}

plot_volcano(res_list$gut_duct, "de_gut_duct", pcol, output_plots)
plot_volcano(res_list$duct_node, "de_duct_node", pcol, output_plots)
plot_volcano(res_list$node_gut, "de_node_gut", pcol, output_plots)

plot_ma(res_list$gut_duct, "de_gut_duct", pcol, output_plots)
plot_ma(res_list$duct_node, "de_duct_node", pcol, output_plots)
plot_ma(res_list$node_gut, "de_node_gut", pcol, output_plots)

# Boxplots for top 10 DE genes (scaled values)
sig_all <- Reduce(`|`, lapply(res_list, function(x) x[[pcol]] < pthresh & abs(x$log2FoldChange) > lfc_thresh))
sig_gene_symbols <- res_list$gut_duct$GeneSymbol[sig_all]
sig_gene_symbols <- unique(sig_gene_symbols)

min_p_sig <- min_p_all
if (length(sig_gene_symbols) > 0) {
  min_p_sig <- min_p_sig[names(min_p_sig) %in% sig_gene_symbols]
}
min_p_sig <- min_p_sig[order(min_p_sig, na.last = TRUE)]

top10 <- names(min_p_sig)[seq_len(min(10, length(min_p_sig)))]

box_mat <- scaled_mat[top10, , drop = FALSE]
box_df <- data.frame(
  Gene = rep(rownames(box_mat), times = ncol(box_mat)),
  Sample = rep(colnames(box_mat), each = nrow(box_mat)),
  Value = as.vector(box_mat),
  Group = rep(col_data$group, each = nrow(box_mat))
)

box_plot <- ggplot(box_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ Gene, ncol = 5, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Top 10 DE genes (scaled)", x = "Group", y = "Scaled expression")

ggsave(file.path(output_plots, "boxplots_top10.png"), box_plot, width = 8, height = 12, dpi = 300)

combined <- do.call(
  rbind,
  lapply(names(res_list), function(name) {
    df <- res_list[[name]]
    df$Contrast <- name
    df
  })
)

combined <- combined[is.finite(combined[[pcol]]) & is.finite(combined$log2FoldChange), ]
combined <- combined[combined[[pcol]] < pthresh & abs(combined$log2FoldChange) > lfc_thresh, ]

combined_up <- combined[order(-combined$log2FoldChange), ]
combined_down <- combined[order(combined$log2FoldChange), ]

combined_up <- head(combined_up, 10)
combined_down <- head(combined_down, 10)

cols <- c("GeneSymbol", "log2FoldChange", "pvalue", "padj", "Contrast")
write.table(
  combined_up[, cols],
  file.path(output_plots, "top10_up.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_down[, cols],
  file.path(output_plots, "top10_down.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
