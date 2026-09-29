# tEPOR gene editing study: erythroid RNA-seq analysis
# Set RNASEQ_COUNTS and RNASEQ_OUTPUT or edit Settings below.
# Input must contain Geneid and sample count columns only (no featureCounts annotation columns).
# Original models, filters, markers, and plot parameters are retained.
############################################################
# tEPOR erythroid maturation RNA-seq analysis — v2
#
# Input:
#   gene_counts_featureCounts_d0_d7_d11_d14.txt
#
# Experimental structure:
#   _1/_2/_3 = replicates from the SAME donor
#
# Main comparisons:
#   D7  : tEPOR vs MOCK
#   D11 : tEPOR vs MOCK
#   D14 : tEPOR vs MOCK
#
# D0:
#   MOCK only -> visualization baseline only
#
# D14_diluted_tEPOR:
#   excluded from main analysis
#
# Outputs:
#   - normalized counts
#   - marker trajectory plots
#   - per-day DESeq2
#   - day x treatment interaction
#   - erythroid maturation heatmap
############################################################


############################################################
# 0. Packages
############################################################

# CRAN packages
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "pheatmap"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "DESeq2",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})


############################################################
# 1. Settings
############################################################

input_file <- Sys.getenv("RNASEQ_COUNTS", "gene_counts_featureCounts_d0_d7_d11_d14.txt")

output_dir <- Sys.getenv("RNASEQ_OUTPUT", "tEPOR_erythroid_RNA_analysis_v2")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


############################################################
# 2. Marker gene sets
############################################################

# Direct RNA equivalents of flow markers
flow_marker_genes <- c(
  "GYPA",       # GPA / CD235a
  "TFRC",       # CD71
  "CD36",       # CD36
  "ITGA4"       # CD49d / integrin alpha4
)

# Terminal erythroid maturation genes
terminal_marker_genes <- c(
  "SLC4A1",     # Band3 / CD233
  "AHSP",
  "ANK1",
  "SPTB"
)

# All markers
marker_genes <- unique(c(
  flow_marker_genes,
  terminal_marker_genes
))


############################################################
# 3. Read raw featureCounts file
############################################################

counts_raw <- read.delim(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("INPUT FILE\n")
cat("========================================\n")

cat("Dimensions:\n")
print(dim(counts_raw))

cat("\nColumn names:\n")
print(colnames(counts_raw))


############################################################
# 4. Remove D14 diluted tEPOR
############################################################

diluted_columns <- grepl(
  "diluted_tEPOR",
  colnames(counts_raw),
  ignore.case = TRUE
)

if (any(diluted_columns)) {

  cat("\nRemoving diluted tEPOR samples:\n")

  print(
    colnames(counts_raw)[diluted_columns]
  )

  counts_raw <- counts_raw[
    ,
    !diluted_columns,
    drop = FALSE
  ]
}


############################################################
# 5. Detect/remove genes with NA raw counts
############################################################

count_cols <- setdiff(
  colnames(counts_raw),
  "Geneid"
)

# Check that all count columns are numeric/integer
cat("\nCount column classes:\n")
print(
  table(
    sapply(
      counts_raw[, count_cols, drop = FALSE],
      class
    )
  )
)

# Find rows containing at least one NA
bad_rows <- which(
  apply(
    counts_raw[, count_cols, drop = FALSE],
    1,
    anyNA
  )
)

if (length(bad_rows) > 0) {

  cat("\n========================================\n")
  cat("WARNING: genes containing NA counts\n")
  cat("========================================\n")

  cat(
    "Number of affected genes:",
    length(bad_rows),
    "\n"
  )

  print(
    counts_raw[
      bad_rows,
      "Geneid",
      drop = FALSE
    ]
  )

  # Save them before removal
  write.csv(
    counts_raw[bad_rows, ],
    file.path(
      output_dir,
      "removed_genes_with_NA_counts.csv"
    ),
    row.names = FALSE
  )

  # Remove
  counts_raw <- counts_raw[
    -bad_rows,
    ,
    drop = FALSE
  ]
}

# Final NA check
count_cols <- setdiff(
  colnames(counts_raw),
  "Geneid"
)

if (
  anyNA(
    counts_raw[, count_cols, drop = FALSE]
  )
) {
  stop(
    "ERROR: NA values remain in raw counts."
  )
}

cat(
  "\nNA check passed.\nRemaining genes:",
  nrow(counts_raw),
  "\n"
)


############################################################
# 6. Prepare count matrix
############################################################

gene_ids <- counts_raw$Geneid

# Remove Ensembl version suffix if present
# ENSG00000123456.12 -> ENSG00000123456
gene_ids <- sub(
  "\\..*$",
  "",
  gene_ids
)

count_df <- counts_raw[
  ,
  count_cols,
  drop = FALSE
]

# Safe numeric conversion
count_df[] <- lapply(
  count_df,
  function(x) {
    as.numeric(as.character(x))
  }
)

# Check again after conversion
if (anyNA(count_df)) {
  stop(
    "ERROR: numeric conversion generated NA values."
  )
}

# DESeq2 requires integer counts
count_matrix <- as.matrix(count_df)

storage.mode(count_matrix) <- "integer"

rownames(count_matrix) <- gene_ids


############################################################
# 7. Handle duplicated Ensembl IDs
############################################################

duplicated_ids <- duplicated(
  rownames(count_matrix)
)

if (any(duplicated_ids)) {

  cat(
    "\nRemoving",
    sum(duplicated_ids),
    "duplicated Ensembl IDs.\n"
  )

  count_matrix <- count_matrix[
    !duplicated_ids,
    ,
    drop = FALSE
  ]
}


############################################################
# 8. Count matrix QC
############################################################

cat("\n========================================\n")
cat("COUNT MATRIX QC\n")
cat("========================================\n")

cat(
  "Genes:",
  nrow(count_matrix),
  "\n"
)

cat(
  "Samples:",
  ncol(count_matrix),
  "\n"
)

cat(
  "Total NA:",
  sum(is.na(count_matrix)),
  "\n"
)

cat(
  "Minimum count:",
  min(count_matrix),
  "\n"
)

cat(
  "Maximum count:",
  max(count_matrix),
  "\n"
)

cat("\nLibrary sizes:\n")
print(
  colSums(count_matrix)
)

if (any(count_matrix < 0)) {
  stop(
    "ERROR: negative counts detected."
  )
}


############################################################
# 9. Build metadata
############################################################

sample_names <- colnames(
  count_matrix
)

metadata <- data.frame(
  sample = sample_names,
  stringsAsFactors = FALSE
)

metadata$day <- dplyr::case_when(
  grepl(
    "^D0_",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "D0",

  grepl(
    "^D7_",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "D7",

  grepl(
    "^D11_",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "D11",

  grepl(
    "^D14_",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "D14",

  TRUE ~ NA_character_
)

metadata$treatment <- dplyr::case_when(

  grepl(
    "MOCK",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "MOCK",

  grepl(
    "tEPOR",
    metadata$sample,
    ignore.case = TRUE
  ) ~ "tEPOR",

  TRUE ~ NA_character_
)

metadata$replicate <- sub(
  ".*_",
  "",
  metadata$sample
)

metadata$day <- factor(
  metadata$day,
  levels = c(
    "D0",
    "D7",
    "D11",
    "D14"
  )
)

metadata$treatment <- factor(
  metadata$treatment,
  levels = c(
    "MOCK",
    "tEPOR"
  )
)

rownames(metadata) <- metadata$sample


############################################################
# 10. Metadata QC
############################################################

if (any(is.na(metadata$day))) {

  print(
    metadata[
      is.na(metadata$day),
      ,
      drop = FALSE
    ]
  )

  stop(
    "ERROR: some samples could not be assigned to day."
  )
}

if (any(is.na(metadata$treatment))) {

  print(
    metadata[
      is.na(metadata$treatment),
      ,
      drop = FALSE
    ]
  )

  stop(
    "ERROR: some samples could not be assigned to treatment."
  )
}

if (
  !all(
    rownames(metadata) ==
      colnames(count_matrix)
  )
) {
  stop(
    "ERROR: metadata/sample order mismatch."
  )
}

cat("\n========================================\n")
cat("SAMPLE METADATA\n")
cat("========================================\n")

print(metadata)

cat("\nSample numbers:\n")

print(
  table(
    metadata$day,
    metadata$treatment
  )
)

write.csv(
  metadata,
  file.path(
    output_dir,
    "sample_metadata.csv"
  ),
  row.names = FALSE
)


############################################################
# 11. Ensembl -> gene symbol
############################################################

gene_annotation <- data.frame(
  ENSEMBL = rownames(count_matrix),
  stringsAsFactors = FALSE
)

gene_annotation$SYMBOL <-
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = gene_annotation$ENSEMBL,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

cat("\n========================================\n")
cat("ANNOTATION\n")
cat("========================================\n")

cat(
  "Genes with symbol:",
  sum(!is.na(gene_annotation$SYMBOL)),
  "\n"
)

cat(
  "Genes without symbol:",
  sum(is.na(gene_annotation$SYMBOL)),
  "\n"
)

cat("\nAnnotation example:\n")
print(
  head(gene_annotation)
)

write.csv(
  gene_annotation,
  file.path(
    output_dir,
    "Ensembl_gene_annotation.csv"
  ),
  row.names = FALSE
)


############################################################
# 12. DESeq2 normalization
############################################################

dds_norm <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ 1
)

# Basic low-count filtering
keep_gene <- rowSums(
  counts(dds_norm)
) >= 10

cat(
  "\nGenes before low-count filter:",
  nrow(dds_norm),
  "\n"
)

dds_norm <- dds_norm[
  keep_gene,
]

cat(
  "Genes after low-count filter:",
  nrow(dds_norm),
  "\n"
)

dds_norm <- estimateSizeFactors(
  dds_norm
)

normalized_counts <- counts(
  dds_norm,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  file.path(
    output_dir,
    "DESeq2_normalized_counts.csv"
  )
)


############################################################
# 13. Add gene symbols
############################################################

norm_df <- as.data.frame(
  normalized_counts
)

norm_df <- tibble::rownames_to_column(
  norm_df,
  "ENSEMBL"
)

norm_df <- dplyr::left_join(
  norm_df,
  gene_annotation,
  by = "ENSEMBL"
)

write.csv(
  norm_df,
  file.path(
    output_dir,
    "DESeq2_normalized_counts_with_symbols.csv"
  ),
  row.names = FALSE
)


############################################################
# 14. Check requested markers
############################################################

available_markers <- marker_genes[
  marker_genes %in%
    norm_df$SYMBOL
]

missing_markers <- setdiff(
  marker_genes,
  available_markers
)

cat("\n========================================\n")
cat("ERYTHROID MARKERS\n")
cat("========================================\n")

cat("\nAvailable:\n")
print(available_markers)

cat("\nMissing:\n")
print(missing_markers)


############################################################
# 15. Extract marker expression
############################################################

marker_df <- dplyr::filter(
  norm_df,
  SYMBOL %in% available_markers
)

# Calculate mean expression across all samples.
# Explicit dplyr::select prevents AnnotationDbi conflict.
marker_df$mean_expression <- rowMeans(
  dplyr::select(
    marker_df,
    dplyr::all_of(sample_names)
  ),
  na.rm = TRUE
)

# In case more than one Ensembl ID maps to same symbol,
# retain the most highly expressed entry.
marker_df <- marker_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(
    order_by = mean_expression,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()


############################################################
# 16. Convert marker table to long format
############################################################

marker_long <- marker_df %>%
  dplyr::select(
    SYMBOL,
    dplyr::all_of(sample_names)
  ) %>%
  tidyr::pivot_longer(
    cols = -SYMBOL,
    names_to = "sample",
    values_to = "normalized_count"
  )

metadata_join <- metadata %>%
  tibble::rownames_to_column(
    "metadata_row"
  )

marker_long <- dplyr::left_join(
  marker_long,
  metadata_join,
  by = "sample"
)

marker_long <- marker_long %>%
  dplyr::mutate(
    log2_expression =
      log2(normalized_count + 1)
  )

write.csv(
  marker_long,
  file.path(
    output_dir,
    "erythroid_marker_normalized_expression.csv"
  ),
  row.names = FALSE
)


############################################################
# 17. Flow-corresponding RNA marker plot
############################################################

flow_long <- dplyr::filter(
  marker_long,
  SYMBOL %in% flow_marker_genes
)

p_flow <- ggplot(
  flow_long,
  aes(
    x = day,
    y = log2_expression,
    group = treatment
  )
) +

  geom_point(
    aes(shape = treatment),
    size = 2.5,
    position =
      position_jitter(
        width = 0.06,
        height = 0
      )
  ) +

  stat_summary(
    aes(linetype = treatment),
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +

  stat_summary(
    aes(shape = treatment),
    fun = mean,
    geom = "point",
    size = 3
  ) +

  facet_wrap(
    ~ SYMBOL,
    scales = "free_y",
    ncol = 2
  ) +

  theme_classic(
    base_size = 13
  ) +

  labs(
    x = "Day of erythroid differentiation",
    y = "log2(DESeq2 normalized count + 1)",
    shape = "Treatment",
    linetype = "Treatment",
    title =
      "RNA expression of flow-corresponding erythroid markers"
  ) +

  theme(
    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),

    legend.position =
      "right"
  )

ggsave(
  file.path(
    output_dir,
    "Flow_marker_RNA_trajectory.pdf"
  ),
  p_flow,
  width = 8,
  height = 6
)

ggsave(
  file.path(
    output_dir,
    "Flow_marker_RNA_trajectory.png"
  ),
  p_flow,
  width = 8,
  height = 6,
  dpi = 300
)


############################################################
# 18. Marker mean +/- SEM
############################################################

marker_summary <- marker_long %>%
  dplyr::group_by(
    SYMBOL,
    day,
    treatment
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),

    mean =
      mean(
        log2_expression,
        na.rm = TRUE
      ),

    sd =
      sd(
        log2_expression,
        na.rm = TRUE
      ),

    sem =
      sd / sqrt(n),

    .groups = "drop"
  )

write.csv(
  marker_summary,
  file.path(
    output_dir,
    "erythroid_marker_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 19. All marker trajectories
############################################################

p_all <- ggplot(
  marker_summary,
  aes(
    x = day,
    y = mean,
    group = treatment,
    linetype = treatment,
    shape = treatment
  )
) +

  geom_line(
    linewidth = 0.8
  ) +

  geom_point(
    size = 2.5
  ) +

  geom_errorbar(
    aes(
      ymin = mean - sem,
      ymax = mean + sem
    ),
    width = 0.1,
    linewidth = 0.5
  ) +

  facet_wrap(
    ~ SYMBOL,
    scales = "free_y",
    ncol = 4
  ) +

  theme_classic(
    base_size = 11
  ) +

  labs(
    x = "Day of erythroid differentiation",
    y = "log2(DESeq2 normalized count + 1)",
    linetype = "Treatment",
    shape = "Treatment",
    title =
      "Erythroid maturation RNA expression kinetics"
  ) +

  theme(
    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "All_erythroid_markers_trajectory.pdf"
  ),
  p_all,
  width = 12,
  height = 6
)

ggsave(
  file.path(
    output_dir,
    "All_erythroid_markers_trajectory.png"
  ),
  p_all,
  width = 12,
  height = 6,
  dpi = 300
)


############################################################
# 20. Per-day differential expression
#     tEPOR vs MOCK
#     D7 / D11 / D14
############################################################

comparison_days <- c("D7", "D11", "D14")

all_DE_results <- list()

for (current_day in comparison_days) {

  cat(
    "\n========================================\n",
    current_day, ": tEPOR vs MOCK\n",
    "========================================\n",
    sep = ""
  )

  # Select samples for current day
  samples_day <- rownames(metadata)[
    metadata$day == current_day
  ]

  counts_day <- count_matrix[
    ,
    samples_day,
    drop = FALSE
  ]

  metadata_day <- metadata[
    samples_day,
    ,
    drop = FALSE
  ]

  # Remove unused factor levels
  metadata_day$treatment <- droplevels(
    metadata_day$treatment
  )

  # Set MOCK as reference
  metadata_day$treatment <- relevel(
    metadata_day$treatment,
    ref = "MOCK"
  )

  # Check sample numbers
  cat("\nSamples:\n")
  print(metadata_day[, c("sample", "day", "treatment", "replicate")])

  cat("\nTreatment counts:\n")
  print(table(metadata_day$treatment))

  # Create DESeq2 object
  dds_day <- DESeqDataSetFromMatrix(
    countData = counts_day,
    colData = metadata_day,
    design = ~ treatment
  )

  # Remove very low-count genes
  keep <- rowSums(counts(dds_day)) >= 10
  dds_day <- dds_day[keep, ]

  # Run DESeq2
  dds_day <- DESeq(
    dds_day,
    quiet = TRUE
  )

  # tEPOR vs MOCK
  res_day <- results(
    dds_day,
    contrast = c(
      "treatment",
      "tEPOR",
      "MOCK"
    ),
    alpha = 0.05
  )

  # Convert to data frame
  res_df <- as.data.frame(res_day)

  res_df <- tibble::rownames_to_column(
    res_df,
    "ENSEMBL"
  )

  # Add gene symbols
  res_df <- dplyr::left_join(
    res_df,
    gene_annotation,
    by = "ENSEMBL"
  )

  # Add day
  res_df$day <- current_day

  # Sort by adjusted p value
  res_df <- dplyr::arrange(
    res_df,
    padj
  )

  # Store result
  all_DE_results[[current_day]] <- res_df

  # Save individual day result
  write.csv(
    res_df,
    file.path(
      output_dir,
      paste0(
        current_day,
        "_tEPOR_vs_MOCK_DESeq2.csv"
      )
    ),
    row.names = FALSE
  )

}

cat("\nPer-day DESeq2 analysis completed.\n")


############################################################
# 21. Combine per-day DE results
############################################################

DE_combined <- dplyr::bind_rows(
  all_DE_results
)

write.csv(
  DE_combined,
  file.path(
    output_dir,
    "All_days_tEPOR_vs_MOCK_DESeq2.csv"
  ),
  row.names = FALSE
)


############################################################
# 22. Marker-specific per-day statistics
############################################################

marker_stats <- DE_combined %>%

  dplyr::filter(
    SYMBOL %in% marker_genes
  ) %>%

  dplyr::select(
    day,
    SYMBOL,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  ) %>%

  dplyr::arrange(
    SYMBOL,
    day
  )

write.csv(
  marker_stats,
  file.path(
    output_dir,
    "Erythroid_marker_tEPOR_vs_MOCK_statistics.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("ERYTHROID MARKER STATISTICS\n")
cat("========================================\n")

print(marker_stats)


############################################################
# 23. Day x treatment interaction
#
# IMPORTANT:
# D0 excluded because there is no D0 tEPOR.
#
# Full:
# ~ day + treatment + day:treatment
#
# Reduced:
# ~ day + treatment
#
# LRT asks whether temporal trajectory differs
# between tEPOR and MOCK.
############################################################

interaction_samples <-
  rownames(metadata)[
    metadata$day %in%
      c(
        "D7",
        "D11",
        "D14"
      )
  ]

counts_interaction <-
  count_matrix[
    ,
    interaction_samples,
    drop = FALSE
  ]

metadata_interaction <-
  metadata[
    interaction_samples,
    ,
    drop = FALSE
  ]

metadata_interaction$day <-
  droplevels(
    metadata_interaction$day
  )

metadata_interaction$treatment <-
  droplevels(
    metadata_interaction$treatment
  )

metadata_interaction$day <-
  relevel(
    metadata_interaction$day,
    ref = "D7"
  )

metadata_interaction$treatment <-
  relevel(
    metadata_interaction$treatment,
    ref = "MOCK"
  )

dds_interaction <-
  DESeqDataSetFromMatrix(
    countData =
      counts_interaction,

    colData =
      metadata_interaction,

    design =
      ~ day +
        treatment +
        day:treatment
  )

keep_interaction <-
  rowSums(
    counts(
      dds_interaction
    )
  ) >= 10

dds_interaction <-
  dds_interaction[
    keep_interaction,
  ]

dds_interaction <- DESeq(
  dds_interaction,
  test = "LRT",
  reduced =
    ~ day + treatment,
  quiet = TRUE
)

interaction_results <-
  results(
    dds_interaction,
    alpha = 0.05
  )

interaction_df <-
  as.data.frame(
    interaction_results
  )

interaction_df <-
  tibble::rownames_to_column(
    interaction_df,
    "ENSEMBL"
  )

interaction_df <-
  dplyr::left_join(
    interaction_df,
    gene_annotation,
    by = "ENSEMBL"
  )

interaction_df <-
  dplyr::arrange(
    interaction_df,
    padj
  )

write.csv(
  interaction_df,
  file.path(
    output_dir,
    "Day_by_tEPOR_interaction_DESeq2_LRT.csv"
  ),
  row.names = FALSE
)


############################################################
# 24. Marker interaction statistics
############################################################

interaction_marker_results <-
  interaction_df %>%

  dplyr::filter(
    SYMBOL %in%
      marker_genes
  ) %>%

  dplyr::select(
    ENSEMBL,
    SYMBOL,
    baseMean,
    stat,
    pvalue,
    padj
  ) %>%

  dplyr::arrange(
    padj
  )

write.csv(
  interaction_marker_results,
  file.path(
    output_dir,
    "Erythroid_marker_day_by_tEPOR_interaction.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("DAY x tEPOR INTERACTION — MARKERS\n")
cat("========================================\n")

print(
  interaction_marker_results
)


############################################################
# 25. VST for heatmap
############################################################

vst_object <- vst(
  dds_norm,
  blind = TRUE
)

vst_counts <- assay(
  vst_object
)

vst_df <- as.data.frame(
  vst_counts
)

vst_df <-
  tibble::rownames_to_column(
    vst_df,
    "ENSEMBL"
  )

vst_df <-
  dplyr::left_join(
    vst_df,
    gene_annotation,
    by = "ENSEMBL"
  )


############################################################
# 26. Prepare erythroid heatmap
############################################################

heatmap_df <-
  dplyr::filter(
    vst_df,
    SYMBOL %in%
      marker_genes
  )

# Calculate variance across samples
heatmap_df$gene_variance <-
  apply(
    heatmap_df[
      ,
      sample_names,
      drop = FALSE
    ],
    1,
    var
  )

# Resolve duplicated symbols
heatmap_df <-
  heatmap_df %>%

  dplyr::group_by(
    SYMBOL
  ) %>%

  dplyr::slice_max(
    order_by =
      gene_variance,

    n = 1,

    with_ties =
      FALSE
  ) %>%

  dplyr::ungroup()


############################################################
# 27. Create heatmap matrix
############################################################

heatmap_matrix <-
  as.matrix(
    heatmap_df[
      ,
      sample_names,
      drop = FALSE
    ]
  )

rownames(
  heatmap_matrix
) <- heatmap_df$SYMBOL


############################################################
# 28. Heatmap sample annotation
############################################################

annotation_col <-
  data.frame(
    Day =
      metadata$day,

    Treatment =
      metadata$treatment
  )

rownames(
  annotation_col
) <- rownames(
  metadata
)


############################################################
# 29. Order heatmap samples
############################################################

sample_order <-
  order(
    metadata$day,
    metadata$treatment,
    metadata$replicate
  )

heatmap_matrix <-
  heatmap_matrix[
    ,
    sample_order,
    drop = FALSE
  ]

annotation_col <-
  annotation_col[
    sample_order,
    ,
    drop = FALSE
  ]


############################################################
# 30. Heatmap PDF
############################################################

pdf(
  file.path(
    output_dir,
    "Erythroid_maturation_heatmap.pdf"
  ),
  width = 11,
  height = 8
)

pheatmap::pheatmap(
  heatmap_matrix,

  scale = "row",

  cluster_rows = TRUE,

  cluster_cols = FALSE,

  annotation_col =
    annotation_col,

  show_colnames = TRUE,

  fontsize_row = 10,

  fontsize_col = 7,

  border_color = NA,

  main =
    "Erythroid maturation transcriptional program"
)

dev.off()


############################################################
# 31. Heatmap PNG
############################################################

png(
  file.path(
    output_dir,
    "Erythroid_maturation_heatmap.png"
  ),
  width = 3300,
  height = 1800,
  res = 300
)

pheatmap::pheatmap(
  heatmap_matrix,

  scale = "row",

  cluster_rows = TRUE,

  cluster_cols = FALSE,

  annotation_col =
    annotation_col,

  show_colnames = TRUE,

  fontsize_row = 10,

  fontsize_col = 7,

  border_color = NA,

  main =
    "Erythroid maturation transcriptional program"
)

dev.off()


############################################################
# 32. Flow/RNA correspondence statistics
############################################################

flow_marker_stats <-
  marker_stats %>%

  dplyr::filter(
    SYMBOL %in%
      c(
        "GYPA",
        "TFRC",
        "CD36",
        "ITGA4",
        "SLC4A1"
      )
  ) %>%

  dplyr::mutate(

    Flow_marker =
      dplyr::case_when(

        SYMBOL == "GYPA" ~
          "GPA / CD235a",

        SYMBOL == "TFRC" ~
          "CD71",

        SYMBOL == "CD36" ~
          "CD36",

        SYMBOL == "ITGA4" ~
          "CD49d / integrin alpha4",

        SYMBOL == "SLC4A1" ~
          "Band3 / CD233",

        TRUE ~ SYMBOL
      )
  ) %>%

  dplyr::select(
    day,
    Flow_marker,
    SYMBOL,
    log2FoldChange,
    pvalue,
    padj
  )

write.csv(
  flow_marker_stats,
  file.path(
    output_dir,
    "Flow_vs_RNA_marker_statistics.csv"
  ),
  row.names = FALSE
)


############################################################
# 33. Save session info
############################################################

sink(
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

sessionInfo()

sink()


############################################################
# 34. FINISHED
############################################################

cat("\n\n")
cat("========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nResults saved in:\n",
  output_dir,
  "\n",
  sep = ""
)

cat("\nKey outputs:\n")

cat(
  "1. Flow_marker_RNA_trajectory.pdf\n"
)

cat(
  "2. All_erythroid_markers_trajectory.pdf\n"
)

cat(
  "3. Erythroid_maturation_heatmap.pdf\n"
)

cat(
  "4. Erythroid_marker_tEPOR_vs_MOCK_statistics.csv\n"
)

cat(
  "5. Erythroid_marker_day_by_tEPOR_interaction.csv\n"
)

cat(
  "6. Flow_vs_RNA_marker_statistics.csv\n"
)

cat(
  "7. Day_by_tEPOR_interaction_DESeq2_LRT.csv\n"
)

cat(
  "8. All_days_tEPOR_vs_MOCK_DESeq2.csv\n"
)

cat("\nDone.\n")
