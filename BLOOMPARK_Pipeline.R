# =============================================================================
# BLOOMPARK
# BLOOd Multi-omics for PARKinson Disease's early detection
#
# Comprehensive Multi-Transcriptomic Analysis of Peripheral Blood to Identify
# Key Genes in Parkinson's Disease
#
# Analysis overview:
#   01. Count matrix assembly (RNA-seq & miRNA-seq)
#   02. Metadata preparation and blood chemistry integration
#   03. Data cleaning and sample harmonization
#   04. Gene annotation and symbol mapping
#   05. Differential expression analysis (DEGs & DEmiRNAs per stage)
#   06. miRNA target retrieval and miRNA–mRNA correlation analysis
# =============================================================================


# ── Dependencies ──────────────────────────────────────────────────────────────

library(readr)
library(data.table)
library(readxl)
library(edgeR)
library(dplyr)
library(tidyr)
library(R.utils)
library(ggplot2)
library(AnnotationHub)
library(ensembldb)
library(AnnotationDbi)
library(HGNChelper)
library(MIRit)
library(doParallel)
library(BiocParallel)
library(patchwork)
library(UpSetR)


# ── Working directory & shared output log ─────────────────────────────────────

setwd("/home/dangelo/BLOOMPARK/")

STATS_FILE <- "General_statistics.txt"   # Flat key–value log appended throughout


# =============================================================================
# 01. COUNT MATRIX ASSEMBLY
# =============================================================================

# ── RNA-seq (featureCounts) ───────────────────────────────────────────────────

rnaseq_counts_dir <- "/PPMI_RNAseq_IR3_Analysis/counts"

fc_files <- list.files(
  rnaseq_counts_dir,
  pattern  = "featureCounts",
  full.names = TRUE
)

# Build a Unique_ID from filename components: <PATNO>_<EVENT_ID>
make_unique_id_rnaseq <- function(filepath) {
  parts <- strsplit(basename(filepath), "\\.")[[1]]
  paste0(parts[2], "_", parts[3])
}

file_unique_ids        <- vapply(fc_files, make_unique_id_rnaseq, character(1))
names(fc_files)        <- file_unique_ids

# Gene IDs from the first file (column 1); counts from column 7 (featureCounts layout)
gene_ids <- read.table(
  fc_files[1],
  header       = TRUE,
  sep          = "\t",
  comment.char = "#",
  quote        = ""
)[[1]]

read_counts_column <- function(filepath) {
  fread(filepath, select = 7, skip = "#", colClasses = "integer")[[1]]
}

rnaseq_counts_raw <- do.call(cbind, lapply(fc_files, read_counts_column))
rnaseq_counts_raw <- cbind(Geneid = gene_ids, rnaseq_counts_raw)

cat("RNAseq_samples", ncol(rnaseq_counts_raw) - 1, "\n", file = STATS_FILE, append = TRUE)
cat("RNAseq_genes",   nrow(rnaseq_counts_raw),       "\n", file = STATS_FILE, append = TRUE)


# ── miRNA-seq (pre-assembled matrix) ─────────────────────────────────────────

mirna_counts_dir <- "/PPMI_sncRNAcounts/counts"

mirna_counts_raw <- fread(file.path(mirna_counts_dir, "mirna_quantification_matrix_raw.csv"))

# Rename sample columns to <PATNO>_<EVENT_ID> to match RNA-seq convention
make_unique_id_mirna <- function(colname) {
  parts <- strsplit(colname, "\\.")[[1]]
  paste0(parts[2], "_", parts[3])
}

sample_cols_mirna        <- colnames(mirna_counts_raw)[-1]
unique_ids_mirna         <- vapply(sample_cols_mirna, make_unique_id_mirna, character(1))
colnames(mirna_counts_raw)[-1] <- unique_ids_mirna

mirna_counts_raw <- as.data.frame(mirna_counts_raw)

cat("miRNA_samples", ncol(mirna_counts_raw) - 1, "\n", file = STATS_FILE, append = TRUE)
cat("miRNA_genes",   nrow(mirna_counts_raw),       "\n", file = STATS_FILE, append = TRUE)


# =============================================================================
# 02. METADATA PREPARATION
# =============================================================================

# ── Core clinical metadata ────────────────────────────────────────────────────

metadata_raw <- readxl::read_excel(
  "PPMI_Curated_Data_Cut_Public_20250714.xlsx"
)

metadata_raw <- metadata_raw %>%
  mutate(Unique_ID = paste0(PATNO, "_", EVENT_ID)) %>%
  as.data.frame()

rownames(metadata_raw) <- metadata_raw$Unique_ID

cat("Metadata_samples", nrow(metadata_raw), "\n", file = STATS_FILE, append = TRUE)


# ── NSD stage filtering & collapsing ─────────────────────────────────────────
# Keep records that are NSD-positive, plus HC (COHORT 2) and Prodromal (COHORT 4)
# where NSD_STAGE is NA.

metadata_nsd <- dplyr::filter(
  metadata_raw,
  NSD_STAGE != "Not NSD" | (is.na(NSD_STAGE) & COHORT %in% c(2, 4))
)

# Map granular NSD stages to four analysis-ready severity categories
metadata_nsd <- metadata_nsd %>%
  mutate(New_Stage = case_when(
    COHORT == 2                          ~ "Healthy",
    COHORT == 4                          ~ "Pre-clinical",
    NSD_STAGE %in% c("1a", "1b", "2a") ~ "Pre-clinical",
    NSD_STAGE == "2b"                   ~ "Mild",
    NSD_STAGE == "3"                    ~ "Moderate",
    NSD_STAGE %in% c("4", "5", "6")    ~ "Severe",
    TRUE                                 ~ NA_character_
  ))


# ── Blood chemistry integration ───────────────────────────────────────────────
# Only rows with GI/L units are used (cellular blood counts: WBC, Platelets, etc.)

blood_chem_raw <- read.delim(
  "Blood_Chemistry___Hematology-Archived_05Sep2025.csv",
  sep    = ",",
  header = TRUE
)

blood_chem_raw$Unique_ID <- paste0(blood_chem_raw$PATNO, "_", blood_chem_raw$EVENT_ID)

blood_chem_gil <- blood_chem_raw[blood_chem_raw$LSIUNIT == "GI/L",
                                   c("Unique_ID", "PATNO", "EVENT_ID", "LTSTNAME", "LSIRES")]

# Reshape from long to wide: one column per analyte
blood_chem_wide <- pivot_wider(
  blood_chem_gil,
  names_from  = LTSTNAME,
  values_from = LSIRES
)

# Join the four key analytes to the NSD-filtered metadata
metadata_with_chem <- left_join(
  metadata_nsd,
  blood_chem_wide[, c("Unique_ID", "WBC", "Platelets", "Neutrophils", "Lymphocytes")],
  by = "Unique_ID"
)

# Coerce list-columns produced by pivot_wider to plain numeric vectors
metadata_with_chem <- metadata_with_chem %>%
  mutate(across(
    c(WBC, Platelets, Neutrophils, Lymphocytes),
    ~ as.numeric(sapply(.x, function(v) if (is.null(v) || length(v) == 0) NA else v[[1]]))
  ))


# ── Missing value imputation (hierarchical: patient → stage → global) ─────────
# Strategy: within-patient mean first, then within-stage mean, then global mean.

impute_hierarchical <- function(df, vars) {
  # Pass 1: within-patient mean
  df <- df %>%
    group_by(PATNO) %>%
    mutate(across(all_of(vars), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
    ungroup()

  # Pass 2: within-stage mean
  df <- df %>%
    group_by(New_Stage) %>%
    mutate(across(all_of(vars), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
    ungroup()

  # Pass 3: global mean
  df <- df %>%
    mutate(across(all_of(vars), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))

  df
}

blood_vars <- c("WBC", "Platelets", "Neutrophils", "Lymphocytes", "BMI")
metadata_imputed <- impute_hierarchical(metadata_with_chem, blood_vars)

# Derived ratio: neutrophil-to-lymphocyte ratio (NLR)
metadata_imputed$NLR <- metadata_imputed$Neutrophils / metadata_imputed$Lymphocytes

# Log imputation counts
for (var in blood_vars) {
  n_imputed <- sum(is.na(metadata_with_chem[[var]]) & !is.na(metadata_imputed[[var]]))
  cat(paste0("Imputed_", var, "_values"), n_imputed, "\n", file = STATS_FILE, append = TRUE)
}


# ── Pre-filtering summary statistics ──────────────────────────────────────────

write.csv(
  as.data.frame(table(metadata_imputed$NSD_STAGE, useNA = "ifany")),
  file      = "Stats_Metadata_PreFiltering_NSD_Stage_Counts.csv",
  row.names = FALSE
)
write.csv(
  as.data.frame(table(metadata_imputed$COHORT, useNA = "ifany")),
  file      = "Stats_Metadata_PreFiltering_Cohort_Counts.csv",
  row.names = FALSE
)
write.csv(
  as.data.frame(table(metadata_imputed$New_Stage, useNA = "ifany")),
  file      = "Stats_Metadata_PreFiltering_New_Stage_Counts.csv",
  row.names = FALSE
)


# =============================================================================
# 03. DATA CLEANING & SAMPLE HARMONIZATION
# =============================================================================

# Keep only samples present in all three data sources
common_samples <- Reduce(
  intersect,
  list(
    colnames(rnaseq_counts_raw)[-1],
    colnames(mirna_counts_raw)[-1],
    metadata_imputed$Unique_ID
  )
)

rnaseq_counts <- as.data.frame(rnaseq_counts_raw[, c("Geneid", common_samples)])
mirna_counts  <- as.data.frame(mirna_counts_raw[,  c("miRNA",  common_samples)])
metadata      <- as.data.frame(metadata_imputed) %>%
  dplyr::filter(Unique_ID %in% common_samples)

# ── Slim metadata for downstream modelling ────────────────────────────────────

COVARIATES <- c("Unique_ID", "PATNO", "EVENT_ID", "NSD_STAGE", "New_Stage",
                "COHORT", "SEX", "age", "WBC", "Platelets", "Neutrophils",
                "Lymphocytes", "BMI", "NLR", "LEDD")

metadata_model <- metadata[, COVARIATES]

# Impute LEDD:
#   • Healthy controls (COHORT 2) are not on levodopa → set to 0
#   • PD patients: within-patient → within-stage → global mean
metadata_model <- metadata_model %>%
  mutate(LEDD = ifelse(COHORT == 2 & is.na(LEDD), 0, LEDD))

metadata_model <- impute_hierarchical(metadata_model, "LEDD")

n_imputed_ledd <- sum(is.na(metadata[["LEDD"]]) & !is.na(metadata_model[["LEDD"]]))
cat("Imputed_LEDD_values", n_imputed_ledd, "\n", file = STATS_FILE, append = TRUE)

# ── Final sample counts ───────────────────────────────────────────────────────

cat("Final_filtered_samples", nrow(metadata), "\n", file = STATS_FILE, append = TRUE)

sex_counts <- table(metadata$SEX)
cat("Filtered_male_samples",   sex_counts["1"], "\n", file = STATS_FILE, append = TRUE)
cat("Filtered_female_samples", sex_counts["0"], "\n", file = STATS_FILE, append = TRUE)

# ── Post-filtering summary statistics ─────────────────────────────────────────

write.csv(
  as.data.frame(table(metadata$NSD_STAGE, useNA = "ifany")),
  file      = "Stats_Metadata_PostFiltering_NSD_Stage_Counts.csv",
  row.names = FALSE
)
write.csv(
  as.data.frame(table(metadata$New_Stage, useNA = "ifany")),
  file      = "Stats_Metadata_PostFiltering_New_Stage_Counts.csv",
  row.names = FALSE
)
write.csv(
  as.data.frame(table(metadata$COHORT, useNA = "ifany")),
  file      = "Stats_Metadata_PostFiltering_Cohort_Counts.csv",
  row.names = FALSE
)

# ── Persist harmonized objects ────────────────────────────────────────────────

write.csv(metadata,       file = "metadata_filtered.csv",  row.names = FALSE)
write.csv(metadata_model, file = "metadata_useful.csv",    row.names = FALSE)
fwrite(rnaseq_counts, file = "counts_matrix_RNAseq_filtered.csv")
fwrite(mirna_counts,  file = "counts_matrix_miRNA_filtered.csv")

# Tidy workspace before the annotation-heavy step
rm(list = setdiff(ls(), c("rnaseq_counts", "mirna_counts", "metadata", "metadata_model",
                           "STATS_FILE", "impute_hierarchical")))


# =============================================================================
# 04. GENE ANNOTATION, SYMBOL MAPPING & BIOTYPE FILTERING
# =============================================================================

# ── Prepare raw count matrix ──────────────────────────────────────────────────

# Round miRNA counts downward (non-integer values from some quantifiers)
mirna_numeric_cols <- sapply(mirna_counts, is.numeric)
mirna_counts[mirna_numeric_cols] <- floor(mirna_counts[mirna_numeric_cols])

# RNA-seq count matrix as a numeric matrix with ENSG rownames
counts_mat <- as.matrix(
  lapply(rnaseq_counts[, -1], function(x) {
    if (is.character(x) || is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  }) %>% as.data.frame()
)
rownames(counts_mat) <- rnaseq_counts$Geneid

# Strip Ensembl version suffixes (e.g. ENSG00000001234.7 → ENSG00000001234)
ensg_ids <- sub("\\..*$", "", rownames(counts_mat))


# ── Fetch gene annotations from the latest GRCh38 EnsDb ──────────────────────

ah  <- AnnotationHub()
edb <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38"))[[length(query(ah, c("EnsDb", "Homo sapiens", "GRCh38")))]]

gene_annot <- as.data.frame(
  genes(edb, return.type = "DataFrame")[, c("gene_id", "gene_name", "gene_biotype")]
)
colnames(gene_annot) <- c("ensembl_id", "symbol_raw", "biotype")

# Map ENSG IDs in the count matrix to annotation table
annot_idx  <- match(ensg_ids, gene_annot$ensembl_id)
symbol_raw <- gene_annot$symbol_raw[annot_idx]
biotype    <- gene_annot$biotype[annot_idx]


# ── Harmonise gene symbols via HGNChelper ─────────────────────────────────────
# Replaces deprecated aliases with approved HGNC symbols

hgnc_check   <- HGNChelper::checkGeneSymbols(symbol_raw, unmapped.as.na = TRUE)
symbol_clean <- hgnc_check$Suggested.Symbol
symbol_clean[!nzchar(symbol_clean)] <- NA

n_symbols_corrected <- sum(
  hgnc_check$x != hgnc_check$Suggested.Symbol & !is.na(hgnc_check$Suggested.Symbol),
  na.rm = TRUE
)


# ── Build annotation table and resolve duplicate symbols ─────────────────────
# When multiple ENSG IDs map to the same symbol, keep the one with the highest
# mean expression (most likely to be the canonical locus).

mean_expr <- rowMeans(counts_mat, na.rm = TRUE)

annot_table <- data.frame(
  ensg      = ensg_ids,
  symbol    = symbol_clean,
  biotype   = biotype,
  mean_expr = mean_expr,
  stringsAsFactors = FALSE
)

# Drop rows without a valid symbol
annot_valid   <- annot_table[!is.na(annot_table$symbol), ]
counts_valid  <- counts_mat[!is.na(annot_table$symbol), , drop = FALSE]

# Sort by symbol, then descending mean expression, then deduplicate
ord             <- order(annot_valid$symbol, -annot_valid$mean_expr)
annot_sorted    <- annot_valid[ord, ]
counts_sorted   <- counts_valid[ord, , drop = FALSE]
keep_top        <- !duplicated(annot_sorted$symbol)

annot_final     <- annot_sorted[keep_top, ]
counts_symbol   <- counts_sorted[keep_top, , drop = FALSE]
rownames(counts_symbol) <- annot_final$symbol


# ── Biotype-specific count matrices ──────────────────────────────────────────

counts_protein_coding <- counts_symbol[annot_final$biotype == "protein_coding", , drop = FALSE]

nc_mask <- annot_final$biotype %in% c("snRNA", "snoRNA") |
           grepl("lncRNA|lincRNA", annot_final$biotype, ignore.case = TRUE)
counts_non_coding <- counts_symbol[nc_mask, , drop = FALSE]


# ── Annotation summary & persistence ─────────────────────────────────────────

mapping_table <- annot_final[, c("ensg", "symbol", "biotype", "mean_expr")]
write.csv(mapping_table, file = "gene_annotation_mapping_table.csv", row.names = FALSE)

sink("Gene_Annotation_Summary.txt")
cat("=== Gene annotation summary ===\n")
cat("Total input genes          :", nrow(counts_mat), "\n")
cat("Genes with valid symbols   :", nrow(annot_valid),
    sprintf("(%.1f%%)\n", 100 * nrow(annot_valid) / nrow(counts_mat)))
cat("Unique symbols (final)     :", nrow(annot_final), "\n")
cat("Protein-coding genes       :", sum(annot_final$biotype == "protein_coding"),
    sprintf("(%.1f%%)\n", 100 * mean(annot_final$biotype == "protein_coding")))
cat("Selected non-coding genes  :", sum(nc_mask),
    sprintf("(%.1f%%)\n", 100 * mean(nc_mask)))
sink()

# Convert back to data.frame for pipeline consistency
rnaseq_all_symbols <- data.frame(Geneid = rownames(counts_symbol),    counts_symbol,    check.names = FALSE)
rnaseq_pc          <- data.frame(Geneid = rownames(counts_protein_coding), counts_protein_coding, check.names = FALSE)
rnaseq_nc          <- data.frame(Geneid = rownames(counts_non_coding), counts_non_coding, check.names = FALSE)

fwrite(rnaseq_all_symbols, file = "counts_matrix_RNAseq_symbols.csv")
fwrite(rnaseq_pc,          file = "counts_matrix_protein_coding.csv")
fwrite(rnaseq_nc,          file = "counts_matrix_non_coding.csv")


# =============================================================================
# 05. DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

# ── Prepare MIRit input matrices ──────────────────────────────────────────────

# miRNA matrix: rows = miRNA IDs, columns = sample Unique_IDs
rownames(mirna_counts)      <- mirna_counts$miRNA
mirna_mat                   <- as.matrix(mirna_counts[, -1])

# All-gene RNA-seq matrix (protein-coding + selected ncRNA)
rnaseq_mat                  <- as.matrix(rnaseq_all_symbols[, -1])
rownames(rnaseq_mat)        <- rnaseq_all_symbols$Geneid
colnames(rnaseq_mat)        <- sub("^X", "", colnames(rnaseq_mat))   # remove leading X added by R


# ── Sample metadata for MIRit ─────────────────────────────────────────────────
# Columns must align with the miRNA matrix column order.

sample_order <- colnames(mirna_mat)
meta_idx     <- match(sample_order, metadata_model$Unique_ID)

meta_mirit <- data.frame(
  primary  = sample_order,
  mirnaCol = sample_order,
  geneCol  = sample_order,
  disease  = gsub("-", "_", metadata_model$New_Stage[meta_idx]),
  age      = metadata_model$age[meta_idx],
  SEX      = as.factor(metadata_model$SEX[meta_idx]),
  BMI      = metadata_model$BMI[meta_idx],
  WBC      = metadata_model$WBC[meta_idx],
  Platelets= metadata_model$Platelets[meta_idx],
  NLR      = metadata_model$NLR[meta_idx],
  LEDD     = metadata_model$LEDD[meta_idx],
  patient  = metadata_model$Unique_ID[meta_idx]
)


# ── MirnaExperiment object ────────────────────────────────────────────────────

experiment_base <- MirnaExperiment(
  mirnaExpr       = mirna_mat,
  geneExpr        = rnaseq_mat,
  samplesMetadata = meta_mirit,
  pairedSamples   = TRUE
)

# MDS plots to visualise expression variability
mds_genes <- plotDimensions(experiment_base, assay = "genes",    condition = "SEX", title = "MDS – all genes")
mds_mirna <- plotDimensions(experiment_base, assay = "microRNA", condition = "SEX", title = "MDS – miRNA")

mds_combined <- ggpubr::ggarrange(mds_genes, mds_mirna,
                                   nrow = 1, labels = "AUTO", common.legend = TRUE)
ggsave("MDS_plot_all_genes_and_miRNA.png", mds_combined, width = 15, height = 6, dpi = 300)
ggsave("MDS_plot_all_genes_and_miRNA.pdf", mds_combined, width = 15, height = 6)
ggsave("MDS_plot_all_genes_and_miRNA.svg", mds_combined, width = 15, height = 6)


# ── Shared model formula ──────────────────────────────────────────────────────
# All DGE/DEmiRNA calls use this design to control for known confounders.

DE_MODEL <- ~ disease + age + SEX + BMI + WBC + Platelets + NLR + LEDD

# Shared filterByExpr thresholds
FILTER_ARGS <- list(
  min.count       = 5,
  min.total.count = 10,
  large.n         = 10,
  min.prop        = 0.3
)

DE_CUTOFF_PADJ  <- 0.05
DE_CUTOFF_LOGFC <- 0.58    # ~1.5-fold change


# ── Helper: run gene + miRNA DE for one contrast ──────────────────────────────

run_de_analysis <- function(experiment, contrast) {
  gene_exp <- performGeneDE(
    experiment,
    group    = "disease",
    contrast = contrast,
    design   = DE_MODEL,
    pCutoff  = DE_CUTOFF_PADJ,
    logFC    = DE_CUTOFF_LOGFC,
    method   = "edgeR",
    filterByExpr.args = FILTER_ARGS
  )
  mirna_exp <- performMirnaDE(
    experiment,
    group    = "disease",
    contrast = contrast,
    design   = DE_MODEL,
    pCutoff  = DE_CUTOFF_PADJ,
    logFC    = DE_CUTOFF_LOGFC,
    method   = "edgeR",
    filterByExpr.args = FILTER_ARGS
  )
  list(genes = gene_exp, mirna = mirna_exp)
}


# ── Run DE for all four severity contrasts (vs Healthy) ──────────────────────

de_preclinical <- run_de_analysis(experiment_base, "Pre_clinical-Healthy")
de_mild        <- run_de_analysis(experiment_base, "Mild-Healthy")
de_moderate    <- run_de_analysis(experiment_base, "Moderate-Healthy")
de_severe      <- run_de_analysis(experiment_base, "Severe-Healthy")

# Convenience named list for iteration
de_stages <- list(
  Preclinical = de_preclinical,
  Mild        = de_mild,
  Moderate    = de_moderate,
  Severe      = de_severe
)


# ── Extract significant DEGs / DEmiRNAs ───────────────────────────────────────

for (stage in names(de_stages)) {
  deg   <- geneDE(de_stages[[stage]]$genes)
  demi  <- mirnaDE(de_stages[[stage]]$mirna)

  cat(paste0("DE_all_genes_",  stage, "_vs_Healthy"), nrow(deg),  "\n", file = STATS_FILE, append = TRUE)
  cat(paste0("DE_miRNAs_all_", stage, "_vs_Healthy"), nrow(demi), "\n", file = STATS_FILE, append = TRUE)

  fwrite(as.data.frame(deg),  file = paste0("DE_all_genes_",  stage, "_vs_Healthy.csv"))
  fwrite(as.data.frame(demi), file = paste0("DE_miRNAs_all_", stage, "_vs_Healthy.csv"))
}

# Full (non-significant-filtered) result tables – needed for volcano plots
for (stage in names(de_stages)) {
  full_genes <- geneDE(de_stages[[stage]]$genes, onlySignificant = FALSE)
  full_mirna <- mirnaDE(de_stages[[stage]]$mirna, onlySignificant = FALSE)

  fwrite(as.data.frame(full_genes), file = paste0("DGE_results_all_genes_combined_",  stage, "_vs_Healthy.csv"))
  fwrite(as.data.frame(full_mirna), file = paste0("DGE_results_all_miRNAs_combined_", stage, "_vs_Healthy.csv"))
}


# ── Volcano plots: 2 × 4 grid (mRNA rows × severity columns) ─────────────────

STAGES     <- c("Preclinical", "Mild", "Moderate", "Severe")
FC_THR     <- DE_CUTOFF_LOGFC
PADJ_THR   <- DE_CUTOFF_PADJ

read_de_table <- function(filepath, stage, omics_layer) {
  df <- read_csv(filepath, show_col_types = FALSE)

  required <- c("ID", "logFC", "adj.P.Val")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s", filepath, paste(missing, collapse = ", ")))
  }

  df %>%
    transmute(
      ID        = as.character(ID),
      logFC     = as.numeric(logFC),
      adj.P.Val = as.numeric(`adj.P.Val`),
      condition = stage,
      omics     = omics_layer
    ) %>%
    dplyr::filter(!is.na(logFC), !is.na(adj.P.Val))
}

annotate_de_direction <- function(df) {
  df %>%
    mutate(
      neg_log10_padj = -log10(adj.P.Val),
      direction = case_when(
        adj.P.Val < PADJ_THR & logFC >  FC_THR ~ "Up",
        adj.P.Val < PADJ_THR & logFC < -FC_THR ~ "Down",
        TRUE ~ "NotSig"
      )
    )
}

make_volcano <- function(df, title = NULL) {
  ggplot(df, aes(x = logFC, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.65, size = 1.1) +
    geom_vline(xintercept = c(-FC_THR, FC_THR), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PADJ_THR),   linetype = "dashed", linewidth = 0.4) +
    scale_color_manual(values = c(NotSig = "grey70", Down = "blue", Up = "red")) +
    labs(
      title = title,
      x     = expression(log[2]("Fold change")),
      y     = expression(-log[10]("adj.P.Val"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(face = "bold", size = 10, hjust = 0.5),
      panel.grid.minor   = element_blank()
    )
}

mrna_volcano_list <- lapply(STAGES, function(s) {
  filepath <- paste0("DGE_results_all_genes_combined_",  s, "_vs_Healthy.csv")
  read_de_table(filepath, s, "mRNA") %>% annotate_de_direction() %>% make_volcano(title = s)
})

mirna_volcano_list <- lapply(STAGES, function(s) {
  filepath <- paste0("DGE_results_all_miRNAs_combined_", s, "_vs_Healthy.csv")
  read_de_table(filepath, s, "miRNA") %>% annotate_de_direction() %>% make_volcano(title = s)
})

volcano_grid <- (wrap_plots(mrna_volcano_list,  nrow = 1) /
                 wrap_plots(mirna_volcano_list, nrow = 1))

ggsave("volcano_grid_mrna_mirna.pdf", volcano_grid, width = 14, height = 6.5)
ggsave("volcano_grid_mrna_mirna.png", volcano_grid, width = 14, height = 6.5, dpi = 300)
ggsave("volcano_grid_mrna_mirna.svg", volcano_grid, width = 14, height = 6.5)

# Per-stage volcano panels (gene + miRNA side by side)
for (stage in STAGES) {
  panel <- ggpubr::ggarrange(
    mrna_volcano_list[[which(STAGES == stage)]],
    mirna_volcano_list[[which(STAGES == stage)]],
    nrow = 1, labels = "AUTO", common.legend = TRUE
  )
  ggsave(paste0("Volcano_plot_AllGenes_", stage, "_vs_Healthy.png"),
         panel, width = 15, height = 6, dpi = 300)
}


# ── UpSet plot: DEG overlap across severity stages ───────────────────────────

degs_by_stage <- list(
  `Pre-clinical` = rownames(geneDE(de_preclinical$genes)),
  Mild           = rownames(geneDE(de_mild$genes)),
  Moderate       = rownames(geneDE(de_moderate$genes)),
  Severe         = rownames(geneDE(de_severe$genes))
)

upset_data <- fromList(degs_by_stage)
upset_plot <- upset(
  upset_data,
  nsets        = 4,
  nintersects  = NA,
  order.by     = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "skyblue",
  text.scale     = c(2, 2, 2, 1.5, 2, 2)
)

png("UpSet_plot_Genes_DEGs.png", width = 10, height = 6, units = "in", res = 300)
print(upset_plot)
dev.off()

svg("UpSet_plot_Genes_DEGs.svg", width = 10, height = 6)
print(upset_plot)
dev.off()


# =============================================================================
# 06. miRNA TARGET RETRIEVAL & miRNA–mRNA CORRELATION
# =============================================================================

# ── Retrieve predicted/validated miRNA targets ────────────────────────────────

de_preclinical$mirna <- getTargets(de_preclinical$mirna)
de_mild$mirna        <- getTargets(de_mild$mirna)
de_moderate$mirna    <- getTargets(de_moderate$mirna)
de_severe$mirna      <- getTargets(de_severe$mirna)

# Persist target–miRNA pair tables
fwrite(as.data.frame(mirnaTargets(de_preclinical$mirna)), file = "miRNA_targets_AllGenes_Preclinical_vs_Healthy.csv")
fwrite(as.data.frame(mirnaTargets(de_mild$mirna)),        file = "miRNA_targets_AllGenes_Mild_vs_Healthy.csv")
fwrite(as.data.frame(mirnaTargets(de_moderate$mirna)),    file = "miRNA_targets_AllGenes_Moderate_vs_Healthy.csv")
fwrite(as.data.frame(mirnaTargets(de_severe$mirna)),      file = "miRNA_targets_AllGenes_Severe_vs_Healthy.csv")

cat("miRNA_target_pairs_AllGenes",
    nrow(mirnaTargets(de_preclinical$mirna)), "\n",
    file = STATS_FILE, append = TRUE)


# ── Parallelisation (Linux multi-core) ───────────────────────────────────────

register(MulticoreParam(workers = 16))

BATCH_COVARIATES <- c("age", "BMI", "WBC", "Platelets", "NLR", "LEDD")


# ── Helper: assemble joint experiment, correct batch effects, run correlation ─
# Merges miRNA and gene DE results, removes technical variation, then tests
# Spearman correlation between each DEmiRNA and its predicted target genes.

run_mirna_mrna_integration <- function(mirna_exp, gene_exp, bpparam = MulticoreParam(workers = 16)) {

  # Transfer gene DE results into the miRNA experiment object
  joint_exp <- mirna_exp
  slot(joint_exp, "geneDE") <- slot(gene_exp, "geneDE")

  # Batch correction for miRNA assay
  joint_exp <- batchCorrection(
    joint_exp,
    assay       = "microRNA",
    batch       = "patient",
    batch2      = "SEX",
    covariates  = BATCH_COVARIATES,
    includeWsva = TRUE,
    n.sv        = 2,
    weight.by.sd = TRUE
  )

  # Batch correction for gene assay
  joint_exp <- batchCorrection(
    joint_exp,
    assay       = "genes",
    batch       = "patient",
    batch2      = "SEX",
    covariates  = BATCH_COVARIATES,
    includeWsva = TRUE,
    n.sv        = 2,
    weight.by.sd = TRUE
  )

  # Spearman correlation integration
  joint_exp <- mirnaIntegration(joint_exp, test = "correlation", BPPARAM = bpparam)

  integration(joint_exp)
}


# ── Run integration for each severity stage ───────────────────────────────────

integration_preclinical <- run_mirna_mrna_integration(de_preclinical$mirna, de_preclinical$genes)
fwrite(as.data.frame(integration_preclinical), file = "miRNA_mRNA_correlation_AllGenes_Preclinical_vs_Healthy.csv")

integration_mild <- run_mirna_mrna_integration(de_mild$mirna, de_mild$genes)
fwrite(as.data.frame(integration_mild), file = "miRNA_mRNA_correlation_AllGenes_Mild_vs_Healthy.csv")

# Moderate: SerialParam used to avoid memory issues with larger sample set
integration_moderate <- run_mirna_mrna_integration(de_moderate$mirna, de_moderate$genes, bpparam = SerialParam())
fwrite(as.data.frame(integration_moderate), file = "miRNA_mRNA_correlation_AllGenes_Moderate_vs_Healthy.csv")

integration_severe <- run_mirna_mrna_integration(de_severe$mirna, de_severe$genes)
fwrite(as.data.frame(integration_severe), file = "miRNA_mRNA_correlation_AllGenes_Severe_vs_Healthy.csv")


# ── Log correlation pair counts ───────────────────────────────────────────────

cat("miRNA_mRNA_correlation_pairs_AllGenes_Preclinical_vs_Healthy", nrow(integration_preclinical), "\n", file = STATS_FILE, append = TRUE)
cat("miRNA_mRNA_correlation_pairs_AllGenes_Mild_vs_Healthy",        nrow(integration_mild),        "\n", file = STATS_FILE, append = TRUE)
cat("miRNA_mRNA_correlation_pairs_AllGenes_Moderate_vs_Healthy",    nrow(integration_moderate),    "\n", file = STATS_FILE, append = TRUE)
cat("miRNA_mRNA_correlation_pairs_AllGenes_Severe_vs_Healthy",      nrow(integration_severe),      "\n", file = STATS_FILE, append = TRUE)


# =============================================================================
# END OF ANALYSIS
# =============================================================================
