suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(ggrepel)
})

DATA_DIR    <- "/home/SexyThighs/data/harmonized"
RESULTS_DIR <- file.path(DATA_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE)

cat("Loading data...\n")
rna_raw     <- readRDS(file.path(DATA_DIR, "rna_matrix.rds"))
gene_map    <- readRDS(file.path(DATA_DIR, "gene_symbol_map.rds"))
survival_df <- readRDS(file.path(DATA_DIR, "survival_df.rds"))

## 1. Ensembl -> HGNC
cat("Mapping gene symbols...\n")
rownames(rna_raw) <- sub("\\..*", "", rownames(rna_raw))

map_cols <- colnames(gene_map)
ens_col  <- map_cols[grep("ensembl|ensg", map_cols, ignore.case = TRUE)][1]
sym_col  <- map_cols[grep("symbol|hgnc|gene_name", map_cols, ignore.case = TRUE)][1]
cat(sprintf("  Map columns: '%s' -> '%s'\n", ens_col, sym_col))

gmap <- gene_map[, c(ens_col, sym_col)]
colnames(gmap) <- c("ensembl_id", "symbol")
gmap <- gmap[!is.na(gmap$symbol) & gmap$symbol != "", ]

keep     <- intersect(rownames(rna_raw), gmap$ensembl_id)
rna_sub  <- rna_raw[keep, ]
gmap_sub <- gmap[match(keep, gmap$ensembl_id), ]

vars     <- apply(rna_sub, 1, var)
ord      <- order(vars, decreasing = TRUE)
rna_sub  <- rna_sub[ord, ]
gmap_sub <- gmap_sub[ord, ]
keep2    <- !duplicated(gmap_sub$symbol) & !is.na(gmap_sub$symbol)
rna_sub  <- rna_sub[keep2, ]
rownames(rna_sub) <- gmap_sub$symbol[keep2]
cat(sprintf("  Final matrix: %d genes x %d patients\n", nrow(rna_sub), ncol(rna_sub)))

## 2. Hallmark gene sets
cat("Fetching MSigDB Hallmark gene sets...\n")
h_df   <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- split(h_df$gene_symbol, h_df$gs_name)
cat(sprintf("  %d Hallmark gene sets\n", length(h_list)))

## 3. ssGSEA
cat("Running ssGSEA...\n")
es <- tryCatch({
  params <- ssgseaParam(rna_sub, h_list)
  gsva(params, verbose = FALSE)
}, error = function(e) {
  gsva(rna_sub, h_list, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
})
cat(sprintf("  Result: %d pathways x %d patients\n", nrow(es), ncol(es)))

## 4. Shannon entropy per patient
cat("Computing entropy...\n")
shannon_entropy <- function(x) {
  x <- x - min(x) + 1e-9
  p <- x / sum(x)
  -sum(p * log2(p))
}
entropy_vec <- apply(es, 2, shannon_entropy)

## 5. Merge
common_pts <- intersect(names(entropy_vec), survival_df$patient)
cat(sprintf("  Overlapping patients: %d\n", length(common_pts)))

entropy_df <- data.frame(patient = names(entropy_vec), entropy = entropy_vec,
                         stringsAsFactors = FALSE)
full_df <- survival_df |>
  filter(patient %in% common_pts) |>
  left_join(entropy_df, by = "patient")

## 6. Correlation
cat("\n── Bottleneck Score vs Entropy ──\n")
ct <- cor.test(full_df$bottleneck_score, full_df$entropy, method = "spearman")
cat(sprintf("  Spearman rho = %.3f, p = %.4g\n", ct$estimate, ct$p.value))

## 7. Cox: entropy ~ survival
cat("\n── Cox: entropy + age + stage ──\n")
cox_e <- coxph(Surv(OS_days, OS_status) ~ entropy + age + stage_clean, data = full_df)
print(summary(cox_e)$coefficients)

## 8. Per-pathway disorder
cat("\n── Per-pathway disorder ──\n")
es_common <- es[, common_pts]
bs_common <- full_df$bottleneck_score[match(common_pts, full_df$patient)]

pathway_disorder <- apply(es_common, 1, function(ps) {
  cor(abs(ps - median(ps)), bs_common, method = "spearman")
})

disorder_df <- data.frame(pathway = names(pathway_disorder),
                          spearman_r = pathway_disorder,
                          stringsAsFactors = FALSE) |>
  arrange(desc(spearman_r))

cat("Top 10 disordered pathways in high-bottleneck patients:\n")
print(head(disorder_df, 10))

## 9. Save
write.csv(full_df |> select(patient, bottleneck_score, score_group, entropy,
                              OS_days, OS_status, age, stage_clean),
          file.path(RESULTS_DIR, "hallmark_entropy_full.csv"), row.names = FALSE)
write.csv(disorder_df, file.path(RESULTS_DIR, "hallmark_pathway_disorder.csv"), row.names = FALSE)
saveRDS(es, file.path(RESULTS_DIR, "ssgsea_hallmark_matrix.rds"))

## 10. Figures
r_label <- sprintf("Spearman rho = %.2f\np = %.3g", ct$estimate, ct$p.value)

p1 <- ggplot(full_df, aes(x = bottleneck_score, y = entropy, color = score_group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, label = r_label, hjust = 1.1, vjust = 1.5, size = 3.5) +
  scale_color_manual(values = c("High" = "#c0392b", "Low" = "#2980b9")) +
  labs(x = "Composite Bottleneck Score", y = "Transcriptome Entropy (bits)",
       color = "Score Group",
       title = "Bottleneck Score vs Pathway-Level Transcriptome Entropy") +
  theme_classic(base_size = 12)

ggsave(file.path(RESULTS_DIR, "fig_entropy_vs_bottleneck.pdf"), p1, width = 6, height = 4.5)

top_n  <- head(disorder_df, 10)
bot_n  <- tail(disorder_df, 5)
plot_df <- bind_rows(top_n, bot_n) |>
  mutate(pathway_clean = gsub("HALLMARK_", "", pathway),
         direction = ifelse(spearman_r > 0, "More disordered (high bottleneck)",
                                            "More ordered (low bottleneck)"))

p2 <- ggplot(plot_df, aes(x = reorder(pathway_clean, spearman_r),
                           y = spearman_r, fill = direction)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("More disordered (high bottleneck)" = "#c0392b",
                                "More ordered (low bottleneck)"    = "#2980b9")) +
  labs(x = NULL, y = "Spearman r (pathway disorder ~ bottleneck score)",
       fill = NULL, title = "Hallmark Pathways Associated with Transcriptome Disorder") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(RESULTS_DIR, "fig_pathway_disorder.pdf"), p2, width = 8, height = 6)

cat("\nDone. Results in", RESULTS_DIR, "\n")
