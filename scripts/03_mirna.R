# =============================================================================
# miRNA analysis (Figure 4a-c, Supplementary Tables S2-S3)
# Paired DESeq2 (~ subject + condition) with apeglm shrinkage, and multiMiR
# targets among the meta-DEGs. Run meta_final.R first.
# =============================================================================
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(multiMiR)

out_dir <- "output_revised/mirna/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- "output_revised/mirna/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# Inputs from the meta-analysis
meta_degs_file <- "output_revised/meta/tables/meta_degs.csv"
study_tab_file <- "output_revised/meta/tables/study_cohort_toptab.tsv"
stopifnot(file.exists(meta_degs_file), file.exists(study_tab_file))

COL_BASELINE <- "#D4A96A"
COL_WEEK4    <- "tomato3"
COL_UP       <- "#C1121F"
COL_DOWN     <- "#1D4E89"
COL_NS       <- "#ADB5BD"

# ── DESeq2 ────────────────────────────────────────────────────────────────────
# The file starts with a byte-order mark
miRNA_counts <- read.table("data/mirna/counts_paired.csv", header=T, sep=",", fileEncoding="UTF-8-BOM")
stopifnot("Feature" %in% colnames(miRNA_counts))
rownames(miRNA_counts) <- miRNA_counts$Feature
miRNA_counts$Feature   <- NULL

# Timepoint and subject from the column names
condition <- sub("^X?\\d+_", "", colnames(miRNA_counts))
subject   <- sub("^X?(\\d+)_.*", "\\1", colnames(miRNA_counts))
stopifnot(all(condition %in% c("baseline", "treatment")),
          sum(condition == "baseline") == 11, sum(condition == "treatment") == 11,
          all(table(subject) == 2))
col_data  <- data.frame(subject=subject, condition=condition)

dds <- DESeqDataSetFromMatrix(countData=miRNA_counts, colData=col_data,
                              design=~subject + condition)
dds <- DESeq(dds)

res_shrunk <- as.data.frame(lfcShrink(dds,
                                      coef="condition_treatment_vs_baseline", type="apeglm"))
res_shrunk <- res_shrunk[complete.cases(res_shrunk) & res_shrunk$baseMean > 50, ]
res_shrunk$sig <- ifelse(res_shrunk$log2FoldChange > 0  & res_shrunk$padj < 0.05, "Up",
                         ifelse(res_shrunk$log2FoldChange < 0  & res_shrunk$padj < 0.05, "Down", "NS"))

deg_list   <- res_shrunk[abs(res_shrunk$log2FoldChange) > 0.58 & res_shrunk$padj < 0.05, ]
sig_mirnas <- rownames(deg_list)

meta_degs <- read.csv(meta_degs_file)
study_tab <- read.delim(study_tab_file, sep="\t")

# ── Panel A: ComplexHeatmap ───────────────────────────────────────────────────
vsd     <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_mat <- assay(vsd)

heat_mat           <- vsd_mat[sig_mirnas, ]
rownames(heat_mat) <- gsub("hsa-", "", sig_mirnas)
heat_mat           <- t(scale(t(heat_mat)))

col_order  <- c(which(col_data$condition=="baseline"),
                which(col_data$condition=="treatment"))
heat_mat   <- heat_mat[, col_order]
cond_order <- col_data$condition[col_order]
col_fun    <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "white", COL_UP))

ht <- Heatmap(heat_mat,
              col = col_fun,
              top_annotation = HeatmapAnnotation(
                Timepoint = cond_order,
                col       = list(Timepoint=c(baseline=COL_BASELINE, treatment=COL_WEEK4)),
                annotation_name_gp      = gpar(fontsize=11),
                annotation_legend_param = list(title_gp=gpar(fontsize=11),
                                               labels_gp=gpar(fontsize=10))
              ),
              right_annotation = rowAnnotation(
                Direction = ifelse(deg_list$log2FoldChange > 0, "Up", "Down"),
                col       = list(Direction=c(Up=COL_UP, Down=COL_DOWN)),
                annotation_name_gp      = gpar(fontsize=11),
                annotation_legend_param = list(title_gp=gpar(fontsize=11),
                                               labels_gp=gpar(fontsize=10)),
                width = unit(0.8, "cm")
              ),
              cluster_columns   = FALSE,
              cluster_rows      = TRUE,
              show_column_names = FALSE,
              row_names_gp      = gpar(fontsize=11),
              column_split      = factor(cond_order, levels=c("baseline","treatment")),
              column_title      = NULL,
              heatmap_legend_param = list(title="Z-score",
                                          title_gp=gpar(fontsize=11),
                                          labels_gp=gpar(fontsize=10)),
              border = FALSE)

draw(ht, merge_legend=TRUE)
# Save Figure 4a
pdf(file.path(fig_dir, "Fig4a_heatmap.pdf"), width = 8, height = 4); draw(ht, merge_legend=TRUE); dev.off()
png(file.path(fig_dir, "Fig4a_heatmap.png"), width = 8, height = 4, units = "in", res = 300); draw(ht, merge_legend=TRUE); dev.off()

# ── Panel B: Volcano (tRF style) ─────────────────────────────────────────────
res_df_plot       <- res_shrunk
res_df_plot$mirna <- gsub("hsa-", "", rownames(res_df_plot))
res_df_plot$label <- ifelse(res_df_plot$sig != "NS", res_df_plot$mirna, NA)

pVolc <- ggplot(res_df_plot, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(shape=1, size=1.8, color="grey60", stroke=0.4) +
  geom_point(data=subset(res_df_plot, sig=="Up"),
             shape=16, size=2.5, color=COL_UP) +
  geom_point(data=subset(res_df_plot, sig=="Down"),
             shape=16, size=2.5, color=COL_DOWN) +
  geom_text_repel(aes(label=label), size=3.2,
                  na.rm=TRUE, max.overlaps=20,
                  segment.color="grey60", segment.size=0.3) +
  geom_vline(xintercept=c(-0.58, 0.58), linetype="dashed",
             color="grey50", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             color="grey50", linewidth=0.4) +
  labs(x=expression("log"[2]*"FC (Week 4 vs Baseline)"),
       y=expression("-log"[10]*"(adj. p-value)")) +
  theme_bw(base_size=13) +
  theme(axis.title       = element_text(),
        axis.text        = element_text(),
        panel.grid.major = element_line(color="white", linewidth=0.5),
        panel.background = element_rect(fill="white"))

print(pVolc)
ggsave(file.path(fig_dir, "Fig4b_volcano.pdf"), pVolc, width = 6, height = 5)
ggsave(file.path(fig_dir, "Fig4b_volcano.png"), pVolc, width = 6, height = 5, dpi = 300)

# ── multiMiR targets ──────────────────────────────────────────────────────────
targets_pred <- get_multimir(mirna=sig_mirnas, table="predicted", summary=FALSE)
targets_val  <- get_multimir(mirna=sig_mirnas, table="validated", summary=FALSE)

pred_df          <- targets_pred@data[, c("mature_mirna_id","target_symbol")]
pred_df$evidence <- "Predicted"
val_df           <- targets_val@data[, c("mature_mirna_id","target_symbol")]
val_df$evidence  <- "Validated"

targets_df <- rbind(pred_df, val_df)

# Databases supporting each miRNA-target pair
db_all <- rbind(targets_pred@data[, c("mature_mirna_id","target_symbol","database")],
                targets_val@data[,  c("mature_mirna_id","target_symbol","database")])
db_summary <- aggregate(database ~ mature_mirna_id + target_symbol, data = db_all,
                        FUN = function(x) paste(sort(unique(x)), collapse = ";"))
cat("Predicted-target databases queried:",
    paste(sort(unique(targets_pred@data$database)), collapse = ", "), "\n")

targets_df <- targets_df[order(targets_df$evidence == "Validated", decreasing=TRUE), ]
targets_df <- targets_df[!duplicated(targets_df[, c("mature_mirna_id","target_symbol")]), ]

targets_meta             <- targets_df[targets_df$target_symbol %in% meta_degs$hgnc_symbol, ]
targets_meta$mirna_clean <- gsub("hsa-", "", targets_meta$mature_mirna_id)

# ── Panel C: Tile plot ────────────────────────────────────────────────────────
tile_df            <- targets_meta[, c("mirna_clean","target_symbol","evidence")]
tile_df            <- tile_df[!duplicated(tile_df[, c("mirna_clean","target_symbol")]), ]
tile_df$evidence_f <- factor(tile_df$evidence, levels=c("Validated","Predicted"))

# order rows by total hits, cols by miRNA lfc
gene_order  <- names(sort(table(tile_df$target_symbol), decreasing=TRUE))
mirna_order <- gsub("hsa-", "", rownames(deg_list)[order(deg_list$log2FoldChange)])

tile_df$target_symbol <- factor(tile_df$target_symbol, levels=gene_order)
tile_df$mirna_clean   <- factor(tile_df$mirna_clean,   levels=mirna_order)

# pastel validated/predicted colors
COL_VAL  <- "#F4A7A7"
COL_PRED <- "#C5D8F0"

pTile <- ggplot(tile_df, aes(x=mirna_clean, y=target_symbol, fill=evidence_f)) +
  geom_tile(color="white", linewidth=0.7) +
  scale_fill_manual(values=c(Validated=COL_VAL, Predicted=COL_PRED), name="Evidence") +
  labs(x=NULL, y=NULL) +
  theme_classic(base_size=12) +
  theme(axis.text.x     = element_text(angle=40, hjust=1, size=11),
        axis.text.y     = element_text(size=10),
        legend.position = "top",
        legend.title    = element_text(size=11),
        legend.text     = element_text(size=10),
        panel.border    = element_rect(color="black", fill=NA, linewidth=0.5))

print(pTile)
ggsave(file.path(fig_dir, "Fig4c_targets.pdf"), pTile, width = 6, height = 7)
ggsave(file.path(fig_dir, "Fig4c_targets.png"), pTile, width = 6, height = 7, dpi = 300)

# ── miRNA-mRNA pairs ──────────────────────────────────────────────────────────
mirna_lfc      <- data.frame(mature_mirna_id=rownames(deg_list),
                             mirna_lfc=deg_list$log2FoldChange)
study_lfc      <- study_tab[, c("hgnc_symbol","logFC")]
colnames(study_lfc)[2] <- "mrna_lfc"

inv_df <- merge(targets_meta, mirna_lfc,  by="mature_mirna_id")
inv_df <- merge(inv_df,       study_lfc,  by.x="target_symbol", by.y="hgnc_symbol")
inv_df <- inv_df[!duplicated(inv_df[, c("mature_mirna_id","target_symbol")]), ]
inv_df$mirna_clean <- gsub("hsa-", "", inv_df$mature_mirna_id)
inv_df$inverse     <- (inv_df$mirna_lfc > 0 & inv_df$mrna_lfc < 0) |
  (inv_df$mirna_lfc < 0 & inv_df$mrna_lfc > 0)

inv_only <- inv_df[inv_df$inverse, ]

cat(sprintf("Total pairs: %d | Inverse: %d (%.1f%%)\n",
            nrow(inv_df), nrow(inv_only), 100*nrow(inv_only)/nrow(inv_df)))

norm_counts          <- counts(dds, normalized=TRUE)
log_counts           <- log2(norm_counts[sig_mirnas, , drop=FALSE] + 1)
rownames(log_counts) <- gsub("hsa-", "", rownames(log_counts))

expr_long            <- as.data.frame(t(log_counts))
expr_long$subject    <- factor(col_data$subject)
expr_long$timepoint  <- factor(col_data$condition,
                               levels=c("baseline","treatment"),
                               labels=c("Baseline","Week 4"))
expr_long <- pivot_longer(expr_long, -c(subject, timepoint),
                          names_to="miRNA", values_to="expr")

pSlope <- ggplot(expr_long, aes(x=timepoint, y=expr, group=subject)) +
  geom_line(color="grey60", linewidth=0.5) +
  geom_point(aes(color=timepoint), size=2.5) +
  facet_wrap(~miRNA, scales="free_y", nrow=2,
             labeller=labeller(miRNA=label_wrap_gen(15))) +
  scale_color_manual(values=c(Baseline=COL_BASELINE, "Week 4"=COL_WEEK4)) +
  labs(x=NULL, y=expression(Log[2]~Normalized~Counts), color=NULL) +
  theme_classic(base_size=11) +
  theme(strip.text       = element_text(face="italic", size=9),
        strip.background = element_blank(),
        legend.position  = "top",
        panel.border     = element_rect(color="black", fill=NA, linewidth=0.5))

print(pSlope)
ggsave(file.path(fig_dir, "SuppFig_mirna_paired_slopes.pdf"), pSlope, width = 9, height = 5)
ggsave(file.path(fig_dir, "SuppFig_mirna_paired_slopes.png"), pSlope, width = 9, height = 5, dpi = 300)

# ── Write miRNA DE table ──────────────────────────────────────────────────────
mirna_table <- res_shrunk
mirna_table$miRNA <- gsub("hsa-", "", rownames(mirna_table))
mirna_table <- mirna_table[, c("miRNA", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "sig")]
mirna_table <- mirna_table[order(mirna_table$padj), ]

write.csv(mirna_table, file.path(out_dir, "mirna_de_results.csv"), row.names=FALSE)

# ── Write miRNA-mRNA pairs table ──────────────────────────────────────────────
pairs_table <- merge(targets_meta, mirna_lfc, by="mature_mirna_id")
pairs_table <- merge(pairs_table, study_lfc, by.x="target_symbol", by.y="hgnc_symbol")
pairs_table$inverse <- (pairs_table$mirna_lfc > 0 & pairs_table$mrna_lfc < 0) |
  (pairs_table$mirna_lfc < 0 & pairs_table$mrna_lfc > 0)
# Add the meta-analysis log2FC
pairs_table <- merge(pairs_table, meta_degs[, c("hgnc_symbol","logFC_meta")],
                     by.x="target_symbol", by.y="hgnc_symbol")
stopifnot(all(sign(pairs_table$mrna_lfc) == sign(pairs_table$logFC_meta)))
# Supporting database(s) for each pair
pairs_table <- merge(pairs_table, db_summary, by=c("mature_mirna_id","target_symbol"), all.x=TRUE)
pairs_table$mirna_clean <- gsub("hsa-", "", pairs_table$mature_mirna_id)
pairs_table <- pairs_table[, c("mirna_clean", "target_symbol", "evidence", "database",
                               "mirna_lfc", "mrna_lfc", "logFC_meta", "inverse")]
colnames(pairs_table)[colnames(pairs_table) == "mrna_lfc"]   <- "mrna_lfc_study_cohort"
colnames(pairs_table)[colnames(pairs_table) == "logFC_meta"] <- "mrna_lfc_meta"
pairs_table <- pairs_table[order(pairs_table$evidence == "Validated",
                                 decreasing=TRUE), ]
pairs_table$inverse <- ifelse(pairs_table$inverse, "Yes", "No")
cat(sprintf("Supp Table 3: %d pairs, %d validated, %d inverse\n", nrow(pairs_table),
            sum(pairs_table$evidence == "Validated"), sum(pairs_table$inverse == "Yes")))
write.csv(pairs_table, file.path(out_dir, "mirna_mrna_pairs.csv"), row.names=FALSE)

