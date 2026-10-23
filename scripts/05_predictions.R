# =============================================================================
# EARLY5 k-TSP prediction of week-4 culture conversion (Figure 3, Table 2,
# Supplementary Tables S5-S7, Supplementary Figures S1-S2)
# k-TSP with LOOCV, comparison with Sweeney3 and RISK6, logistic regression
# with MGIT TTP, and Cox models for time to conversion.
# =============================================================================
suppressPackageStartupMessages({
  library(switchBox); library(pROC); library(limma); library(edgeR)
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(ggpubr); library(rstatix)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(TBSignatureProfiler); library(SummarizedExperiment)
  library(survival); library(patchwork)
  library(logistf)
})
set.seed(42)

dir.create("output_revised/predictions/tables",  recursive=TRUE, showWarnings=FALSE)
dir.create("output_revised/predictions/figures", recursive=TRUE, showWarnings=FALSE)

COL_UP   <- "#C1121F"   # k-TSP LOOCV -- primary, validated result
COL_DOWN <- "#1D4E89"   # Sweeney_OD_3
COL_R6   <- "#7B8B6F"   # PennNich_RISK_6
COL_FAST <- "#D95F5F"   # fast converters / predicted fast
COL_SLOW <- "#7F7F7F"   # slow converters / predicted slow
week_colors <- c("Week 4"="#D95F5F","Week 8"="#E8A87C","Week 12"="#7FBFCF","Week 24"="#9B72CF")

SCORE_BREAKS <- seq(0, 1, 0.25)   # shared y-axis for every EARLY5-score panel
SCORE_LAB    <- "EARLY5 score"

pub_theme <- theme_classic(base_size=14) +
  theme(axis.title=element_text(color="black"), axis.text=element_text(color="black"),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.6),
        panel.grid.major=element_line(color="grey92", linewidth=0.3), panel.grid.minor=element_blank())

# p-value format for figure annotations
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# =============================================================================
# 1. EARLY5 gene pairs
# =============================================================================

EARLY5 <- data.frame(
  gene1 = c("DIS3L","COL9A2","PIEZO1","PID1","PDE4DIP"),
  gene2 = c("DOCK4","FAM241A","H2BC12","LYPLA1","GRB10")
)
EARLY5$pair <- paste0(EARLY5$gene1," > ",EARLY5$gene2)

# =============================================================================
# 2. CLINICAL METADATA -> DEVELOPMENT COHORT
# =============================================================================

metadata <- read.table("metadata/gse89403/SraRunTable.csv", sep=",", header=TRUE)
metadata$tech_group <- paste(metadata$subject, metadata$time, sep="_")

clinical <- metadata %>%
  group_by(tech_group) %>%
  summarise(subject=subject[1], timepoint=time[1],
            timetonegativity=timetonegativity[1], treatmentresult=treatmentresult[1],
            mgit=mean(mgit, na.rm=TRUE), xpert=mean(xpert, na.rm=TRUE), .groups="drop") %>%
  mutate(mgit=ifelse(is.nan(mgit), NA, mgit), xpert=ifelse(is.nan(xpert), NA, xpert))

baseline <- clinical %>%
  filter(timepoint=="DX", !is.na(timetonegativity)) %>%
  mutate(fast_responder=factor(ifelse(timetonegativity=="Week04","Fast","Slow"), levels=c("Slow","Fast")))

cat(sprintf("Baseline patients: %d\n", nrow(baseline))); print(table(baseline$fast_responder))

# =============================================================================
# 3. RNA-SEQ COUNTS -> VOOM-NORMALIZED EXPRESSION
# =============================================================================

count_files <- list.files("data/mrna/gse89403/counts", full.names=TRUE)
dat <- lapply(count_files, function(f) {
  df <- fread(f); df_sub <- df[, c(1,7)]
  colnames(df_sub) <- c("GeneID", gsub(".*_([0-9]+)_.*","\\1",colnames(df)[7]))
  df_sub
})
mat <- do.call(cbind, lapply(dat, function(x) x[,-1]))
# Order columns numerically; columns are matched by SRR run ID below
mat <- as.data.frame(mat)[, order(suppressWarnings(as.numeric(colnames(mat))))]
rownames(mat) <- dat[[1]]$GeneID
colnames(mat) <- gsub(".*/(SRR[0-9]+)\\.bam","\\1",colnames(mat))

meta_sub   <- metadata[, c("Run","subject","time")]
has_both   <- with(meta_sub, table(subject,time)) >= 1
valid_subj <- rownames(has_both[rowSums(has_both)==2, ])
meta_sub   <- meta_sub[meta_sub$subject %in% valid_subj, ]
meta_sub$tech_group <- paste(meta_sub$subject, meta_sub$time, sep="_")

summed_counts <- t(rowsum(t(mat[, meta_sub$Run]), group=meta_sub$tech_group))

sample_info <- data.frame(sample_name=colnames(summed_counts)) %>%
  separate(sample_name, into=c("subject","timepoint"), sep="_", extra="merge")
sample_info$sample_name <- colnames(summed_counts)
sample_info$timepoint   <- factor(sample_info$timepoint, levels=c("DX","week_4"))

keep   <- rowSums(cpm(summed_counts) > 1) >= 10
dge    <- calcNormFactors(DGEList(counts=summed_counts[keep, ]))
design <- model.matrix(~subject + timepoint, data=sample_info)
v      <- voom(dge, design, plot=FALSE)

# =============================================================================
# 4. ENSEMBL -> HGNC
# =============================================================================

ens_ids  <- sub("\\..*","", rownames(v$E))
hgnc_vec <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=ens_ids, column="SYMBOL",
                                  keytype="ENSEMBL", multiVals="first")
keep_rows <- !is.na(hgnc_vec) & hgnc_vec != "" & !duplicated(hgnc_vec)
expr_full <- v$E[keep_rows, ]; rownames(expr_full) <- hgnc_vec[keep_rows]

# =============================================================================
# 5. EXTRACT BASELINE (DX) SAMPLES
# =============================================================================

dx_cols <- grep("_DX$", colnames(expr_full), value=TRUE)
subj_dx <- trimws(sub("_DX$","", dx_cols))
baseline$subject <- trimws(baseline$subject)

idx <- match(baseline$subject, subj_dx)
if (any(is.na(idx))) {
  cat("Baseline subjects with no matching DX sample:\n"); print(baseline$subject[is.na(idx)])
}
baseline <- baseline[!is.na(idx), ]
expr_dx_full <- expr_full[, dx_cols[idx[!is.na(idx)]], drop=FALSE]
colnames(expr_dx_full) <- baseline$subject
stopifnot(identical(colnames(expr_dx_full), baseline$subject))

expr_dx <- expr_dx_full[rowMeans(expr_dx_full) >= 2, , drop=FALSE]
cat(sprintf("Patients: %d | Genes (k-TSP): %d | Genes (full): %d\n",
            ncol(expr_dx_full), nrow(expr_dx), nrow(expr_dx_full)))

missing <- setdiff(unique(c(EARLY5$gene1, EARLY5$gene2)), rownames(expr_dx_full))
if (length(missing) > 0) stop("Missing EARLY5 genes: ", paste(missing, collapse=", "))

outcome <- factor(baseline$fast_responder, levels=c("Slow","Fast"))
n <- length(outcome)

# =============================================================================
# 6. k-TSP LOOCV -- PRIMARY RESULT
# =============================================================================

loocv_score <- numeric(n); loocv_label <- character(n)
loocv_pairs <- vector("list", n)

for (i in seq_len(n)) {
  cat(sprintf("LOOCV fold %d/%d\r", i, n))
  clf   <- SWAP.Train.KTSP(expr_dx[,-i,drop=FALSE], outcome[-i], krange=1:5)
  stats <- SWAP.KTSP.Statistics(inputMat=expr_dx[,i,drop=FALSE], classifier=clf)
  loocv_score[i] <- mean(stats$comparisons)
  loocv_label[i] <- as.character(SWAP.KTSP.Classify(expr_dx[,i,drop=FALSE], clf))
  loocv_pairs[[i]] <- clf$TSPs
}
cat("\n")

loocv_df <- data.frame(subject=baseline$subject, true_label=as.character(outcome),
                       pred_label=loocv_label, ktsp_score=loocv_score)
loocv_acc <- mean(loocv_df$true_label == loocv_df$pred_label)
roc_ktsp  <- roc(loocv_df$true_label, loocv_df$ktsp_score,
                 levels=c("Slow","Fast"), direction="auto", quiet=TRUE)
auc_ktsp    <- as.numeric(auc(roc_ktsp))
auc_ktsp_ci <- ci.auc(roc_ktsp, method="delong")
cat(sprintf("\nPRIMARY: k-TSP LOOCV AUC = %.3f (95%% CI %.3f-%.3f) | Accuracy = %.3f\n",
            auc_ktsp, auc_ktsp_ci[1], auc_ktsp_ci[3], loocv_acc))

# Pairs used per fold
k_per_fold <- sapply(loocv_pairs, nrow)
cat("\nPairs selected per LOOCV fold (k):\n"); print(table(k_per_fold))
# Save pairs-per-fold counts
dir.create("output_revised/predictions/tables", recursive=TRUE, showWarnings=FALSE)
write.csv(as.data.frame(table(k_pairs = k_per_fold), responseName = "n_folds"),
          "output_revised/predictions/tables/SuppTable_pairs_per_fold.csv", row.names=FALSE)

# Out-of-fold classification summary. SWAP.KTSP.Classify uses a majority vote;
# check that this equals "score > 0.5 = Fast" (ties = Slow).
stopifnot(all((loocv_df$ktsp_score > 0.5) == (loocv_df$pred_label == "Fast")))
cm <- table(true = factor(loocv_df$true_label, levels=c("Fast","Slow")),
            pred = factor(loocv_df$pred_label, levels=c("Fast","Slow")))
print(cm)
class_df <- data.frame(
  decision_rule  = "Fast if > 50% of pair votes are for fast conversion (score > 0.5; ties = Slow)",
  n = sum(cm), true_fast = sum(cm["Fast", ]), true_slow = sum(cm["Slow", ]),
  TP = cm["Fast","Fast"], FN = cm["Fast","Slow"], FP = cm["Slow","Fast"], TN = cm["Slow","Slow"],
  sensitivity = round(cm["Fast","Fast"] / sum(cm["Fast", ]), 3),
  specificity = round(cm["Slow","Slow"] / sum(cm["Slow", ]), 3),
  PPV         = round(cm["Fast","Fast"] / sum(cm[, "Fast"]), 3),
  NPV         = round(cm["Slow","Slow"] / sum(cm[, "Slow"]), 3),
  accuracy    = round(sum(diag(cm)) / sum(cm), 3))
print(class_df, row.names=FALSE)
write.csv(class_df, "output_revised/predictions/tables/SuppTable_classification_metrics.csv", row.names=FALSE)

# Pair-selection stability across folds
all_pairs <- unlist(lapply(loocv_pairs, function(x)
  apply(x, 1, function(z) paste0(z[1]," > ",z[2]))))
pair_counts  <- sort(table(all_pairs), decreasing=TRUE)
stability_df <- data.frame(pair=names(pair_counts), selected_n=as.integer(pair_counts),
                           selected_fraction=as.integer(pair_counts)/n)
cat("\nPair-selection stability across LOOCV folds:\n"); print(head(stability_df, 20), row.names=FALSE)

# =============================================================================
# 7. BENCHMARKS -- Sweeney_OD_3 / PennNich_RISK_6
# =============================================================================

# RISK6 lists FCGR1B under its old symbol; org.Hs.eg.db returns the current
# HGNC symbol FCGR1BP. A copy of the matrix is used for signature scoring only,
# with FCGR1BP renamed to FCGR1B; the k-TSP data are unchanged.
expr_sig <- expr_dx_full
if (!"FCGR1B" %in% rownames(expr_sig) && "FCGR1BP" %in% rownames(expr_sig)) {
  rownames(expr_sig)[rownames(expr_sig) == "FCGR1BP"] <- "FCGR1B"
  cat("RISK6: FCGR1BP renamed to FCGR1B for signature scoring\n")
}

# Stop if any comparator gene is missing
for (sig_name in c("Sweeney_OD_3","PennNich_RISK_6")) {
  missing_sig_genes <- setdiff(TBsignatures[[sig_name]], rownames(expr_sig))
  if (length(missing_sig_genes) > 0)
    stop(sprintf("%s missing genes: %s", sig_name, paste(missing_sig_genes, collapse=", ")))
}

se_dx <- SummarizedExperiment(assays=list(voom_logCPM=expr_sig),
                              colData=DataFrame(subject=colnames(expr_sig)))

# update_genes = FALSE so the profiler does not rename FCGR1B back to FCGR1BP
se_dx <- runTBsigProfiler(input=se_dx, useAssay="voom_logCPM",
                          signatures=TBsignatures[c("Sweeney_OD_3","PennNich_RISK_6")],
                          algorithm="ssGSEA", combineSigAndAlgorithm=TRUE, parallel.sz=1,
                          update_genes=FALSE)

score_cols <- colnames(colData(se_dx))
s3_col <- grep("Sweeney_OD_3",    score_cols, value=TRUE)
r6_col <- grep("PennNich_RISK_6", score_cols, value=TRUE)
stopifnot(length(s3_col)==1, length(r6_col)==1)

sweeney3_score <- colData(se_dx)[[s3_col]]; risk6_score <- colData(se_dx)[[r6_col]]
names(sweeney3_score) <- names(risk6_score) <- rownames(colData(se_dx))

sw_idx <- match(colnames(expr_dx_full), names(sweeney3_score))
r6_idx <- match(colnames(expr_dx_full), names(risk6_score))
if (any(is.na(sw_idx)) || any(is.na(r6_idx))) {
  stop("Sample name mismatch between expr_dx_full and se_dx colData: ",
       paste(colnames(expr_dx_full)[is.na(sw_idx) | is.na(r6_idx)], collapse=", "))
}
sweeney3_score <- sweeney3_score[sw_idx]; risk6_score <- risk6_score[r6_idx]

roc_s3 <- roc(as.character(outcome), sweeney3_score, levels=c("Slow","Fast"), direction="auto", quiet=TRUE)
roc_r6 <- roc(as.character(outcome), risk6_score,    levels=c("Slow","Fast"), direction="auto", quiet=TRUE)
auc_s3 <- as.numeric(auc(roc_s3)); s3_ci <- ci.auc(roc_s3, method="delong")
auc_r6 <- as.numeric(auc(roc_r6)); r6_ci <- ci.auc(roc_r6, method="delong")
cat(sprintf("Sweeney_OD_3 AUC    = %.3f (95%% CI %.3f-%.3f)\n", auc_s3, s3_ci[1], s3_ci[3]))
cat(sprintf("PennNich_RISK_6 AUC = %.3f (95%% CI %.3f-%.3f)\n", auc_r6, r6_ci[1], r6_ci[3]))

# ROC direction chosen by direction = "auto": "<" = higher score in Fast,
# ">" = higher score in Slow (expected for Sweeney3 and RISK6).
cat(sprintf("ROC directions (auto): EARLY5 %s | Sweeney3 %s | RISK6 %s\n",
            roc_ktsp$direction, roc_s3$direction, roc_r6$direction))

stopifnot(roc_ktsp$direction == "<", roc_s3$direction == ">", roc_r6$direction == ">")

# Paired DeLong tests on direction-aligned curves (comparator scores negated,
# so all three curves use "<"); AUCs are unchanged.
roc_s3_al <- roc(as.character(outcome), -sweeney3_score, levels=c("Slow","Fast"), direction="<", quiet=TRUE)
roc_r6_al <- roc(as.character(outcome), -risk6_score,    levels=c("Slow","Fast"), direction="<", quiet=TRUE)
stopifnot(abs(auc(roc_s3_al) - auc_s3) < 1e-12, abs(auc(roc_r6_al) - auc_r6) < 1e-12)
delong_s3 <- roc.test(roc_ktsp, roc_s3_al, method="delong", paired=TRUE)
delong_r6 <- roc.test(roc_ktsp, roc_r6_al, method="delong", paired=TRUE)
cat(sprintf("DeLong EARLY5 vs Sweeney3: p = %.4f | EARLY5 vs RISK6: p = %.4f\n",
            delong_s3$p.value, delong_r6$p.value))

# =============================================================================
# 8. RESULTS TABLE + SAVE
# =============================================================================

results_df <- data.frame(
  Analysis = c("EARLY5 (LOOCV)", "Sweeney_OD_3", "PennNich_RISK_6"),
  n       = n,
  AUC     = round(c(auc_ktsp, auc_s3, auc_r6), 3),
  CI_low  = round(c(auc_ktsp_ci[1], s3_ci[1], r6_ci[1]), 3),
  CI_high = round(c(auc_ktsp_ci[3], s3_ci[3], r6_ci[3]), 3),
  ROC_direction = c(roc_ktsp$direction, roc_s3$direction, roc_r6$direction),
  DeLong_p_vs_EARLY5 = c(NA, signif(delong_s3$p.value, 3), signif(delong_r6$p.value, 3))
)
cat("\n--- AUC SUMMARY ---\n"); print(results_df, row.names=FALSE)

write.csv(loocv_df,     "output_revised/predictions/tables/ktsp_loocv_results.csv", row.names=FALSE)
write.csv(stability_df, "output_revised/predictions/tables/SuppTable_pair_stability_all.csv", row.names=FALSE)
write.csv(results_df,   "output_revised/predictions/tables/GSE89403_signature_AUC_summary.csv", row.names=FALSE)

# =============================================================================
# 9. FIGURE 3 -- PANELS a-e
# =============================================================================

## Panel a: cumulative culture conversion
conv_df <- as.data.frame(baseline) %>%
  filter(!is.na(timetonegativity), !is.na(treatmentresult)) %>%
  mutate(timetonegativity=factor(timetonegativity, levels=c("Week04","Week08","Week12","Week24"),
                                 labels=c("Week 4","Week 8","Week 12","Week 24")),
         treatmentresult=factor(treatmentresult, levels=c("Definite Cure","Probable Cure","Possible Cure")))

counts_df <- as.data.frame(table(conv_df$timetonegativity, conv_df$treatmentresult)) %>%
  setNames(c("timetonegativity","treatmentresult","n")) %>%
  complete(timetonegativity, treatmentresult, fill=list(n=0)) %>%
  arrange(timetonegativity, treatmentresult) %>%
  group_by(treatmentresult) %>% mutate(n_cumulative=cumsum(n)) %>% ungroup() %>%
  group_by(timetonegativity) %>%
  mutate(total_cumulative=sum(n_cumulative), pct_cumulative=total_cumulative/nrow(baseline)*100) %>%
  ungroup()

p_A <- ggplot(counts_df, aes(x=timetonegativity, y=n_cumulative, fill=treatmentresult)) +
  geom_bar(stat="identity", width=0.6, color="black", linewidth=0.4) +
  geom_text(aes(label=sprintf("n=%d\n(%.0f%%)", total_cumulative, pct_cumulative), y=total_cumulative),
            stat="unique", vjust=-0.4, size=3.5, color="black") +
  scale_fill_manual(values=c("Definite Cure"="#8FBC94","Probable Cure"="#F6E5A0","Possible Cure"="#F7C59F"),
                    name="Treatment\nResult") +
  scale_y_continuous(limits=c(0, nrow(baseline)+12), expand=expansion(mult=c(0,0))) +
  labs(x="Sputum Culture Conversion By Week", y="Cumulative N Patients") +
  pub_theme + theme(legend.position="right")
print(p_A)
ggsave("output_revised/predictions/figures/PanelA_conversion.pdf", p_A, width=6, height=5)
ggsave("output_revised/predictions/figures/PanelA_conversion.png", p_A, width=6, height=5, dpi=300)

## Panel b: EARLY5 LOOCV score by week of conversion
b_df <- data.frame(subject=baseline$subject,
                   week=factor(baseline$timetonegativity, levels=c("Week04","Week08","Week12","Week24"),
                               labels=c("Week 4","Week 8","Week 12","Week 24"))) %>%
  left_join(loocv_df[, c("subject","ktsp_score")], by="subject")

stat_b <- b_df %>%
  wilcox_test(ktsp_score ~ week, comparisons=list(c("Week 4","Week 8"), c("Week 4","Week 12"), c("Week 4","Week 24"))) %>%
  adjust_pvalue(method="BH") %>% add_significance("p.adj") %>% add_xy_position(x="week", step.increase=0.08)
stat_b$label <- fmt_p(stat_b$p.adj)

p_B <- ggplot(b_df, aes(x=week, y=ktsp_score)) +
  geom_boxplot(aes(fill=week), alpha=0.7, outlier.shape=NA, width=0.5) +
  geom_jitter(width=0.1, height=0, size=1.8, alpha=0.7) +
  stat_pvalue_manual(stat_b, label="p.adj.signif", tip.length=0.01, size=5) +
  scale_fill_manual(values=week_colors) +
  scale_y_continuous(breaks=SCORE_BREAKS) +
  labs(x="Week of culture conversion", y=SCORE_LAB) +
  pub_theme + theme(legend.position="none")
print(p_B)
ggsave("output_revised/predictions/figures/PanelB_ktsp_by_week.pdf", p_B, width=6, height=5)
ggsave("output_revised/predictions/figures/PanelB_ktsp_by_week.png", p_B, width=6, height=5, dpi=300)
cat("\nPanel b pairwise Wilcoxon on k-TSP LOOCV score (BH-adjusted):\n")
print(stat_b[, c("group1","group2","n1","n2","p","p.adj")])

## Panel c: ROC with DeLong p values
roc_df <- bind_rows(
  data.frame(fpr=1-roc_ktsp$specificities, tpr=roc_ktsp$sensitivities, Signature="k-TSP LOOCV"),
  data.frame(fpr=1-roc_s3$specificities,   tpr=roc_s3$sensitivities,   Signature="Sweeney_OD_3"),
  data.frame(fpr=1-roc_r6$specificities,   tpr=roc_r6$sensitivities,   Signature="PennNich_RISK_6")
) %>% group_by(Signature) %>% arrange(fpr, tpr, .by_group=TRUE) %>% ungroup()

sig_labels <- c(sprintf("EARLY5: %.2f (%.2f\u2013%.2f)",   auc_ktsp, auc_ktsp_ci[1], auc_ktsp_ci[3]),
                              sprintf("Sweeney3: %.2f (%.2f\u2013%.2f)", auc_s3,   s3_ci[1],       s3_ci[3]),
                              sprintf("RISK6: %.2f (%.2f\u2013%.2f)",    auc_r6,   r6_ci[1],       r6_ci[3]))
roc_df$Signature <- factor(roc_df$Signature, levels=c("k-TSP LOOCV","Sweeney_OD_3","PennNich_RISK_6"), labels=sig_labels)

p_C <- ggplot(roc_df, aes(fpr, tpr, color=Signature)) +
  geom_abline(linetype="dotted", color="grey60", linewidth=0.6) +
  geom_step(linewidth=1.1, direction="vh") +
  scale_color_manual(values=setNames(c(COL_UP, COL_DOWN, COL_R6), sig_labels)) +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0)) + scale_y_continuous(limits=c(0,1), expand=c(0.01,0)) +
  labs(x="1 - Specificity", y="Sensitivity", color=NULL) +
  pub_theme + theme(legend.position=c(0.66,0.14), legend.text=element_text(size=9),
                    legend.title=element_text(size=8),
                    legend.background=element_rect(fill="white", color="grey80"))
print(p_C)
ggsave("output_revised/predictions/figures/PanelC_ROC.pdf", p_C, width=6, height=5)
ggsave("output_revised/predictions/figures/PanelC_ROC.png", p_C, width=6, height=5, dpi=300)

## Panel d: gene-pair selection frequency (top 5 shown; full list in supplement)
top5_stability <- head(stability_df, 5)
p_D <- ggplot(top5_stability, aes(x=selected_fraction, y=reorder(pair, selected_fraction))) +
  geom_segment(aes(xend=0, yend=pair), color="grey70", linewidth=0.8) +
  geom_point(size=4.5, color=COL_UP) +
  geom_text(aes(label=sprintf("%d/%d", selected_n, n)), hjust=-0.3, size=3.5) +
  scale_x_continuous(limits=c(0, 1.15), breaks=seq(0, 1, 0.25), labels=scales::percent) +
  labs(x="Selection frequency across LOOCV folds", y=NULL) +
  pub_theme
print(p_D)
ggsave("output_revised/predictions/figures/PanelD_stability.pdf", p_D, width=6, height=4)
ggsave("output_revised/predictions/figures/PanelD_stability.png", p_D, width=6, height=4, dpi=300)

## Panel e: bacterial burden vs EARLY5 LOOCV score (MGIT; Xpert in Supplementary Figure S2)
burden_df <- loocv_df %>% left_join(as.data.frame(baseline)[, c("subject","mgit","xpert")], by="subject")

spearman_row <- function(df, xvar) {
  d  <- df[!is.na(df[[xvar]]), ]
  ct <- suppressWarnings(cor.test(d$ktsp_score, d[[xvar]], method="spearman", exact=FALSE))
  data.frame(Variable=xvar, n=nrow(d), rho=round(ct$estimate, 3), p=signif(ct$p.value, 3))
}
burden_cor <- bind_rows(spearman_row(burden_df, "mgit"), spearman_row(burden_df, "xpert"))
# Xpert Spearman excluding Ct >= 40 (sensitivity analysis)
xpert_lt40 <- spearman_row(burden_df[is.na(burden_df$xpert) | burden_df$xpert < 40, ], "xpert")
xpert_lt40$Variable <- "xpert (Ct < 40 only)"
burden_cor <- bind_rows(burden_cor, xpert_lt40)
# Number of patients at the Xpert Ct ceiling (40)
cat(sprintf("Xpert Ct = 40 (possible 'not detected' coding): %d of %d patients\n",
            sum(burden_df$xpert >= 40, na.rm=TRUE), sum(!is.na(burden_df$xpert))))
cat("\n--- Spearman: EARLY5 LOOCV score vs baseline burden ---\n")
print(burden_cor, row.names=FALSE)
write.csv(burden_cor, "output_revised/predictions/tables/burden_spearman.csv", row.names=FALSE)

burden_panel <- function(df, xvar, xlab, col, cor_row) {
  d <- df[!is.na(df[[xvar]]), ]
  lab <- sprintf("rho = %.2f\n%s\nn = %d", cor_row$rho, fmt_p(cor_row$p), cor_row$n)
  ggplot(d, aes(x=.data[[xvar]], y=ktsp_score)) +
    geom_count(color=col, alpha=0.8) +
    scale_size_area(max_size=7, name="Patients",
                    breaks=function(l) seq(1, max(l))) +
    annotate("text", x=min(d[[xvar]]), y=1.12, hjust=0, vjust=1, size=3.5, label=lab) +
    scale_y_continuous(breaks=SCORE_BREAKS, limits=c(-0.02, 1.14)) +
    labs(x=xlab, y=SCORE_LAB) +
    pub_theme + theme(legend.position="right", legend.title=element_text(size=9))
}

p_E <- burden_panel(burden_df, "mgit",  "Baseline MGIT time-to-positivity (days)", COL_DOWN,
                    burden_cor[burden_cor$Variable=="mgit", ])
p_S_xpert <- burden_panel(burden_df, "xpert", "Baseline Xpert Ct", COL_R6,
                          burden_cor[burden_cor$Variable=="xpert", ])
print(p_E); print(p_S_xpert)
ggsave("output_revised/predictions/figures/PanelE_MGIT.pdf", p_E, width=6, height=5)
ggsave("output_revised/predictions/figures/PanelE_MGIT.png", p_E, width=6, height=5, dpi=300)
ggsave("output_revised/predictions/figures/SuppFig_Xpert.pdf", p_S_xpert, width=6, height=5)
ggsave("output_revised/predictions/figures/SuppFig_Xpert.png", p_S_xpert, width=6, height=5, dpi=300)

# =============================================================================
# 10. TABLE 2 -- week-4 conversion ~ EARLY5 + MGIT TTP (logistic, complete cases)
#     Score per 0.2 (one pair vote), profile-likelihood 95% CIs,
#     likelihood-ratio p values, Firth sensitivity analysis
# =============================================================================

expl_df <- burden_df %>%
  mutate(fast_binary=as.integer(true_label=="Fast"),
         score02=ktsp_score/0.2,
         mgit_iqr=mgit/IQR(mgit, na.rm=TRUE)) %>%
  filter(!is.na(mgit_iqr))
cat(sprintf("\nTable 2 models: n = %d MGIT-complete patients, %d fast converters\n",
            nrow(expl_df), sum(expl_df$fast_binary)))

fit_k_only <- glm(fast_binary ~ score02,            data=expl_df, family="binomial")
fit_m_only <- glm(fast_binary ~ mgit_iqr,           data=expl_df, family="binomial")
fit_k_m    <- glm(fast_binary ~ score02 + mgit_iqr, data=expl_df, family="binomial")
lrt_k_m    <- anova(fit_m_only, fit_k_m, test="Chisq")

# Likelihood-ratio p value for each term (drop1, test = "LRT");
# the Wald p value is kept in a separate column.
or_row <- function(fit, term, model, label) {
  ci <- suppressMessages(confint(fit))[term, ]
  lr <- drop1(fit, test="LRT")
  data.frame(Model=model, Predictor=label, n=nobs(fit),
             OR=round(exp(coef(fit)[term]), 2),
             CI_low=round(exp(ci[1]), 2), CI_high=round(exp(ci[2]), 2),
             p=signif(lr[term, "Pr(>Chi)"], 3),
             p_wald=signif(summary(fit)$coefficients[term, "Pr(>|z|)"], 3), row.names=NULL)
}
table2 <- bind_rows(
  or_row(fit_k_only, "score02",  "EARLY5 only",       "EARLY5 score (per 0.2)"),
  or_row(fit_m_only, "mgit_iqr", "MGIT TTP only",     "MGIT TTP (per IQR)"),
  or_row(fit_k_m,    "score02",  "EARLY5 + MGIT TTP", "EARLY5 score (per 0.2)"),
  or_row(fit_k_m,    "mgit_iqr", "EARLY5 + MGIT TTP", "MGIT TTP (per IQR)"),
  data.frame(Model="LRT: EARLY5 + MGIT vs MGIT only", Predictor="", n=nobs(fit_k_m),
             OR=NA, CI_low=NA, CI_high=NA, p=signif(lrt_k_m$`Pr(>Chi)`[2], 3), p_wald=NA)
)
cat("\n--- TABLE 2: week-4 culture conversion (logistic, MGIT-complete cases) ---\n")
print(table2, row.names=FALSE)
write.csv(table2,    "output_revised/predictions/tables/Table2_logistic_EARLY5_MGIT.csv", row.names=FALSE)
write.csv(burden_df, "output_revised/predictions/tables/supp_burden_scores.csv", row.names=FALSE)

# Firth penalised logistic regression (sensitivity analysis)
fit_firth <- logistf(fast_binary ~ score02 + mgit_iqr, data=expl_df)
firth_tab <- data.frame(Predictor=names(fit_firth$coefficients), n=nrow(expl_df),
                        OR=round(exp(fit_firth$coefficients), 2),
                        CI_low=round(exp(fit_firth$ci.lower), 2),
                        CI_high=round(exp(fit_firth$ci.upper), 2),
                        p=signif(fit_firth$prob, 3), row.names=NULL) %>%
  filter(Predictor != "(Intercept)")
cat("\n--- SUPP: Firth logistic (EARLY5 per 0.2 + MGIT TTP) ---\n"); print(firth_tab, row.names=FALSE)
write.csv(firth_tab, "output_revised/predictions/tables/SuppTable_firth_logistic.csv", row.names=FALSE)

# =============================================================================
# 11. SUPPLEMENTARY -- time to conversion by out-of-fold EARLY5 call (KM + Cox)
#  - Unadjusted and MGIT-adjusted HRs; continuous score per 0.2; PH checks.
#  - Grouping uses loocv_label (majority vote, out-of-fold), no tuned threshold.
#  - Event = culture conversion; all 76 convert by week 24 (no censoring).
#  - Conversion observed only at weeks 4/8/12/24; Efron ties.
# =============================================================================

km_df <- loocv_df %>%
  left_join(as.data.frame(baseline)[, c("subject","timetonegativity","mgit")], by="subject") %>%
  mutate(week_num = as.numeric(sub("Week0?", "", timetonegativity)),
         event    = 1L,
         group    = factor(ifelse(pred_label=="Fast","Predicted fast","Predicted slow"),
                           levels=c("Predicted slow","Predicted fast")),
         score02  = ktsp_score/0.2,
         mgit_iqr = mgit/IQR(mgit, na.rm=TRUE))
stopifnot(!any(is.na(km_df$week_num)))

km_fit  <- survfit(Surv(week_num, event) ~ group, data=km_df)
km_sum  <- summary(km_fit, times=c(0, 4, 8, 12, 24), extend=TRUE)
grp_n   <- table(km_df$group)
km_plot <- data.frame(time=km_sum$time, converted=1 - km_sum$surv, n_risk=km_sum$n.risk,
                      grp=factor(sub("group=", "", as.character(km_sum$strata)), levels=names(grp_n)))
km_plot$group <- factor(km_plot$grp, levels=names(grp_n),
                        labels=sprintf("%s (n=%d)", names(grp_n), as.integer(grp_n)))

lr   <- survdiff(Surv(week_num, event) ~ group, data=km_df)
lr_p <- pchisq(lr$chisq, df=length(lr$n) - 1, lower.tail=FALSE)

cc_km    <- km_df %>% filter(!is.na(mgit))
cox_u    <- coxph(Surv(week_num, event) ~ group,              data=km_df, ties="efron")
cox_u_cc <- coxph(Surv(week_num, event) ~ group,              data=cc_km, ties="efron")
cox_adj  <- coxph(Surv(week_num, event) ~ group + mgit_iqr,   data=cc_km, ties="efron")
cox_cont <- coxph(Surv(week_num, event) ~ score02 + mgit_iqr, data=cc_km, ties="efron")

hr_row <- function(fit, term, model, n_used) {
  s <- summary(fit)
  data.frame(Model=model, Term=term, n=n_used,
             HR=round(s$conf.int[term, "exp(coef)"], 2),
             CI_low=round(s$conf.int[term, "lower .95"], 2),
             CI_high=round(s$conf.int[term, "upper .95"], 2),
             p=signif(s$coefficients[term, "Pr(>|z|)"], 3))
}
cox_table <- bind_rows(
  hr_row(cox_u,    "groupPredicted fast", "Unadjusted",                        nrow(km_df)),
  hr_row(cox_u_cc, "groupPredicted fast", "Unadjusted (MGIT complete cases)",  nrow(cc_km)),
  hr_row(cox_adj,  "groupPredicted fast", "Adjusted for MGIT TTP",             nrow(cc_km)),
  hr_row(cox_adj,  "mgit_iqr",            "Adjusted for MGIT TTP",             nrow(cc_km)),
  hr_row(cox_cont, "score02",             "Continuous score (per 0.2) + MGIT TTP", nrow(cc_km)),
  hr_row(cox_cont, "mgit_iqr",            "Continuous score (per 0.2) + MGIT TTP", nrow(cc_km)))
cat(sprintf("\n--- SUPP KM: log-rank p = %s ---\n", signif(lr_p, 3)))
cat("Cox models (HR > 1 = faster conversion):\n"); print(cox_table, row.names=FALSE)
write.csv(cox_table, "output_revised/predictions/tables/SuppTable_cox_models.csv", row.names=FALSE)

# Cox models with exact handling of ties (sensitivity analysis)
cox_exact <- bind_rows(
  hr_row(coxph(Surv(week_num, event) ~ group,              data=km_df, ties="exact"), "groupPredicted fast", "Unadjusted (exact ties)",              nrow(km_df)),
  hr_row(coxph(Surv(week_num, event) ~ group + mgit_iqr,   data=cc_km, ties="exact"), "groupPredicted fast", "Adjusted for MGIT TTP (exact ties)",   nrow(cc_km)),
  hr_row(coxph(Surv(week_num, event) ~ score02 + mgit_iqr, data=cc_km, ties="exact"), "score02",             "Continuous score + MGIT TTP (exact ties)", nrow(cc_km)))
cat("Cox models with exact ties:\n"); print(cox_exact, row.names=FALSE)
write.csv(cox_exact, "output_revised/predictions/tables/SuppTable_cox_models_exact_ties.csv", row.names=FALSE)

cat("\nProportional-hazards checks (cox.zph):\n")
print(cox.zph(cox_u)); print(cox.zph(cox_adj)); print(cox.zph(cox_cont))

unadj_row <- cox_table[cox_table$Model=="Unadjusted", ]
adj_row   <- cox_table[cox_table$Model=="Adjusted for MGIT TTP" & cox_table$Term=="groupPredicted fast", ]
km_subtitle <- sprintf("Log-rank %s\nUnadjusted HR %.2f (%.2f-%.2f), n = %d\nMGIT-adjusted HR %.2f (%.2f-%.2f), %s, n = %d",
                       fmt_p(lr_p),
                       unadj_row$HR, unadj_row$CI_low, unadj_row$CI_high, unadj_row$n,
                       adj_row$HR, adj_row$CI_low, adj_row$CI_high, fmt_p(adj_row$p), adj_row$n)

km_cols <- setNames(c(COL_SLOW, COL_FAST), levels(km_plot$group))
p_KM <- ggplot(km_plot, aes(x=time, y=converted, color=group)) +
  geom_step(linewidth=1.1) +
  geom_point(data=subset(km_plot, time > 0), size=2.5) +
  scale_color_manual(values=km_cols) +
  scale_x_continuous(breaks=c(0, 4, 8, 12, 24), limits=c(0, 24.5)) +
  scale_y_continuous(limits=c(0, 1), breaks=SCORE_BREAKS, labels=scales::percent) +
  labs(x="Weeks on treatment", y="Culture converted", color="EARLY5 (out-of-fold)",
       subtitle=km_subtitle) +                                  # stats above panel, no overlap
  pub_theme + theme(legend.position=c(0.78, 0.22),              # empty lower-right region
                    legend.title=element_text(size=9), legend.text=element_text(size=9),
                    plot.subtitle=element_text(size=9),
                    legend.background=element_rect(fill="white", color="grey80"))

p_risk <- ggplot(km_plot, aes(x=time, y=grp, label=n_risk, color=grp)) +
  geom_text(size=3.5) +
  scale_color_manual(values=setNames(c(COL_SLOW, COL_FAST), names(grp_n)), guide="none") +
  scale_x_continuous(breaks=c(0, 4, 8, 12, 24), limits=c(0, 24.5)) +
  labs(x=NULL, y=NULL, title="Number at risk") +
  theme_classic(base_size=11) +
  theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(),
        plot.title=element_text(size=10))

supp_km <- p_KM / p_risk + plot_layout(heights=c(4, 1))
print(supp_km)
ggsave("output_revised/predictions/figures/SuppFig_KM_conversion.pdf", supp_km, width=7, height=7)
ggsave("output_revised/predictions/figures/SuppFig_KM_conversion.png", supp_km, width=7, height=7, dpi=300)

# =============================================================================
# 12. ASSEMBLED FIGURE 3 (a-e)
# =============================================================================

fig3 <- (p_A | p_B) / (p_C | (p_D / p_E)) +
  plot_layout(heights=c(1, 1.5)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=20))
ggsave("output_revised/predictions/figures/Figure3_assembled.pdf", fig3, width=14, height=14)
ggsave("output_revised/predictions/figures/Figure3_assembled.png", fig3, width=14, height=14, dpi=300)

cat("\nDone.\n")