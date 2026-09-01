suppressPackageStartupMessages({
  library(switchBox); library(pROC); library(limma); library(edgeR)
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(ggpubr); library(biomaRt); library(rstatix); library(ggpubr)
})
set.seed(42)

COL_UP       <- "#C1121F"
COL_DOWN     <- "#1D4E89"
COL_BASELINE <- "#D4A96A"
COL_WEEK4    <- "tomato3"

dir.create("output/meta/tables", recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# DATA PREP
# =============================================================================

metadata <- read.table("metadata/gse89403/SraRunTable.csv", sep=",", header=TRUE)
metadata$tech_group <- paste(metadata$subject, metadata$time, sep="_")

clinical <- metadata %>%
  group_by(tech_group) %>%
  summarise(subject=subject[1], timepoint=time[1],
            timetonegativity=timetonegativity[1],
            treatmentresult=treatmentresult[1],
            mgit=mean(mgit, na.rm=TRUE),
            xpert=mean(xpert, na.rm=TRUE), .groups="drop")

baseline <- clinical %>%
  filter(timepoint=="DX", !is.na(timetonegativity)) %>%
  mutate(fast_responder=ifelse(timetonegativity=="Week04","Fast","Slow"))

# =============================================================================
# VOOM NORMALIZATION
# =============================================================================

count_files <- list.files("data/mrna/gse89403/counts", full.names=TRUE)
dat <- lapply(count_files, function(f) {
  df <- fread(f); df_sub <- df[, c(1,7)]
  colnames(df_sub) <- c("GeneID", gsub(".*_([0-9]+)_.*","\\1",colnames(df)[7]))
  df_sub
})
mat <- do.call(cbind, lapply(dat, function(x) x[,-1]))
mat <- as.data.frame(mat)[, order(as.numeric(colnames(mat)))]
rownames(mat) <- dat[[1]]$GeneID
colnames(mat) <- gsub(".*/(SRR[0-9]+)\\.bam","\\1",colnames(mat))

meta_sub   <- metadata[, c("Run","subject","time")]
has_both   <- with(meta_sub, table(subject,time)) >= 1
valid_subj <- rownames(has_both[rowSums(has_both)==2, ])
meta_sub   <- meta_sub[meta_sub$subject %in% valid_subj, ]
meta_sub$tech_group <- paste(meta_sub$subject, meta_sub$time, sep="_")

mat_filtered  <- mat[, meta_sub$Run]
summed_counts <- t(rowsum(t(mat_filtered), group=meta_sub$tech_group))

sample_info <- data.frame(sample_name=colnames(summed_counts)) %>%
  separate(sample_name, into=c("subject","timepoint"), sep="_", extra="merge")
sample_info$sample_name <- colnames(summed_counts)
sample_info$timepoint   <- factor(sample_info$timepoint, levels=c("DX","week_4"))

keep   <- rowSums(cpm(summed_counts) > 1) >= 10
dge    <- DGEList(counts=summed_counts[keep, ])
dge    <- calcNormFactors(dge)
design <- model.matrix(~subject + timepoint, data=sample_info)
v      <- voom(dge, design, plot=FALSE)

# =============================================================================
# GENE SYMBOL MAPPING
# =============================================================================

mart    <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl",
                      mirror="useast")
ens_ids <- sub("\\..*","", rownames(v$E))
ens_map <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="ensembl_gene_id", values=unique(ens_ids), mart=mart)

hgnc_vec  <- ens_map$hgnc_symbol[match(ens_ids, ens_map$ensembl_gene_id)]
keep_rows <- !is.na(hgnc_vec) & hgnc_vec != "" & !duplicated(hgnc_vec)
expr_full <- v$E[keep_rows, ]
rownames(expr_full) <- hgnc_vec[keep_rows]

dx_cols <- grep("_DX$", colnames(expr_full), value=TRUE)
subj_dx <- sub("_DX$","", dx_cols)
idx     <- match(baseline$subject, subj_dx)
valid   <- !is.na(idx)
baseline <- baseline[valid, ]

expr_dx_full <- expr_full[, dx_cols][, idx[valid]]
colnames(expr_dx_full) <- baseline$subject

expr_dx <- expr_dx_full[rowMeans(expr_dx_full) >= 2, ]
cat("Genes retained after mean log-CPM >= 2 filter on DX samples:", nrow(expr_dx), "\n")

outcome <- factor(baseline$fast_responder, levels=c("Slow","Fast"))
n       <- ncol(expr_dx)
cat("Class distribution:\n"); print(table(outcome))
cat("Majority-class ceiling:", round(max(table(outcome))/n, 3), "\n\n")

# =============================================================================
# k-TSP LOOCV
# =============================================================================

loocv_score <- numeric(n)
loocv_label <- character(n)

for (i in seq_len(n)) {
  clf   <- SWAP.Train.KTSP(expr_dx[,-i,drop=FALSE], outcome[-i], krange=1:5)
  stats <- SWAP.KTSP.Statistics(inputMat=expr_dx[,i,drop=FALSE], classifier=clf)
  loocv_score[i] <- mean(stats$comparisons)
  loocv_label[i] <- as.character(SWAP.KTSP.Classify(expr_dx[,i,drop=FALSE], clf))
}

loocv_df <- data.frame(subject=baseline$subject, true_label=as.character(outcome),
                       pred_label=loocv_label, ktsp_score=loocv_score,
                       fast_responder=baseline$fast_responder)

loocv_acc <- mean(loocv_df$true_label == loocv_df$pred_label)
roc_ktsp  <- roc(loocv_df$true_label, loocv_df$ktsp_score,
                 levels=c("Slow","Fast"), direction="<", quiet=TRUE)
auc_ci    <- ci.auc(roc_ktsp, method="delong")
obs_auc   <- as.numeric(auc(roc_ktsp))

final_clf <- SWAP.Train.KTSP(expr_dx, outcome, krange=1:5)
tsp_mat   <- final_clf$TSPs

cat("k-TSP pairs (final model):\n"); print(tsp_mat)
cat(sprintf("LOOCV AUC: %.3f (95%% CI %.3f-%.3f) | Accuracy: %.3f\n",
            obs_auc, auc_ci[1], auc_ci[3], loocv_acc))

# =============================================================================
# CLINICAL DATA + LOGISTIC REGRESSION
# =============================================================================

score_clin <- loocv_df %>%
  left_join(as.data.frame(baseline)[, c("subject","mgit","xpert")], by="subject") %>%
  mutate(fast_binary=as.integer(fast_responder=="Fast"),
         mgit_iqr=mgit / IQR(mgit, na.rm=TRUE))

fit_ktsp_only <- glm(fast_binary ~ ktsp_score,             data=score_clin, family="binomial")
fit_mgit_only <- glm(fast_binary ~ mgit_iqr,               data=score_clin, family="binomial")
fit_ktsp_mgit <- glm(fast_binary ~ ktsp_score + mgit_iqr,  data=score_clin, family="binomial")
rt_result    <- anova(fit_mgit_only, fit_ktsp_mgit, test="Chisq")

cat("\n--- k-TSP only ---\n");           print(summary(fit_ktsp_only)$coefficients)
cat("\n--- MGIT TTD (per IQR) only ---\n"); print(summary(fit_mgit_only)$coefficients)
cat("\n--- k-TSP + MGIT TTD ---\n");     print(summary(fit_ktsp_mgit)$coefficients)
cat("\n--- LRT ---\n");                  print(lrt_result)
cat("\n--- ORs 95% CI ---\n")
print(round(exp(cbind(OR=coef(fit_ktsp_mgit), confint(fit_ktsp_mgit))), 3))

round(exp(cbind(OR=coef(fit_ktsp_only), confint(fit_ktsp_only))), 3)
round(exp(cbind(OR=coef(fit_mgit_only), confint(fit_mgit_only))), 3)

cor_mgit <- cor.test(score_clin$ktsp_score, score_clin$mgit, method="spearman")
cat(sprintf("Spearman k-TSP ~ MGIT TTD: rho=%.3f, p=%s\n",
            cor_mgit$estimate, signif(cor_mgit$p.value, 3)))

# =============================================================================
# THEME
# =============================================================================

pub_theme <- theme_classic(base_size=16) +
  theme(axis.title=element_text(size=16, color="black"),
        axis.text=element_text(size=14, color="black"),
        plot.title=element_text(size=17, hjust=0.5, color="black"),
        plot.subtitle=element_text(size=13, hjust=0.5, color="grey40"),
        legend.title=element_text(size=14, color="black"),
        legend.text=element_text(size=13, color="black"),
        axis.line=element_line(color="black", linewidth=0.6),
        axis.ticks=element_line(color="black", linewidth=0.5),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.6),
        panel.grid.major=element_line(color="grey92", linewidth=0.3),
        panel.grid.minor=element_blank(),
        plot.margin=unit(c(10,15,10,10), "pt"))

# =============================================================================
# FIGURES
# =============================================================================

# Culture conversion timing - CUMULATIVE STACKED
conv_df <- as.data.frame(baseline) %>%
  filter(!is.na(timetonegativity), !is.na(treatmentresult)) %>%
  mutate(
    timetonegativity = factor(timetonegativity,
                              levels = c("Week04","Week08","Week12","Week24"),
                              labels = c("Week 4","Week 8","Week 12","Week 24")),
    treatmentresult = factor(treatmentresult,
                             levels = c("Definite Cure","Probable Cure","Possible Cure"))
  )

N_total <- 76

# Build cumulative counts per cure type
counts_df <- as.data.frame(table(conv_df$timetonegativity, conv_df$treatmentresult)) %>%
  setNames(c("timetonegativity","treatmentresult","n")) %>%
  # complete so every timepoint x cure type combo exists (even if 0)
  complete(timetonegativity, treatmentresult, fill = list(n = 0)) %>%
  arrange(timetonegativity, treatmentresult) %>%
  group_by(treatmentresult) %>%
  mutate(n_cumulative = cumsum(n)) %>%
  ungroup() %>%
  # cumulative total across all cure types at each timepoint (for label)
  group_by(timetonegativity) %>%
  mutate(
    total_cumulative = sum(n_cumulative),
    pct_cumulative   = total_cumulative / N_total * 100
  ) %>%
  ungroup()

p_conv <- ggplot(counts_df,
                 aes(x = timetonegativity, y = n_cumulative, fill = treatmentresult)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
  # top label: cumulative n and %
  geom_text(aes(label = sprintf("n = %d\n(%.1f%%)", total_cumulative, pct_cumulative),
                y = total_cumulative),
            stat = "unique", vjust = -0.4, size = 4, color = "black") +
  scale_fill_manual(
    values = c("Definite Cure" = "#8FBC94",
               "Probable Cure" = "#F6E5A0",
               "Possible Cure" = "#F7C59F"),
    name = "Treatment Result"
  ) +
  scale_y_continuous(
    limits = c(0, 95),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Sputum Culture Conversion By Week",
       y = "Cumulative Number of Patients") +
  pub_theme +
  theme(legend.position = "right")

print(p_conv)

# COLORS

week_colors <- c(
  "Week 4"  = "#D95F5F",
  "Week 8"  = "#E8A87C",
  "Week 12" = "#7FBFCF",
  "Week 24" = "#9B72CF"
)

# Supp: clinical covariates (MGIT TTD + Xpert Ct only)
clin_long <- score_clin %>%
  dplyr::select(subject, fast_responder, mgit, xpert) %>%
  pivot_longer(cols=c(mgit, xpert), names_to="variable", values_to="value") %>%
  mutate(variable=factor(variable, levels=c("mgit","xpert"),
                         labels=c("MGIT TTD (days)","Xpert Ct")))

p_clin <- ggplot(clin_long, aes(x=fast_responder, y=value, fill=fast_responder)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA, width=0.5, linewidth=0.6, color="black") +
  scale_fill_manual(values=c("Fast"=COL_WEEK4, "Slow"=COL_BASELINE)) +
  scale_x_discrete(labels=c("Fast"="Fast","Slow"="Slow")) +
  facet_wrap(~variable, scales="free_y", nrow=1) +
  stat_compare_means(method="wilcox.test", label="p.format", size=4.5, color="black") +
  labs(x=NULL, y="Value") +
  pub_theme +
  theme(legend.position="none",
        strip.text=element_text(size=14, color="black"),
        strip.background=element_rect(fill="grey95", color="black", linewidth=0.4))
print(p_clin)

# Supp: k-TSP vs MGIT TTD
p_mgit <- ggplot(score_clin, aes(x=mgit, y=ktsp_score, color=fast_responder)) +
  geom_point(size=3.5, alpha=0.85) +
  geom_smooth(method="lm", se=TRUE, color="grey30", fill="grey85", linewidth=0.8) +
  scale_color_manual(values=c("Fast"=COL_WEEK4, "Slow"=COL_BASELINE), name="Converter") +
  annotate("text", x=min(score_clin$mgit, na.rm=TRUE), y=max(score_clin$ktsp_score)*0.97,
           hjust=0, size=5, color="black",
           label=sprintf("rho = %.3f\np = %s", cor_mgit$estimate,
                         signif(cor_mgit$p.value, 2))) +
  labs(x="Baseline MGIT TTD (days)", y="k-TSP Score") +
  pub_theme +
  theme(legend.position=c(0.15,0.85),
        legend.background=element_rect(fill="white", color="grey80", linewidth=0.4))
print(p_mgit)

# Save figs and tables

write.csv(loocv_df, "output/meta/tables/ktsp_loocv_results.csv", row.names=FALSE)
write.csv(pairs_df, "output/meta/tables/ktsp_pairs.csv",         row.names=FALSE)

lor <- round(exp(cbind(OR=coef(fit_ktsp_mgit), confint(fit_ktsp_mgit))), 3)
write.csv(data.frame(
  Model=c("k-TSP only","k-TSP only","MGIT TTD only","MGIT TTD only",
          "k-TSP+MGIT TTD","k-TSP+MGIT TTD","k-TSP+MGIT TTD","LRT"),
  Predictor=c("Intercept","k-TSP","Intercept","MGIT TTD (per IQR)",
              "Intercept","k-TSP","MGIT TTD (per IQR)","k-TSP+MGIT TTD vs MGIT TTD"),
  OR=c(NA, round(exp(coef(fit_ktsp_only)["ktsp_score"]),3),
       NA, round(exp(coef(fit_mgit_only)["mgit_iqr"]),3),
       NA, lor["ktsp_score","OR"], lor["mgit_iqr","OR"], NA),
  p_value=c(NA, summary(fit_ktsp_only)$coefficients["ktsp_score","Pr(>|z|)"],
            NA, summary(fit_mgit_only)$coefficients["mgit_iqr","Pr(>|z|)"],
            NA, summary(fit_ktsp_mgit)$coefficients["ktsp_score","Pr(>|z|)"],
            summary(fit_ktsp_mgit)$coefficients["mgit_iqr","Pr(>|z|)"],
            lrt_result$`Pr(>Chi)`[2])),
  "output/meta/tables/ktsp_logistic_regression.csv", row.names=FALSE)
