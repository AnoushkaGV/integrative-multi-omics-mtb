# =============================================================================
# Immune deconvolution of GSE89403 (Figure 2d, Supplementary Table S8)
# CIBERSORTx (LM22) fractions and MCP-counter scores; paired Wilcoxon tests
# between baseline and week 4 with Benjamini-Hochberg correction.
# =============================================================================
library(ggplot2)
library(patchwork)
library(MCPcounter)

out_dir <- "output_revised/deconvolution/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- "output_revised/deconvolution/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

COL_BASELINE <- "#D4A96A"
COL_WEEK4    <- "tomato3"

cib_gse  <- read.csv("data/mrna/gse89403/cibersort/CIBERSORTx_Job8_Results.csv", check.names=FALSE)
tpm_gse2 <- read.delim("data/mrna/gse89403/cibersort/cibersortx_input_final.tsv", check.names=FALSE)

# Keep samples with CIBERSORTx p < 0.05
cat(sprintf("CIBERSORTx samples with p < 0.05: %d of %d\n",
            sum(cib_gse$`P-value` < 0.05), nrow(cib_gse)))
cib_gse <- cib_gse[cib_gse$`P-value` < 0.05, ]

cib_gse$Timepoint <- ifelse(grepl("DX", cib_gse$Mixture), "Baseline", "Week4")
cib_gse$PairID    <- sub("_(DX|week_4)$", "", cib_gse$Mixture)
cib_gse_paired    <- cib_gse[cib_gse$PairID %in% intersect(
  cib_gse$PairID[cib_gse$Timepoint=="Baseline"],
  cib_gse$PairID[cib_gse$Timepoint=="Week4"]), ]

mcp_mat           <- as.matrix(tpm_gse2[, -1])
rownames(mcp_mat) <- tpm_gse2[[1]]
mcp_res           <- MCPcounter.estimate(mcp_mat, featuresType="HUGO_symbols")
mcp_gse_df           <- as.data.frame(t(mcp_res))
mcp_gse_df$Mixture   <- rownames(mcp_gse_df)
mcp_gse_df$Timepoint <- ifelse(grepl("DX", mcp_gse_df$Mixture), "Baseline", "Week4")
mcp_gse_df$PairID    <- sub("_(DX|week_4)$", "", mcp_gse_df$Mixture)
mcp_gse_df           <- mcp_gse_df[mcp_gse_df$PairID %in% intersect(
  mcp_gse_df$PairID[mcp_gse_df$Timepoint=="Baseline"],
  mcp_gse_df$PairID[mcp_gse_df$Timepoint=="Week4"]), ]

cell_types <- colnames(cib_gse)[2:23]

cib_labels <- c(
  "B cells naive"                = "B naive",
  "B cells memory"               = "B memory",
  "Plasma cells"                 = "Plasma",
  "T cells CD8"                  = "CD8+ T",
  "T cells CD4 naive"            = "CD4+ naive",
  "T cells CD4 memory resting"   = "CD4+ mem rest",
  "T cells CD4 memory activated" = "CD4+ mem act",
  "T cells follicular helper"    = "Tfh",
  "T cells regulatory (Tregs)"   = "Tregs",
  "T cells gamma delta"          = "Tgd",
  "NK cells resting"             = "NK rest",
  "NK cells activated"           = "NK act",
  "Monocytes"                    = "Mono",
  "Macrophages M0"               = "Mac M0",
  "Macrophages M1"               = "Mac M1",
  "Macrophages M2"               = "Mac M2",
  "Dendritic cells resting"      = "DC rest",
  "Dendritic cells activated"    = "DC act",
  "Mast cells resting"           = "Mast rest",
  "Mast cells activated"         = "Mast act",
  "Eosinophils"                  = "Eos",
  "Neutrophils"                  = "Neutrophils"
)

mcp_labels <- c(
  "T cells"                 = "T cells",
  "CD8 T cells"             = "CD8+ T",
  "Cytotoxic lymphocytes"   = "Cytotoxic",
  "B lineage"               = "B lineage",
  "NK cells"                = "NK",
  "Monocytic lineage"       = "Monocytic",
  "Myeloid dendritic cells" = "mDC",
  "Neutrophils"             = "Neutrophils"
)

sig_label <- function(fdr) {
  ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "")))
}

run_wilcox_paired <- function(df, cts) {
  bl  <- df[df$Timepoint=="Baseline",]; bl  <- bl[order(bl$PairID),]
  wk4 <- df[df$Timepoint=="Week4",];    wk4 <- wk4[order(wk4$PairID),]
  res <- data.frame(CellType=cts, p=NA, mean_bl=NA, mean_wk4=NA)
  for (ct in cts) {
    res$p[res$CellType==ct]        <- wilcox.test(bl[[ct]], wk4[[ct]], paired=TRUE)$p.value
    res$mean_bl[res$CellType==ct]  <- mean(bl[[ct]])
    res$mean_wk4[res$CellType==ct] <- mean(wk4[[ct]])
  }
  res$fdr <- p.adjust(res$p, method="BH")
  res
}

stats_gse <- run_wilcox_paired(cib_gse_paired, cell_types)

keep_cts <- stats_gse$CellType[
  stats_gse$fdr < 0.05 &
    (stats_gse$mean_bl > 0.01 | stats_gse$mean_wk4 > 0.01)
]

make_cib_plot <- function(df, stats, keep, title) {
  long <- do.call(rbind, lapply(keep, function(ct) {
    data.frame(CellType=cib_labels[ct], Proportion=df[[ct]],
               Timepoint=df$Timepoint, PairID=df$PairID)
  }))
  labs_ordered   <- cib_labels[keep]
  long$CellType  <- factor(long$CellType, levels=labs_ordered)
  long$Timepoint <- factor(long$Timepoint, levels=c("Baseline","Week4"))
  
  ymax   <- tapply(long$Proportion, long$CellType, max)
  sig_df <- data.frame(
    CellType = factor(labs_ordered, levels=labs_ordered),
    label    = sig_label(stats$fdr[match(keep, stats$CellType)]),
    y        = ymax[labs_ordered] * 1.12
  )
  
  ggplot(long, aes(x=Timepoint, y=Proportion)) +
    geom_line(aes(group=PairID), color="grey70", linewidth=0.2, alpha=0.3) +
    geom_boxplot(aes(fill=Timepoint), outlier.shape=NA, alpha=0.85,
                 width=0.5, color="grey30", linewidth=0.4) +
    geom_text(data=sig_df, aes(x=1.5, y=y, label=label),
              inherit.aes=FALSE, size=5) +
    scale_fill_manual(values=c(Baseline=COL_BASELINE, Week4=COL_WEEK4)) +
    scale_y_continuous(expand=expansion(mult=c(0.05, 0.25))) +
    facet_wrap(~CellType, scales="free_y", nrow=2) +
    labs(title=title, x=NULL, y="Proportion") +
    theme_classic(base_size=12) +
    theme(legend.position="none",
          strip.background=element_blank(),
          strip.text=element_text(size=10, face="bold"),
          strip.clip="off",
          axis.text.x=element_text(size=9, angle=30, hjust=1),
          axis.text.y=element_text(size=8),
          plot.title=element_text(size=12, face="bold", hjust=0.5),
          panel.spacing=unit(1.5, "lines"),
          plot.margin=margin(10,10,10,10))
}

make_mcp_plot <- function(df, title) {
  cell_cols <- names(mcp_labels)
  stats <- do.call(rbind, lapply(cell_cols, function(ct) {
    b <- df[df$Timepoint=="Baseline",]; w <- df[df$Timepoint=="Week4",]
    b <- b[order(b$PairID),];           w <- w[order(w$PairID),]
    data.frame(CellType=ct, p=wilcox.test(b[[ct]], w[[ct]], paired=TRUE)$p.value)
  }))
  stats$fdr <- p.adjust(stats$p, method="BH")
  
  keep_mcp <- stats$CellType[stats$fdr < 0.05]
  
  long <- do.call(rbind, lapply(keep_mcp, function(ct) {
    data.frame(CellType=mcp_labels[ct], Score=df[[ct]],
               Timepoint=df$Timepoint, PairID=df$PairID)
  }))
  labs_ordered   <- mcp_labels[keep_mcp]
  long$Timepoint <- factor(long$Timepoint, levels=c("Baseline","Week4"))
  long$CellType  <- factor(long$CellType, levels=labs_ordered)
  
  ymax   <- tapply(long$Score, long$CellType, max)
  sig_df <- data.frame(
    CellType = factor(labs_ordered, levels=labs_ordered),
    label    = sig_label(stats$fdr[match(keep_mcp, stats$CellType)]),
    y        = ymax[labs_ordered] * 1.12
  )
  
  ggplot(long, aes(x=Timepoint, y=Score)) +
    geom_line(aes(group=PairID), color="grey70", linewidth=0.2, alpha=0.3) +
    geom_boxplot(aes(fill=Timepoint), outlier.shape=NA, alpha=0.85,
                 width=0.5, color="grey30", linewidth=0.4) +
    geom_text(data=sig_df, aes(x=1.5, y=y, label=label),
              inherit.aes=FALSE, size=5) +
    scale_fill_manual(values=c(Baseline=COL_BASELINE, Week4=COL_WEEK4)) +
    scale_y_continuous(expand=expansion(mult=c(0.05, 0.25))) +
    facet_wrap(~CellType, scales="free_y", nrow=2) +
    labs(title=title, x=NULL, y="MCP-counter Score") +
    theme_classic(base_size=12) +
    theme(legend.position="none",
          strip.background=element_blank(),
          strip.text=element_text(size=10, face="bold"),
          strip.clip="off",
          axis.text.x=element_text(size=9, angle=30, hjust=1),
          axis.text.y=element_text(size=8),
          plot.title=element_text(size=12, face="bold", hjust=0.5),
          panel.spacing=unit(1.5, "lines"),
          plot.margin=margin(10,10,10,10))
}

p_cib_gse <- make_cib_plot(cib_gse_paired, stats_gse, keep_cts,
                           "CIBERSORTx: GSE89403 (n=86 pairs)")
p_mcp_gse <- make_mcp_plot(mcp_gse_df, "MCP-counter: GSE89403 (n=86 pairs)")

print(p_cib_gse)
print(p_mcp_gse)
# Save figure panels
ggsave(file.path(fig_dir, "Fig2e_cibersortx.pdf"), p_cib_gse, width = 10, height = 6)
ggsave(file.path(fig_dir, "Fig2e_cibersortx.png"), p_cib_gse, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Fig2e_mcpcounter.pdf"), p_mcp_gse, width = 8, height = 6)
ggsave(file.path(fig_dir, "Fig2e_mcpcounter.png"), p_mcp_gse, width = 8, height = 6, dpi = 300)

# CIBERSORTx p values
print(stats_gse[, c("CellType", "p", "fdr")])

# MCP-counter p values
cell_cols <- names(mcp_labels)
mcp_stats <- do.call(rbind, lapply(cell_cols, function(ct) {
  b <- mcp_gse_df[mcp_gse_df$Timepoint=="Baseline",]
  w <- mcp_gse_df[mcp_gse_df$Timepoint=="Week4",]
  b <- b[order(b$PairID),]; w <- w[order(w$PairID),]
  data.frame(CellType=ct, p=wilcox.test(b[[ct]], w[[ct]], paired=TRUE)$p.value)
}))
mcp_stats$fdr <- p.adjust(mcp_stats$p, method="BH")
print(mcp_stats)

# Save results for all cell types
cat(sprintf("Pairs analysed: CIBERSORTx %d | MCP-counter %d\n",
            length(unique(cib_gse_paired$PairID)), length(unique(mcp_gse_df$PairID))))
write.csv(data.frame(Method = "CIBERSORTx", stats_gse[, c("CellType","mean_bl","mean_wk4","p","fdr")]),
          file.path(out_dir, "SuppTable_cibersortx_paired_wilcoxon.csv"), row.names = FALSE)
mcp_means <- do.call(rbind, lapply(mcp_stats$CellType, function(ct) data.frame(
  CellType = ct,
  mean_bl  = mean(mcp_gse_df[mcp_gse_df$Timepoint=="Baseline", ct]),
  mean_wk4 = mean(mcp_gse_df[mcp_gse_df$Timepoint=="Week4",    ct]))))
write.csv(data.frame(Method = "MCP-counter", merge(mcp_means, mcp_stats, by = "CellType")),
          file.path(out_dir, "SuppTable_mcpcounter_paired_wilcoxon.csv"), row.names = FALSE)
