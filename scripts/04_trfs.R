library(DESeq2); library(ggplot2); library(dplyr)
library(ComplexHeatmap); library(circlize)
library(patchwork); library(Biostrings)
library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38); library(tidyr)

theme_set(theme_classic(base_size = 13, base_family = "Helvetica"))


# ── Colors ────────────────────────────────────────────────────────────────────
COL_BASELINE <- "#D4A96A"; COL_WEEK4 <- "tomato3"
COL_UP <- "#C1121F"; COL_DOWN <- "#1D4E89"; COL_NS <- "#ADB5BD"
clean_name <- function(x) gsub("_", " ", x)

# ── Load & DESeq2 ─────────────────────────────────────────────────────────────
trf_data <- read.table("data/trf_counts.csv", header=T, sep=",", row.names=1)
baseline_cols  <- paste0("X", c(1,2,3,11,13,14,15,17,20,21,22))
treatment_cols <- paste0("X", c(23,24,25,33,35,36,37,39,42,43,44))
trf_counts <- trf_data[, c(baseline_cols, treatment_cols)]

condition <- c(rep("baseline", 11), rep("treatment", 11))
subject   <- rep(c(1,2,3,11,13,14,15,17,20,21,22), 2)
col_data  <- data.frame(subject = factor(subject), condition = factor(condition),
                        row.names = colnames(trf_counts))

dds  <- DESeqDataSetFromMatrix(trf_counts, col_data, ~ subject + condition)
keep <- rowSums(counts(dds) >= 10) >= 11
dds  <- dds[keep, ]
dds  <- DESeq(dds)
res  <- results(dds, name="condition_treatment_vs_baseline")
res_df   <- as.data.frame(res[complete.cases(res) & res$baseMean > 50, ])
DEGs_n11 <- res_df[abs(res_df$log2FoldChange) > 0.58 & res_df$padj < 0.05, ]

# ── Entity mapping ────────────────────────────────────────────────────────────
entity_map <- read.csv("trf-clusters.csv", header=TRUE)
entity_map_clean <- entity_map
entity_map_clean[,3] <- apply(entity_map[,3:ncol(entity_map)], 1, function(x) {
  valid <- x[!is.na(x) & x != ""]
  if(length(valid) > 0) paste(valid, collapse=",") else ""
})
entity_map_clean <- entity_map_clean[,1:3]
colnames(entity_map_clean) <- c("entity_name", "trf_type", "clusters")
entity_map_clean$is_mito <- grepl("trnaMT_", entity_map_clean$clusters)

results <- merge(data.frame(tRF=rownames(DEGs_n11), DEGs_n11),
                 entity_map_clean[,c("entity_name","trf_type","is_mito")],
                 by.x="tRF", by.y="entity_name", all.x=T)

trfs <- read.csv("trf_clusters_clean.csv", sep=",", header=T)
results_df <- as.data.frame(results)

# ── Shared objects ────────────────────────────────────────────────────────────
aa_labels <- gsub("_", " ", sub("_Comb.*", "", results_df$tRF))
aa_unique       <- unique(aa_labels)
aa_palette      <- setNames(
  colorRampPalette(c("#E07B54","#6A9FB5","#8FBC8F","#BDA0CB",
                     "#F0C987","#7EC8C8","#D98880","#A9CCE3"))(length(aa_unique)),
  aa_unique
)
clean_row_names <- clean_name(results_df$tRF)

# ── Cluster mapping ───────────────────────────────────────────────────────────
de_trf_clusters <- trfs |>
  filter(entity_name %in% results_df$tRF) |>
  dplyr::select(entity_name, clusters)

de_trf_expanded <- de_trf_clusters |>
  rowwise() |>
  mutate(tRNA_name = list(strsplit(clusters, ",")[[1]])) |>
  unnest(tRNA_name) |>
  mutate(tRNA_name = trimws(tRNA_name),
         tRNA_name = sub("_Cluster_.*", "", tRNA_name))

# ── Merge reports ─────────────────────────────────────────────────────────────
report_files <- list.files("data/trfs", pattern="\\.tRFs\\.report\\.tsv$", full.names=TRUE)
reports <- lapply(report_files, function(f) {
  snum <- sub(".*DQui_(\\d+)_S.*", "\\1", basename(f))
  read.delim(f, header=TRUE, stringsAsFactors=FALSE) |> mutate(sample = snum)
})
report_all <- bind_rows(reports)
colnames(report_all) <- c("tRNA_name","tRNA_seq","tRF_seq","mismatch",
                          "tRF_type","coordinate","read_count","rp100k","sample")
report_all <- report_all |>
  mutate(Timepoint = ifelse(sample %in% c("1","2","3","11","13","14","15","17","20","21","22"),
                            "Baseline", "Treatment"))

# ── tRF sequences ─────────────────────────────────────────────────────────────
trf_seqs <- report_all |>
  inner_join(de_trf_expanded, by = "tRNA_name") |>
  group_by(entity_name) |>
  slice_max(read_count, n = 1, with_ties = FALSE) |>
  dplyr::select(entity_name, tRF_seq, tRF_type, tRNA_name)

# ── Length summary ────────────────────────────────────────────────────────────
length_summary <- report_all |>
  inner_join(de_trf_expanded, by = "tRNA_name") |>
  mutate(length = nchar(tRF_seq)) |>
  group_by(entity_name, length) |>
  summarise(total_reads = sum(read_count), .groups = "drop") |>
  left_join(results_df |> dplyr::select(tRF, log2FoldChange),
            by = c("entity_name" = "tRF")) |>
  mutate(AminoAcid    = sub("_Comb.*", "", entity_name),
         entity_label = clean_name(entity_name))

# ── Biogenesis count ──────────────────────────────────────────────────────────
bio_count <- data.frame(trf_type = results_df$trf_type) |>
  dplyr::count(trf_type) |>
  mutate(pct   = n / sum(n),
         label = paste0(trf_type, "\nn=", n))

# ── 7-mer seed matching ───────────────────────────────────────────────────────
meta_degs <- read.csv("output/meta/tables/meta_degs.csv")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
utrs <- threeUTRsByTranscript(txdb, use.names = TRUE)

tx2eg  <- AnnotationDbi::select(txdb, keys = names(utrs), columns = "GENEID", keytype = "TXNAME")
tx2sym <- AnnotationDbi::select(org.Hs.eg.db,
                                keys    = unique(tx2eg$GENEID[!is.na(tx2eg$GENEID)]),
                                columns = "SYMBOL", keytype = "ENTREZID")
tx_map <- tx2eg |>
  left_join(tx2sym, by = c("GENEID" = "ENTREZID")) |>
  dplyr::select(transcript = TXNAME, gene = SYMBOL) |>
  filter(!is.na(gene), gene %in% meta_degs$hgnc_symbol)

utrs_filt <- utrs[names(utrs) %in% tx_map$transcript]
utr_seqs  <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, utrs_filt)

seed_matches <- lapply(seq_len(nrow(trf_seqs)), function(i) {
  trf_name <- trf_seqs$entity_name[i]
  seq_str  <- trf_seqs$tRF_seq[i]
  if (nchar(seq_str) < 8) return(NULL)
  seed_rc <- as.character(reverseComplement(DNAString(substr(seq_str, 2, 8))))
  hits <- lapply(names(utr_seqs), function(tx) {
    m <- matchPattern(DNAString(seed_rc), utr_seqs[[tx]], max.mismatch = 0)
    if (length(m) > 0) data.frame(tRF = trf_name, transcript = tx) else NULL
  })
  bind_rows(hits)
})

seed_match_df <- bind_rows(seed_matches) |>
  inner_join(tx_map, by = "transcript") |>
  distinct(tRF, gene) |>
  left_join(meta_degs |> dplyr::select(hgnc_symbol, sig, logFC_meta),
            by = c("gene" = "hgnc_symbol")) |>
  left_join(results_df |> dplyr::select(tRF, log2FoldChange), by = "tRF")

# ── Anticodon fold change ─────────────────────────────────────────────────────
codon_reads <- report_all |>
  inner_join(de_trf_expanded, by = "tRNA_name") |>
  mutate(anticodon = sub(".*tRNA-[A-Za-z]+-([A-Z]+)-.*", "\\1", tRNA_name)) |>
  group_by(entity_name, anticodon, Timepoint) |>
  summarise(total_reads = sum(read_count), .groups = "drop")

codon_totals <- codon_reads |>
  group_by(entity_name, anticodon) |>
  summarise(ac_total = sum(total_reads), .groups = "drop") |>
  group_by(entity_name) |>
  mutate(ac_pct = ac_total / sum(ac_total)) |>
  ungroup() |>
  filter(ac_pct > 0.05)

codon_fc <- codon_reads |>
  inner_join(codon_totals |> dplyr::select(entity_name, anticodon),
             by = c("entity_name", "anticodon")) |>
  tidyr::pivot_wider(names_from = Timepoint,
                     values_from = total_reads, values_fill = 0) |>
  mutate(log2FC       = log2((Treatment + 1) / (Baseline + 1)),
         entity_label = clean_name(entity_name))

all_combos <- expand.grid(
  entity_label = unique(codon_fc$entity_label),
  anticodon    = unique(codon_fc$anticodon),
  stringsAsFactors = FALSE
)
codon_fc_complete <- all_combos |>
  left_join(codon_fc |> dplyr::select(entity_label, anticodon, log2FC),
            by = c("entity_label", "anticodon"))

# ════════════════════════════════════════════════════════════════════════════
# PLOTS
# ════════════════════════════════════════════════════════════════════════════

# ── Panel A: Heatmap ──────────────────────────────────────────────────────────
norm_mat <- log2(counts(dds, normalized = TRUE) + 1)
mat_plot <- norm_mat[results_df$tRF, ]
rownames(mat_plot) <- clean_row_names
mat_scaled <- t(scale(t(mat_plot)))
col_fun    <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "white", COL_UP))
ht_opt(TITLE_PADDING = unit(3, "mm"))
aa_labels <- gsub("_", " ", sub("_Comb.*", "", results_df$tRF))



ht <- Heatmap(mat_scaled,
              col = col_fun,
              top_annotation = HeatmapAnnotation(
                Timepoint = col_data$condition,
                col = list(Timepoint = c(baseline = COL_BASELINE, treatment = COL_WEEK4)),
                annotation_name_gp = gpar(fontsize = 13),
                annotation_legend_param = list(title_gp = gpar(fontsize=13, fontfamily = "Helvetica"),
                                               labels_gp = gpar(fontsize=12))
              ),
              right_annotation = rowAnnotation(
                "Amino Acid" = aa_labels,
                Direction = ifelse(results_df$log2FoldChange > 0, "Up", "Down"),
                col = list("Amino Acid"= aa_palette,
                           Direction = c(Up = COL_UP, Down = COL_DOWN)),
                annotation_name_gp = gpar(fontsize=13),
                annotation_legend_param = list(title_gp = gpar(fontsize=13),
                                               labels_gp = gpar(fontsize=12)),
                width = unit(1.2, "cm")
              ),
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize=13),
              show_column_names = FALSE,
              column_split = factor(col_data$condition, levels=c("baseline","treatment")),
              column_title = NULL,
              heatmap_legend_param = list(title="Z-score",
                                          title_gp=gpar(fontsize=13),
                                          labels_gp=gpar(fontsize=12)),
              border = FALSE
)

pA <- grid.grabExpr(draw(ht, merge_legend = TRUE))
grid.newpage(); grid.draw(pA)

# ── Panel B: Length lollipop ──────────────────────────────────────────────────
pB <- ggplot(length_summary, aes(x = length, y = total_reads, color = AminoAcid)) +
  geom_segment(aes(xend = length, yend = 0), linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ entity_label, ncol = 2, scales = "free") +
  scale_color_manual(values = aa_palette) +
  scale_y_continuous(
    labels = function(x) ifelse(x >= 1000, paste0(round(x/1000,0),"K"), as.character(x)),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "tRF length (nt)", y = "Read count") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size=12),
        axis.text        = element_text(size=11),
        axis.title       = element_text(size=13),
        panel.spacing    = unit(1, "lines"),
        legend.position  = "none")
print(pB)

# ── Panel C: Biogenesis donut ─────────────────────────────────────────────────
pC <- ggplot(bio_count, aes(x = 2, y = n, fill = trf_type)) +
  geom_col(color = "white", linewidth = 1.5) +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = c("tRF-3"="#AED6E8","3' Half"="#F5D5A0","i-tRF"= "#D4B8E0")) +
  theme_void(base_size = 13) +
  theme(legend.position = "none")
print(pC)

# ── Panel D: Seed match dotplot ───────────────────────────────────────────────
seed_plot_df <- seed_match_df |>
  mutate(tRF_label = clean_name(tRF),
         gene_label = gene,
         AminoAcid  = sub(" Comb.*", "", clean_name(tRF)))

pD <- ggplot(seed_plot_df, aes(x = tRF_label, y = gene_label,
                               color = AminoAcid, size = abs(logFC_meta))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = aa_palette) +
  scale_size_continuous(range = c(2,8), breaks = c(1.0,1.2,1.4,1.6), limits = c(1.0,1.7)) +
  labs(y = "Meta-DEGs", color = "Amino acid",
       size  = expression("|log"[2]*"FC|"),
       title = "7-mer seed matches to Meta-DEG 3\u2019UTRs") +
  theme_classic(base_size = 13) +
  theme(axis.text.x      = element_text(angle=40, hjust=1, size=12),
        axis.text.y      = element_text(size=12, face="italic"),
        axis.title       = element_text(size=13),
        panel.grid.major = element_line(color="grey92", linewidth=0.3),
        plot.title       = element_text(size=13, hjust=0.5),
        legend.title     = element_text(size=12),
        legend.text      = element_text(size=11))
print(pD)

# ── Supplementary: Anticodon heatmap ─────────────────────────────────────────
pSupp <- ggplot(codon_fc_complete, aes(x = anticodon, y = entity_label, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(data = codon_fc_complete |> filter(!is.na(log2FC)),
            aes(label = round(log2FC, 1)),
            size = 4.5, color = "grey20") +
  scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP,
                       midpoint = 0, name = expression("log"[2]*"FC"),
                       na.value = "grey93") +
  labs(x = "Anticodon", y = "DE-tRF",
       title = "Anticodon-level abundance change (Treatment vs Baseline)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=11),
        axis.text.y  = element_text(size=11),
        axis.title   = element_text(size=13),
        plot.title   = element_text(hjust=0.5, size=12),
        legend.title = element_text(),
        panel.grid   = element_blank())
print(pSupp)

###
ggplot(anticodon_pct_tp, aes(x = anticodon, y = entity_label,
                             size = pct, color = Timepoint)) +
  geom_point(alpha = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(Baseline = COL_BASELINE, Treatment = COL_WEEK4)) +
  scale_size_continuous(range = c(2, 12), name = "% reads",
                        breaks = c(25, 50, 75, 100)) +
  labs(x = "Anticodon", y = "DE-tRF", color = "Timepoint",
       title = "Anticodon contribution per DE-tRF (Baseline vs Treatment)") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y  = element_text( size = 11),
    axis.title   = element_text( size = 13),
    plot.title   = element_text( hjust = 0.5, size = 12),
    legend.title = element_text(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

pVolc <- ggplot(res_df_plot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 1, size = 1.8, color = "grey60", stroke = 0.4) +
  geom_point(data = subset(res_df_plot, sig == "Up"),
             shape = 16, size = 2.5, color = COL_UP) +
  geom_point(data = subset(res_df_plot, sig == "Down"),
             shape = 16, size = 2.5, color = COL_DOWN) +
  geom_text_repel(aes(label = label), size = 3.2, fontface = "italic",
                  na.rm = TRUE, max.overlaps = 20,
                  segment.color = "grey60", segment.size = 0.3) +
  labs(x = expression("log"[2]*"FC (Week 4 vs Baseline)"),
       y = expression("-log"[10]*"(adj. p-value)")) +
  coord_fixed() +
  theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(axis.title      = element_text(face = "bold"),
        axis.text       = element_text(face = "bold"),
        panel.grid.major = element_line(color = "white", linewidth = 0.5),
        panel.background = element_rect(fill = "white"))

print(pVolc)

###

# ── Write tables ──────────────────────────────────────────────────────────────
trf_de_table <- merge(data.frame(tRF=rownames(DEGs_n11), DEGs_n11),
                      entity_map_clean[, c("entity_name","trf_type","is_mito")],
                      by.x="tRF", by.y="entity_name", all.x=TRUE)
trf_de_table <- trf_de_table[order(trf_de_table$padj), ]
write.csv(trf_de_table, "output/trf/tables/trf_de_results.csv", row.names=FALSE)

write.csv(seed_match_df, "output/trf/tables/trf_seed_matches.csv", row.names=FALSE)

write.csv(trf_seqs, "output/trf/tables/trf_sequences.csv", row.names=FALSE)

anticodon_bar <- codon_reads |>
  group_by(entity_name, anticodon) |>
  mutate(total_across_tp = sum(total_reads)) |>
  ungroup() |>
  filter(total_across_tp >= 100) |>
  left_join(dominant_anticodon, by = "entity_name") |>
  mutate(
    entity_label = clean_name(entity_name),
    is_dominant  = anticodon == dominant_anticodon,
    Timepoint    = factor(Timepoint, levels = c("Baseline", "Treatment"),
                          labels = c("Baseline", "Week 4"))
  )

pAnticodon <- ggplot(anticodon_bar,
                     aes(x = anticodon, y = total_reads, fill = Timepoint,
                         alpha = is_dominant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ entity_label, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Baseline" = COL_BASELINE, "Week 4" = COL_WEEK4)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
  scale_y_continuous(
    labels = function(x) ifelse(x >= 1000, paste0(round(x/1000, 0), "K"), x),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = "Anticodon", y = "Read count", fill = "Timepoint") +
  theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 12),
    axis.text.x      = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 13),
    panel.spacing    = unit(1.2, "lines"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 12)
  )

print(pAnticodon)