library(limma)
library(edgeR)
library(data.table)
library(tidyr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(GEOquery)
library(oligo)

# =============================================================================
# GSE31348
# =============================================================================

id<-'GSE31348'
data<-getGEO(id, GSEMatrix=F)
ids=c(); titles=c()
for (g in data@gsms) {
  ids<-c(ids, g@header$geo_accession)
  titles<-c(titles, g@header$title)
}
sample_metadata <- data.frame(ids,titles)
sample_metadata$condition[grepl("diagnosis", sample_metadata$titles)] <- "DX"
sample_metadata$condition[grepl("week4", sample_metadata$titles)] <- "Week4"
sample_metadata<- sample_metadata[!is.na(sample_metadata$condition),]

celfiles <- list.celfiles("data/mrna/gse31348/GSE31348/CEL", listGzipped=TRUE, full=TRUE)
cel_data <- read.celfiles(celfiles)
norm<-rma(cel_data)
expres_matrix<-exprs(norm)

gpl<-getGEO('GPL570')
probe_map <- gpl@dataTable@table[,c("ID", "ENTREZ_GENE_ID", "Gene Symbol")]
ent <- probe_map$ENTREZ_GENE_ID[match(rownames(expres_matrix), probe_map$ID)]
gene_expr = avereps(expres_matrix, ID = ent)
colnames(gene_expr) <- sub("_.*", "", colnames(gene_expr))
sample_metadata <- sample_metadata[order(sample_metadata$ids),]
gene_expr <- gene_expr[,order(colnames(gene_expr))]
gene_expr<- gene_expr[,colnames(gene_expr) %in% sample_metadata$ids]
stopifnot(all.equal(colnames(gene_expr), sample_metadata$ids))

f <- factor(sample_metadata$condition, levels=c("DX", "Week4"))
f <- relevel(f, ref = "DX")
design <- model.matrix(~ f)
fit <- lmFit(gene_expr, design)
fit2 <- eBayes(fit)
top_table= topTable(fit2, coef = "fWeek4", n=Inf, confint = T)
affy<-as.data.frame(top_table)
affy <- tibble::rownames_to_column(affy, "entrez_id")
write.table(affy,"output/meta/tables/gse31348_toptab.tsv", sep = "\t", row.names = F, quote = F)

mart_affy <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
entrez_map <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                    filters = "entrezgene_id", values = affy$entrez_id, mart = mart_affy)
affy$entrez_id <- as.integer(affy$entrez_id)
affy$hgnc_symbol<- entrez_map$hgnc_symbol[match(affy$entrez_id, entrez_map$entrezgene_id)]

# =============================================================================
# GSE89403
# =============================================================================

counts_dir <- "data/mrna/gse89403/counts"
count_files <- list.files(path=counts_dir, full.names = T)
dat<- lapply(count_files, function(file){
  df<- fread(file)
  df_sub<- df[,c(1,7)]
  colnames(df_sub)<- c("GeneID", gsub(".*_([0-9]+)_.*", "\\1", colnames(df)[7]))
  return(df_sub)
})

stopifnot(all(sapply(dat[-1], function(df) all.equal(df$GeneID, dat[[1]]$GeneID))))
mat<- Reduce(function(x,y) merge(x,y,by="GeneID", all=T), dat)
mat <- do.call(cbind, lapply(dat, function(df) df[, -1]))
mat<- as.data.frame(mat)
mat<-mat[,order(as.numeric(colnames(mat)))]
rownames(mat)<- dat[[1]]$GeneID
colnames(mat) <- gsub(".*/(SRR[0-9]+)\\.bam", "\\1", colnames(mat))

metadata<- read.table("metadata/gse89403/SraRunTable.csv", sep=",", header=T)
metadata<- metadata[,c(1,24,25,28,30)]
stopifnot(all.equal(metadata$Run, colnames(mat)))
has_both <- with(metadata, table(subject, time)) >= 1
valid_subjects <- rownames(has_both[rowSums(has_both) == 2, ])
metadata <- metadata[metadata$subject %in% valid_subjects, ]

metadata$tech_group <- paste(metadata$subject, metadata$time, sep = "_")
mat_filtered <- mat[, metadata$Run]
summed_counts <- rowsum(t(mat_filtered), group = metadata$tech_group)
summed_counts <- t(summed_counts)

sample_info_gse <- data.frame(sample_name = colnames(summed_counts)) |>
  separate(sample_name, into = c("subject", "timepoint"), sep = "_", extra = "merge")
sample_info_gse <- as.data.frame(sample_info_gse)
sample_info_gse$sample_name <- colnames(summed_counts)
sample_info_gse$timepoint <- factor(sample_info_gse$timepoint, levels = c("DX", "week_4"))

keep <- rowSums(cpm(summed_counts) > 1) >= 10
filtered_counts_gse <- summed_counts[keep, ]
dge_gse <- DGEList(counts = filtered_counts_gse)
dge_gse <- calcNormFactors(dge_gse)
design_gse <- model.matrix(~ subject + timepoint, data = sample_info_gse)
v_gse <- voom(dge_gse, design_gse, plot = FALSE)
fit_gse <- lmFit(v_gse, design_gse)
fit_gse <- eBayes(fit_gse)
tab_89403 <- topTable(fit_gse, coef = "timepointweek_4", number = Inf)

# =============================================================================
# STUDY COHORT
# =============================================================================

counts_dir <- "data/mrna/study_cohort/counts_paired"
count_files <- list.files(path=counts_dir, full.names = T)
dat<- lapply(count_files, function(file){
  df<- fread(file)
  df_sub<- df[,c(1,7)]
  colnames(df_sub)<- c("GeneID", gsub(".*_([0-9]+)_.*", "\\1", colnames(df)[7]))
  return(df_sub)
})

stopifnot(all(sapply(dat[-1], function(df) all.equal(df$GeneID, dat[[1]]$GeneID))))
mat<- Reduce(function(x,y) merge(x,y,by="GeneID", all=T), dat)
mat <- do.call(cbind, lapply(dat, function(df) df[, -1]))
mat<- as.data.frame(mat)
mat<-mat[,order(as.numeric(colnames(mat)))]
rownames(mat)<- dat[[1]]$GeneID

sample_names<- colnames(mat)
condition<- ifelse(as.numeric(sample_names)<=22,"DX","week_4")
subject<- c(rep(1:11), rep(1:11))
coldata<- data.frame(row.names=sample_names, subject=factor(subject),condition=factor(condition))
counts<- mat
keep <- rowSums(cpm(counts) > 1) >= 3
counts <- counts[keep, ]
dge_study <- DGEList(counts)
dge_study <- calcNormFactors(dge_study)

coldata <- coldata[match(colnames(counts), rownames(coldata)), ]
sample <- factor(coldata$subject)
timepoint <- factor(coldata$condition, levels = c("DX", "week_4"))
design_study <- model.matrix(~sample+timepoint)
v_study <- voom(dge_study, design_study, plot = FALSE)
fit_study <- lmFit(v_study, design_study)
fit_study <- eBayes(fit_study)
tab_study <- topTable(fit_study, coef = "timepointweek_4", number = Inf)

# =============================================================================
# GENE SYMBOL MAPPING
# =============================================================================

tab_89403$ensembl_id <- sub("\\..*", "", rownames(tab_89403))
tab_study$ensembl_id <- sub("\\..*", "", rownames(tab_study))

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ens_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = unique(c(tab_89403$ensembl_id, tab_study$ensembl_id)),
                 mart = mart)

tab_89403$hgnc_symbol<- ens_map$hgnc_symbol[match(tab_89403$ensembl_id, ens_map$ensembl_gene_id)]
tab_study$hgnc_symbol<- ens_map$hgnc_symbol[match(tab_study$ensembl_id, ens_map$ensembl_gene_id)]
tab_89403<- tab_89403[!is.na(tab_89403$hgnc_symbol) & !duplicated(tab_89403$hgnc_symbol),]
tab_study<- tab_study[!is.na(tab_study$hgnc_symbol) & !duplicated(tab_study$hgnc_symbol),]

# =============================================================================
# STOUFFER'S META-ANALYSIS
# =============================================================================

tab_89403$z_score <- (-1) * qnorm(0.5 * tab_89403$P.Value) * sign(tab_89403$logFC)
tab_study$z_score <- (-1) * qnorm(0.5 * tab_study$P.Value) * sign(tab_study$logFC)

all_studies <- rbind(
  data.frame(tab_89403[c("hgnc_symbol", "logFC", "z_score")], study = "GSE89403"),
  data.frame(tab_study[c("hgnc_symbol", "logFC", "z_score")], study = "study_cohort")
)

weights<- c("GSE89403" = sqrt(86), "study_cohort"= sqrt(11))

meta_res<- by(all_studies, all_studies$hgnc_symbol, function(x){
  w<- weights[x$study]; z<- x$z_score
  combined_z<- sum(w*z)/sqrt(sum(w^2))
  combined_p<- 2*pnorm(abs(combined_z), lower.tail = F)
  meta_logFC<-sum(w*x$logFC)/sum(w)
  dir<- ifelse(all(sign(x$logFC)==1) | all(sign(x$logFC)==-1),"Y","N")
  data.frame(hgnc_symbol=x$hgnc_symbol[1], n_studies=length(z),
             combined_z=combined_z, pval_meta=combined_p,
             logFC_meta=meta_logFC, dir_consistent=dir)
})

res<-do.call(rbind,meta_res)
res<- res[!res$hgnc_symbol=="",]
res_filt <- res[res$n_studies==2,]
res_filt$padj<- p.adjust(res_filt$pval_meta, method="fdr")
res_filt$sig <- ifelse(
  res_filt$dir_consistent == "Y" & res_filt$padj < 0.05 & abs(res_filt$logFC_meta) >= 1,
  ifelse(res_filt$logFC_meta > 0, "Upregulated", "Downregulated"), "Not significant")
degs<- res_filt[res_filt$padj<0.05 & abs(res_filt$logFC_meta)>=1 & res_filt$dir_consistent=="Y",]

write.csv(degs, "output/meta/tables/meta_degs.csv", row.names = F, quote=F)
write.table(tab_89403,"output/meta/tables/gse89403_toptab.tsv", quote = FALSE, sep="\t")
write.table(tab_study,"output/meta/tables/study_cohort_toptab.tsv", quote=F, sep="\t")

# =============================================================================
# FIGURES
# =============================================================================
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(ggplotify)

# ── Project colors ────────────────────────────────────────────────────────────
COL_GSE89403 <- "lightblue"
COL_COHORT   <- "lightsalmon"
COL_GSE31348 <- "#81B29A"
COL_BASELINE <- "bisque1"
COL_WEEK4    <- "tomato3"
COL_UP       <- "#C1121F"
COL_DOWN     <- "#1D4E89"
COL_NS       <- "#ADB5BD"

# ── Shared theme ──────────────────────────────────────────────────────────────
pub_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title    = element_text(size = 13),
    axis.text     = element_text(size = 11, color = "black"),
    plot.title    = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 10),
    axis.line     = element_line(color = "black", linewidth = 0.6),
    axis.ticks    = element_line(color = "black", linewidth = 0.5),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin   = margin(8, 10, 8, 8)
  )

#Panel A
library(ComplexHeatmap)
library(circlize)
library(ggplotify)

# ── Expression matrices ───────────────────────────────────────────────────────
gse_expr   <- v_gse$E
study_expr <- v_study$E

rownames(gse_expr)   <- tab_89403$hgnc_symbol[match(rownames(gse_expr),   rownames(tab_89403))]
rownames(study_expr) <- tab_study$hgnc_symbol[match(rownames(study_expr), rownames(tab_study))]

gse_expr   <- gse_expr[rownames(gse_expr)     %in% degs$hgnc_symbol, ]
study_expr <- study_expr[rownames(study_expr) %in% degs$hgnc_symbol, ]

common_genes <- intersect(rownames(gse_expr), rownames(study_expr))
gse_z        <- t(scale(t(gse_expr[common_genes, ])))
study_z      <- t(scale(t(study_expr[common_genes, ])))
combined_z   <- cbind(gse_z, study_z)

# ── HPA verification ──────────────────────────────────────────────────────────
hpa      <- read.table("data/rna_immune_cell.tsv", sep="\t", header=T)
hpa_filt <- hpa[hpa$Gene.name %in% rownames(combined_z), ]
hpa_top  <- do.call(rbind, lapply(split(hpa_filt, hpa_filt$Gene.name), function(x) {
  x <- x[order(x$nTPM, decreasing=TRUE), ]
  data.frame(gene      = x$Gene.name[1],
             top_cells = paste(x$Immune.cell[1:min(2,nrow(x))],
                               round(x$nTPM[1:min(2,nrow(x))],1),
                               sep="=", collapse=" | "))
}))
print(hpa_top[order(hpa_top$gene), ], row.names=FALSE)

# ── Cell type map (HPA-validated) ─────────────────────────────────────────────
celltype_map <- c(
  SOCS3="Neutrophil", CXCL8="Neutrophil", SERPING1="Neutrophil",
  BATF2="Neutrophil", GBP6="Neutrophil", FCGR1BP="Neutrophil",
  CD177="Neutrophil", ANXA3="Neutrophil", ANKRD22="Neutrophil",
  KREMEN1="Neutrophil", GALNT14="Neutrophil", SMTNL1="Neutrophil",
  CNGB1="Neutrophil", FGF13="Neutrophil", SHOC1="Neutrophil",
  SLC26A8="Neutrophil", FAM20A="Neutrophil",
  
  TLR5="Monocyte", CASP5="Monocyte",
  C1QA="Monocyte", C1QB="Monocyte", C1QC="Monocyte",
  C2="Monocyte", FCGR1A="Monocyte", MCEMP1="Monocyte",
  MYO10="Monocyte", PDCD1LG2="Monocyte", CD274="Monocyte",
  
  AIM2="B cell", HRK="B cell",
  ZDHHC19="B cell", `IGLV7-46`="B cell",
  
  NRN1="T cell", SEPTIN4="T cell",
  
  IL5RA="Eosinophil", SIGLEC8="Eosinophil", PRSS33="Eosinophil",
  SPNS3="Eosinophil", PMP22="Eosinophil", PLAAT5="Eosinophil",
  ALOX15="Eosinophil",
  
  CACNG6="Basophil",
  
  PLVAP="Dendritic cell",
  
  LINC02073="lncRNA", LINC02972="lncRNA", MEG3="lncRNA",
  LINC01093="lncRNA", LINC02751="lncRNA", LINC02967="lncRNA",
  LINC03078="lncRNA", `BASP1-AS1`="lncRNA"
)

celltype_map["CACNA1E"] <- "Other"
celltype_map["ITGA7"]   <- "Other"
celltype_map["MAOB"]    <- "Other"
celltype_map["LIPM"]    <- "Other"
celltype_map["CAPN13"]  <- "Other"
celltype_map["MGAM2"]   <- "Other"
celltype_map["SFRP2"]   <- "Other"

other_genes <- setdiff(rownames(combined_z), names(celltype_map))
celltype_map[other_genes] <- "Other"

# ── Order cell types by number of genes (descending) ─────────────────────────
ct_counts <- sort(table(celltype_map[rownames(combined_z)]), decreasing=TRUE)
ct_order  <- names(ct_counts)

# ── Print percentages ─────────────────────────────────────────────────────────
total <- sum(ct_counts)
ct_pct <- data.frame(
  CellType   = names(ct_counts),
  N          = as.integer(ct_counts),
  Percentage = round(100 * as.integer(ct_counts) / total, 1)
)
print(ct_pct, row.names=FALSE)

# ── Row ordering: hclust within each cell type group ─────────────────────────
ordered_genes <- unlist(lapply(ct_order, function(ct) {
  g <- intersect(names(celltype_map)[celltype_map == ct], rownames(combined_z))
  if (length(g) <= 2) return(g)
  g[hclust(dist(combined_z[g, ]), method="ward.D2")$order]
}))

combined_z <- combined_z[ordered_genes, ]
row_split  <- factor(celltype_map[ordered_genes], levels = ct_order)

# ── Colors ────────────────────────────────────────────────────────────────────
heatmap_cols <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#1D4E89", "#2E86AB", "white", "#E07A5F", "#7D0A0A")
)

ct_levels <- ct_order[ct_order %in% unique(celltype_map[ordered_genes])]
ct_pal <- setNames(
  colorRampPalette(c("#2A9D8F","#264653","#457B9D","#E9C46A",
                     "#F4A261","#8338EC","#E63946","#B5838D",
                     "#ADB5BD"))(length(ct_levels)),
  ct_levels
)

# ── Column annotation ─────────────────────────────────────────────────────────
col_ha <- HeatmapAnnotation(
  Study     = c(rep("GSE89403",     ncol(gse_z)),
                rep("Study Cohort", ncol(study_z))),
  Timepoint = c(ifelse(sample_info_gse$timepoint == "DX", "Baseline", "Week 4"),
                ifelse(as.numeric(colnames(study_z)) <= 22, "Baseline", "Week 4")),
  col = list(
    Study     = c("GSE89403" = COL_GSE89403, "Study Cohort" = COL_COHORT),
    Timepoint = c("Baseline" = COL_BASELINE, "Week 4" = COL_WEEK4)
  ),
  annotation_name_gp      = gpar(fontsize = 9, fontface = "plain"),
  annotation_legend_param = list(
    Study     = list(title_gp  = gpar(fontsize = 9, fontface = "plain"),
                     labels_gp = gpar(fontsize = 8)),
    Timepoint = list(title_gp  = gpar(fontsize = 9, fontface = "plain"),
                     labels_gp = gpar(fontsize = 8))
  )
)

# ── Row annotation ────────────────────────────────────────────────────────────
row_ha_left <- rowAnnotation(
  `Cell Type` = celltype_map[ordered_genes],
  col = list(
    `Cell Type` = ct_pal
  ),
  show_annotation_name    = FALSE,
  annotation_legend_param = list(
    `Cell Type` = list(title_gp  = gpar(fontsize = 9, fontface = "plain"),
                       labels_gp = gpar(fontsize = 8))
  )
)

# ── Draw ──────────────────────────────────────────────────────────────────────
pA_ch <- Heatmap(
  combined_z,
  name                      = "Z-score",
  col                       = heatmap_cols,
  cluster_rows              = FALSE,
  cluster_columns           = TRUE,
  clustering_method_columns = "ward.D2",
  show_column_dend          = FALSE,
  show_column_names         = FALSE,
  show_row_names            = TRUE,
  row_names_side            = "right",
  row_names_gp              = gpar(fontsize = 7.5),
  row_split                 = row_split,
  row_title                 = NULL,
  row_gap                   = unit(2, "mm"),
  top_annotation            = col_ha,
  left_annotation           = row_ha_left,
  border                    = FALSE,
  heatmap_legend_param      = list(
    title_gp      = gpar(fontsize = 9, fontface = "plain"),
    labels_gp     = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

pA <- as.ggplot(grid.grabExpr(draw(pA_ch)))
pA

#Volcano Plot

top_up   <- degs[degs$sig == "Upregulated", ]
top_up   <- top_up[order(top_up$padj), ][1:5, "hgnc_symbol"]
top_down <- degs[degs$sig == "Downregulated", ]
top_down <- top_down[order(top_down$padj), ][1:5, "hgnc_symbol"]
label_genes <- c(top_up, top_down)

vol_dat <- res_filt
vol_dat$sig <- factor(
  ifelse(vol_dat$padj < 0.05 & abs(vol_dat$logFC_meta) >= 1 & vol_dat$dir_consistent == "Y",
         ifelse(vol_dat$logFC_meta > 0, "Upregulated", "Downregulated"), "NS"),
  levels = c("Upregulated", "Downregulated", "NS")
)
vol_dat$label <- ifelse(vol_dat$hgnc_symbol %in% label_genes, vol_dat$hgnc_symbol, NA)
lab_dn <- subset(vol_dat, !is.na(label) & logFC_meta < 0)
lab_up <- subset(vol_dat, !is.na(label) & logFC_meta > 0)

pVol <- ggplot(vol_dat, aes(x = logFC_meta, y = -log10(padj))) +
  geom_point(data = subset(vol_dat, sig == "NS"),
             shape = 16, size = 0.4,
             color = "grey82", alpha = 0.6) +
  geom_point(data = subset(vol_dat, sig == "Downregulated"),
             shape = 21, size = 3.5, stroke = 0.3,
             fill = COL_DOWN, color = "white") +
  geom_point(data = subset(vol_dat, sig == "Upregulated"),
             shape = 21, size = 3.5, stroke = 0.3,
             fill = COL_UP, color = "white") +
  geom_text_repel(
    data    = lab_dn,
    mapping = aes(x = logFC_meta, y = -log10(padj), label = label),
    color = COL_DOWN, size = 3.8, fontface = "italic",
    max.overlaps = Inf, direction = "y", hjust = 1,
    xlim = c(-Inf, -2.2),
    box.padding = 0.5, point.padding = 0.5,
    segment.size = 0.25, segment.color = "grey60",
    min.segment.length = 0,
    force = 0.5, force_pull = 0,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data    = lab_up,
    mapping = aes(x = logFC_meta, y = -log10(padj), label = label),
    color = COL_UP, size = 3.8, fontface = "italic",
    max.overlaps = Inf, direction = "y", hjust = 0,
    xlim = c(1.8, Inf),
    box.padding = 0.5, point.padding = 0.5,
    segment.size = 0.25, segment.color = "grey60",
    min.segment.length = 0,
    force = 0.5, force_pull = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("Upregulated" = COL_UP, "Downregulated" = COL_DOWN, "NS" = "grey82"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Upregulated" = COL_UP, "Downregulated" = COL_DOWN, "NS" = "grey82"),
    name = NULL
  ) +
  #scale_x_continuous(expand = expansion(mult = c(0.45, 0.6))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(clip = "off", xlim=c(-4,4)) +
  labs(
    title = "Meta-analysis Differential Expression",
    x     = expression(log[2] ~ FC),
    y     = expression(-log[10] ~ adj.p-value)
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_text(size = 14, color = "black"),
    axis.title       = element_text(size = 15),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position  = "top",
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(0.5, "cm"),
    plot.margin      = margin(10, 20, 10, 20),
    aspect.ratio     = 0.9
  ) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 4), nrow = 1),
    color = "none"
  )

pVol

## 

# ═════════════════════════════════════════════════════════════════════════════
# PANEL B — Cross-cohort scatter
# ═════════════════════════════════════════════════════════════════════════════

# ── scatter_dat for pB ────────────────────────────────────────────────────────
lfc_89403 <- tab_89403[tab_89403$hgnc_symbol %in% degs$hgnc_symbol,
                       c("hgnc_symbol", "logFC")]
colnames(lfc_89403)[2] <- "logFC_89403"

lfc_study <- tab_study[tab_study$hgnc_symbol %in% degs$hgnc_symbol,
                       c("hgnc_symbol", "logFC")]
colnames(lfc_study)[2] <- "logFC_study"

scatter_dat <- merge(lfc_89403, lfc_study, by = "hgnc_symbol")
scatter_dat <- merge(scatter_dat, degs[, c("hgnc_symbol", "sig")], by = "hgnc_symbol")
r_val       <- round(cor(scatter_dat$logFC_89403, scatter_dat$logFC_study), 3)

# ── valid_dat for pC ──────────────────────────────────────────────────────────
affy_clean <- affy[!is.na(affy$hgnc_symbol) & affy$hgnc_symbol != "", ]
affy_clean <- affy_clean[!duplicated(affy_clean$hgnc_symbol), c("hgnc_symbol", "logFC")]
colnames(affy_clean)[2] <- "logFC_affy"

valid_dat <- merge(degs[, c("hgnc_symbol", "logFC_meta", "sig")],
                   affy_clean, by = "hgnc_symbol")
valid_dat$concordant <- sign(valid_dat$logFC_meta) == sign(valid_dat$logFC_affy)

n_tested       <- nrow(valid_dat)
n_concordant   <- sum(valid_dat$concordant)
pct_concordant <- round(100 * n_concordant / n_tested, 1)
fisher_p       <- fisher.test(
  matrix(c(n_concordant, n_tested - n_concordant,
           n_tested - n_concordant, n_concordant), nrow = 2)
)$p.value
valid_lab <- paste0(n_concordant, "/", n_tested, " concordant (",
                    pct_concordant, "%)\np = ", signif(fisher_p, 2))

pB <- ggplot(scatter_dat,
             aes(x = logFC_89403, y = logFC_study, fill = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              fill = "grey88", linewidth = 0.7) +
  geom_point(shape = 21, size = 3.5, stroke = 0.3,
             color = "white", alpha = 0.9) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", r_val),
           hjust = 1.15, vjust = -1,
           size = 5, color = "black") +
  scale_fill_manual(
    values = c("Upregulated" = COL_UP, "Downregulated" = COL_DOWN),
    name   = NULL
  ) +
  labs(
    title = "Cross-Cohort Concordance",
    x     = expression(GSE89403 ~ log[2] ~ FC),
    y     = expression("Study Cohort" ~ log[2] ~ FC)
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_text(size = 14, color = "black"),
    axis.title       = element_text(size = 15),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position  = "top",
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(0.5, "cm"),
    plot.margin      = margin(10, 20, 10, 20),
    aspect.ratio     = 1
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4), nrow = 1))

pB


## Panel C
pC <- ggplot(valid_dat,
             aes(x = logFC_meta, y = logFC_affy)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey40", linewidth = 0.6) +
  geom_point(data = subset(valid_dat, concordant == TRUE),
             aes(fill = sig), shape = 21, size = 3.5,
             stroke = 0.3, color = "white", alpha = 0.9) +
  geom_point(data = subset(valid_dat, concordant == FALSE),
             shape = 4, size = 3.5, stroke = 1.2,
             color = "#F4A261") +
  annotate("point", x = 0.3, y = -1.85, shape = 21,
           size = 3.5, fill = "grey50", color = "white") +
  annotate("text", x = 0.5, y = -1.85,
           label = "Concordant", hjust = 0, size = 3.8) +
  annotate("point", x = 0.3, y = -2.05, shape = 4,
           size = 3.5, color = "#F4A261", stroke = 1.2) +
  annotate("text", x = 0.5, y = -2.05,
           label = "Discordant", hjust = 0, size = 3.8) +
  scale_fill_manual(
    values = c("Upregulated" = COL_UP, "Downregulated" = COL_DOWN)
  ) +
  labs(
    title = "Independent Validation (GSE31348)",
    x     = expression("Meta-analysis" ~ log[2] ~ FC),
    y     = expression(GSE31348 ~ log[2] ~ FC)
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_text(size = 14, color = "black"),
    axis.title       = element_text(size = 15),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position  = "none",
    plot.margin      = margin(10, 20, 10, 20),
    aspect.ratio     = 1
  )

pC