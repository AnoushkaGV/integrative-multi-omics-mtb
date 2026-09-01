# =============================================================================
# EARLY5 vs Sweeney3 vs RISK6: signature scoring and ROC comparison
#
# Computes three blood transcriptomic signatures at baseline (DX) and
# compares their ability to predict fast vs. slow sputum culture conversion:
#   - EARLY5   : 5-gene-pair voting score
#   - Sweeney3 : Sweeney_OD_3 signature (via TBSignatureProfiler)
#   - RISK6    : PennNich_RISK_6 signature (via TBSignatureProfiler)
#
# Inputs expected in the environment before running this script:
#   - vsd_mat  : DESeq2 variance-stabilized expression matrix (genes x
#                samples), produced earlier in 05_predictions.R (or wherever
#                the study cohort's DESeq2 object is built)
#   - metadata : sample metadata data frame with columns treatmentresult,
#                timetonegativity, time, tgrv, aligned to colnames(vsd_mat)
#                by rowname/sample ID (IDs suffixed _DX / _week_4)
# =============================================================================

library(ggplot2)

# gene pairs EARLY5/HOPKINS
#DIS3L, DOCK4
# COL9A2, FAM241A 
# PIEZO1, H2BC12 
# PID1, LYPLA1 
# PDE4DIP, GRB10
gp_1 = c("DIS3L", "DOCK4")
gp_2 = c("COL9A2", "FAM241A" )
gp_3 = c("PIEZO1", "H2BC12" )
gp_4 = c("PID1", "LYPLA1" )
gp_5 = c("PDE4DIP", "GRB10")

# create a matrix of zero and one where for each gene pair, if the first gene expression is higher than the second
# is one and if the second is higher than the firs is zero from the gene matrix vsd_mat
gp_matrix = matrix(0, nrow = ncol(vsd_mat), ncol = 5)
colnames(gp_matrix) = c("gp_1", "gp_2", 
                      "gp_3", "gp_4", "gp_5")
rownames(gp_matrix) = colnames(vsd_mat)
for (i in 1:ncol(vsd_mat)) {
  gp_matrix[i, "gp_1"] = ifelse(vsd_mat[gp_1[1], i] > vsd_mat[gp_1[2], i], 1, 0)
  gp_matrix[i, "gp_2"] = ifelse(vsd_mat[gp_2[1], i] > vsd_mat[gp_2[2], i], 1, 0)
  gp_matrix[i, "gp_3"] = ifelse(vsd_mat[gp_3[1], i] > vsd_mat[gp_3[2], i], 1, 0)
  gp_matrix[i, "gp_4"] = ifelse(vsd_mat[gp_4[1], i] > vsd_mat[gp_4[2], i], 1, 0)
  gp_matrix[i, "gp_5"] = ifelse(vsd_mat[gp_5[1], i] > vsd_mat[gp_5[2], i], 1, 0)
}


# create a function to use the previos gp_1, gp_2, gp_3,gp_4 and gp_5 to performt the transformation of 0-1 matrix
# and create a score with the sum of the five gene pairs from a matrix of zero and one where for each gene pair, if the first gene expression is higher than the second gp_1 = c("DIS3L", "DOCK4") gp_2 = c("COL9A2", "FAM241A" ) gp_3 = c("PIEZO1", "H2BC12" ) gp_4 = c("PID1", "LYPLA1" ) gp_5 = c("PDE4DIP", "GRB10"), get the gene value and do the matrix and them sum thee sum

early5 = function(df){

  gene_pairs = list(
    gp_1 = c("DIS3L", "DOCK4"),
    gp_2 = c("COL9A2", "FAM241A"),
    gp_3 = c("PIEZO1", "H2BC12"),
    gp_4 = c("PID1", "LYPLA1"),
    gp_5 = c("PDE4DIP", "GRB10")
  )

  gp_matrix = sapply(gene_pairs, function(pair){
    as.integer(df[pair[1], ] > df[pair[2], ])
  })

  rownames(gp_matrix) = colnames(df)

  score = rowSums(gp_matrix)

  return(score)
}


gp_matrix = as.data.frame(gp_matrix)
gp_matrix$Group = metadata$treatmentresult
gp_matrix$timepoint = metadata$timetonegativity
gp_matrix$conversion = ifelse(gp_matrix$timepoint == 'Week04', 'Fast', "Slow")
gp_matrix$time = metadata$time
gp_matrix$tgrv_iqr = metadata$tgrv / IQR(metadata$tgrv, na.rm=TRUE)

gp_matrix = gp_matrix[complete.cases(gp_matrix$conversion), ]
# keep only ids that have both DX and week_4 the Ids are in the rownames, followed by _DX and _week_4
# se have to treat the IDs and keep only those that have both DX and week_4
ids = unique(gsub("_(DX|week_4)$", "", rownames(gp_matrix)))
ids = ids[ids %in% gsub("_(DX|week_4)$", "", rownames(gp_matrix))]
gp_matrix = gp_matrix[rownames(gp_matrix) %in% paste0(ids, "_DX") | rownames(gp_matrix) %in% paste0(ids, "_week_4"), ]

#retrieve sweeney signature from TBSignatureProfiler package and create a matrix with the expression of the genes in the signature for each sample
sweeey = as.data.frame(t(vsd_mat[TBSignatureProfiler::TBsignatures$Sweeney_OD_3,]))

sweeey$time =  metadata$time

sweeey$conversion = ifelse(metadata$timetonegativity == 'Week04', 'Fast', "Slow")

sweeey = sweeey[complete.cases(sweeey$conversion), ]

sweeney3 = function(df){

  rowMeans(
    df[, c("GBP5","DUSP3")]
  ) - df$KLF2

}

sweeey$Sweeney3 = sweeney3(sweeey)
# retrieve risk6 signarure from TBSignatureProfiler package and create a matrix with the expression of the genes in the signature for each sample
risk6 = as.data.frame(t(vsd_mat[rownames(vsd_mat) %in% TBSignatureProfiler::TBsignatures$PennNich_RISK_6,]))
risk6$time =  metadata$time
risk6$conversion = ifelse(metadata$timetonegativity == 'Week04', 'Fast', "Slow")
risk6 = risk6[complete.cases(risk6$conversion), ]

risk6_fun = function(df){

  down = c("TUBGCP6","TRMT2A","SDR39U1")
  up   = c("GBP2","FCGR1B","SERPING1")

  gm_down = apply(df[, colnames(df) %in% down], 1, mean)
  gm_up   = apply(df[, colnames(df) %in% up], 1, mean)

  risk6_score =  gm_up - gm_down

  return(risk6_score)
}

risk6$RISK6 = risk6_fun(risk6)

# create a roc curve from gp_matrix using the 5 gene pairs from DX timepoint to predict slow and fast conversion
# subjects. the classification is on the conversion colum
library(pROC)
  {
    df = data.frame(outcome = gp_matrix[gp_matrix$time == "DX",]$conversion,
                early5 = rowSums(gp_matrix[gp_matrix$time == "DX",1:5]))
df$outcome = as.factor(df$outcome)
binary_model = glm(outcome ~. , 
                    data = df, 
                    family = "binomial")
name = levels(df$outcome)[levels(df$outcome)  != 0]
gg_df_r= data.frame(pred = fitted(binary_model, 
                                   type = "response"),
                     classified = ifelse(fitted(binary_model,
                                                type = "response") >= .5, 1,0),
                     real = df$outcome) #U can change it

roc_plot = ggplot2::ggplot(gg_df_r, aes(d = real, m = pred)) + 
    plotROC::geom_roc(n.cuts = 10,
                      labels = FALSE,
                      size =2)+
    plotROC::style_roc()+
    ggplot2::geom_abline(slope = 1, 
                         intercept = 0,
                         size = 1,
                         linetype = "dashed",
                         alpha = .4)+
    ggplot2::theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   aspect.ratio = 1)+
    ggplot2::ggtitle("")

par(pty = "s")
roc_EARLY5 = roc(roc_plot$data$real, roc_plot$data$pred)
pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
               col = "#C1121F",
               legacy.axes = T,
               lwd = 3,
               ci=TRUE,
               print.auc=F,
               of="se", 
               specificities=seq(0, 1, .1), 
               ci.type="shape", ci.col=alpha("#C1121F", 0.15),
               asp = NA)

pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
               col = "#C1121F",
               legacy.axes = T,
               print.auc.y = .05,
               lwd = 2,
               ci=TRUE,
               print.auc=TRUE,
               print.auc.pattern = "EARLY5 = %.2f(%.2f-%.2f)",  # = Rounds AUC to 2 digits
               add = T,
               asp = NA)
        
        {
          #sweeney
          df = data.frame(outcome = sweeey[sweeey$time == "DX",]$conversion,
                          Sweeney3 = sweeey[sweeey$time == "DX",]$Sweeney3)
          df$outcome = as.factor(df$outcome)
          binary_model = glm(outcome ~. , 
                        data = df, 
                        family = "binomial")
          name = levels(df$outcome)[levels(df$outcome)  != 0]
          gg_df_r= data.frame(pred = fitted(binary_model, 
                                             type = "response"),
                               classified = ifelse(fitted(binary_model,
                                                          type = "response") >= .5, 1,0),
                               real = df$outcome) #U can change it
          # ROC
          roc_plot = ggplot2::ggplot(gg_df_r, aes(d = real, m = pred)) + 
            plotROC::geom_roc(n.cuts = 10,
                              labels = FALSE,
                              size =2)+
            plotROC::style_roc()+
            ggplot2::geom_abline(slope = 1, 
                                 intercept = 0,
                                 size = 1,
                                 linetype = "dashed",
                                 alpha = .4)+
            ggplot2::theme(panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           aspect.ratio = 1)+
            ggplot2::ggtitle("")
          
          roc_Sweeney3 = roc(roc_plot$data$real, roc_plot$data$pred)
          
          pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
                         col = "#1D4E89",
                         legacy.axes = T,
                         lwd = 3,
                         ci=TRUE,
                         print.auc=F,
                         add = T,
                         of="se", 
                         specificities=seq(0, 1, .1), 
                         ci.type="shape", ci.col=alpha("#1D4E89", 0.15),
                         asp = NA)
          
          pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
                         col = "#1D4E89",
                         print.auc.y = .1,
                         lwd = 2,
                         ci=TRUE,
                         print.auc=TRUE,
                         print.auc.pattern = "Sweeney3 = %.2f(%.2f-%.2f)",  # = Rounds AUC to 2 digits
                         add = T,
                         asp = NA)
        }
        
        {
          #risk6
          df = data.frame(outcome = risk6[risk6$time == "DX",]$conversion,
                          risk6 = risk6[risk6$time == "DX",]$RISK6)
          df$outcome = as.factor(df$outcome)
          binary_model = glm(outcome ~. , 
                        data = df, 
                        family = "binomial")
          name = levels(df$outcome)[levels(df$outcome)  != 0]
          gg_df_r= data.frame(pred = fitted(binary_model, 
                                             type = "response"),
                               classified = ifelse(fitted(binary_model,
                                                          type = "response") >= .5, 1,0),
                               real = df$outcome) #U can change it
          # ROC
          roc_plot = ggplot2::ggplot(gg_df_r, aes(d = real, m = pred)) + 
            plotROC::geom_roc(n.cuts = 10,
                              labels = FALSE,
                              size =2)+
            plotROC::style_roc()+
            ggplot2::geom_abline(slope = 1, 
                                 intercept = 0,
                                 size = 1,
                                 linetype = "dashed",
                                 alpha = .4)+
            ggplot2::theme(panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           aspect.ratio = 1)+
            ggplot2::ggtitle("")
          roc_RISK6 = roc(roc_plot$data$real, roc_plot$data$pred)
          pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
                         col = "#2A9D8F",
                         legacy.axes = T,
                         lwd = 3,
                         ci=TRUE,
                         print.auc=F,
                         add = T,
                         of="se", 
                         specificities=seq(0, 1, .1), 
                         ci.type="shape", ci.col=alpha("#2A9D8F", 0.15),
                         asp = NA)
          
          pROC::plot.roc(roc_plot$data$real, roc_plot$data$pred,
                         col = "#2A9D8F",
                         print.auc.y = .15,
                         lwd = 2,
                         ci=TRUE,
                         print.auc=TRUE,
                         print.auc.pattern = "RISK6 = %.2f(%.2f-%.2f)",  # = Rounds AUC to 2 digits
                         add = T,
                         asp = NA)
  }
        
}
