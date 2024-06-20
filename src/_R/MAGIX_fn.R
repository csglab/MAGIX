

draw_design_mat_heatmap <- function(B, out_prefix, name){
  mid <- min(B) + abs(max(B) - min(B))/2
  col_fun <- colorRamp2( c( min(B), mid, max(B) ), c( "white", "red" ,"black" ) )
  this_height <- 2 + 0.10714 * nrow(B)
  
  pdf( paste0(out_prefix, name ), width = 5, height = 8 )
  draw( ComplexHeatmap::Heatmap( as.matrix(B),
                                 column_title = "",
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 show_row_names = TRUE,
                                 show_column_names = TRUE,
                                 col = col_fun,
                                 use_raster = TRUE,
                                 rect_gp = gpar(col = "grey", lwd = 0.2),
                                 column_names_gp = gpar(fontsize = 6),
                                 row_names_gp = gpar(fontsize = 6),
                                 border = TRUE,
                                 name = " " ))
  dev.off()
}

plot_estimated_library_sizes <- function(model, out_prefix, this.title){
  
  s <- data.frame(names = colnames( model$params$B ),
                  estimated_effective_library_size_factors = model$params$s)
    
  # make V1 an ordered factor
  s$names <- factor(s$names, levels = s$names)
    
  p1 <- ggplot( data = s, aes(x = names, y=estimated_effective_library_size_factors)) + 
                geom_point() + xlab("") + theme_light() + ggtitle(this.title) + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  ggsave( paste0(out_prefix, "estimated_effective_library_size_factors.pdf"), 
          p1, width = (0.4 + nrow(s)*0.1785 ), height = 6 )
    
  write.csv(x = s, row.names = FALSE,
          file = paste0(out_prefix, "estimated_effective_library_size_factors.csv"))
}

coefs_qq_plot <- function( coefs_row, coefs_df, out_prefix_qq ){
  
  # coefs_row <- batch_coeffs_combinations[1,]
  # coefs_df <- coefs
  # out_prefix_qq <- out_prefix_qq
  
  
  colnames(coefs_df) <- gsub(":", "\\.", colnames(coefs_df))
  
  coefs1_name <- gsub( pattern = ":", replacement = "_", x = coefs_row[1] )
  coefs2_name <- gsub( pattern = ":", replacement = "_", x = coefs_row[2] )
  
  pdf( paste0(out_prefix_qq, "qq_plot_", coefs1_name, "_", coefs2_name, ".pdf" ), 
       width = 7, height = 7 )

  smoothScatter( sort( coefs_df[,coefs_row[1]] ), sort( coefs_df[,coefs_row[2]] ),
                 nbin=500,bandwidth=0.01, transformation=function(x)x^0.15,
                 xlab = coefs_row[1], ylab = coefs_row[2], 
                 main = paste0( this.title, "\nQ-Q plots" ))
  abline(0,1)
  dev.off()
  
}



write_ROC_PR <- function( this.coef, merged, out_prefix_chip, this.title, macs50_num ){
  
  # this.coef <- array(batch_coeffs)[1]
  # merged <- merged
  # out_prefix_chip <- out_prefix_chip
  # this.title <- this.title
  # macs50_num <- opt$total_ChIP

  # cat(paste0( colnames(merged), "\n\n" ) )
  # cat(paste0( this.coef, "\n\n" ) )
  
  this.coef <- gsub("-", "\\.", this.coef)
  
  ks_test <- ks.test( x = merged[ merged$overlap_chipPeak == TRUE , this.coef ],
                      y = merged[ merged$overlap_chipPeak == FALSE, this.coef ],
                      alternative = "less" )
  
  # Calculate ROC and PR curves
  roc.res <- roc.curve( scores.class0=merged[, this.coef],
                        weights.class0=merged$overlap_chipPeak, curve=T)
  
  pr.res <- pr.curve( scores.class0=merged[, this.coef],
                      weights.class0=merged$overlap_chipPeak, curve=T)
  
  # Calculate Jaccard index as a function of SELEX score cutoff
  ordered_indices <- order( unlist(merged[,this.coef]), decreasing = T )
  
  P <- sum( merged$overlap_chipPeak )
  N <- sum( ! merged$overlap_chipPeak )
  
  TP <- cumsum( merged[ordered_indices,]$overlap_chipPeak )
  FP <- cumsum( 1-merged[ordered_indices,]$overlap_chipPeak )
  FN <- P - TP
  # jaccard <- TP/(TP+FP+FN)
  jaccard <- TP/(FP+as.integer(macs50_num))
  
  # I'm thinking something simpler. E.g., TP rate at FP rate=0.01
  # Or at a few FP rates; e.g., TP rate at FP rate=0.01, 0.02, 0.05, and 0.1
  TN <- N - FP
  
  sen_spec <- data.frame( hits = merged[ordered_indices,]$overlap_chipPeak, 
                          TP = TP, FP = FP, TN = TN, FN = FN  )
  
  sen_spec$FPR <- sen_spec$FP / N
  # sen_spec$FPR <- sen_spec$FP / ( sen_spec$FP + sen_spec$TN )
  
  sen_spec$TPR <- sen_spec$TP / P
  # sen_spec$TPR <- sen_spec$TP / ( sen_spec$TP + sen_spec$FN )
  
  
  cat( paste0( "Batch: ", this.coef, 
               ", PR AUC: ", pr.res$auc.integral, 
               ", ROC AUC: ", roc.res$auc, 
               ", max Jaccard: ", max(jaccard), 
               ", ks D: ", ks_test$statistic,
               ", ks p-value: ", ks_test$p.value,
               ", TPR_0_01: ", max( sen_spec[ sen_spec$FPR <= 0.01, "TPR" ] ),
               ", TPR_0_02: ", max( sen_spec[ sen_spec$FPR <= 0.02, "TPR" ] ),
               ", TPR_0_1: ", max( sen_spec[ sen_spec$FPR <= 0.1, "TPR" ] ),
               ", TPR_0_2: ", max( sen_spec[ sen_spec$FPR <= 0.2, "TPR" ] ),
               ", peaks_at_MACS_score_50: ", macs50_num,
               "\n" ) )
  
  
  
  pdf( paste0(out_prefix_chip, this.coef, "_ROC_curve.pdf" ), width = 7, height = 7 )
  plot( roc.res, main=paste0(this.title, "\nEnrichment in ", this.coef), color=F )
  abline(0,1)
  dev.off()

  pdf( paste0(out_prefix_chip, this.coef, "_PR_curve.pdf" ), width = 7, height = 7 )
  plot( pr.res, main=paste0(this.title, "\nEnrichment in ", this.coef), color=F )
  dev.off()

  pdf(file = paste0(out_prefix_chip, this.coef, "_jaccard.pdf"), width = 7, height = 7 )
  smoothScatter( unlist(merged[ordered_indices,this.coef]),
                 main = paste0(this.title, "\nEnrichment in ", this.coef),
                 jaccard,xlab="GHT-SELEX score cutoff",
                 nbin=1000,
                 bandwidth = c(0.01,0.001),
                 transformation = function(x)x^0.1)
  dev.off()
}



compare_to_ChIP <- function(){
  
  cat("Check if coefs and stats exist\n")
  if ( ! ( exists("coefs") & exists("stats")  ) ) {
    
    Y <- read.csv( file = paste0(out_prefix, "imputed_log_norm_counts_Y.csv"),
                   row.names = 1 )
    coefs <- read.csv( paste0( out_prefix, "coefs.txt"), sep = "\t", row.names = 1 )

    stats <- read.table(paste0( out_prefix, "peak_stats.txt"),sep="\t", header = TRUE)

    batch_coeffs <- colnames(coefs)
    batch_coeffs <- batch_coeffs[ !grepl(pattern = "Aggregate", x = batch_coeffs) &
                                  !grepl(pattern = "Intercept", x = batch_coeffs) &
                                  !grepl(pattern = "aggregate", x = batch_coeffs) &
                                  !grepl(pattern = "^Batch", x = batch_coeffs) ]
    batch_coeffs <- gsub(":", ".", batch_coeffs )
  }
  
  this.title <- paste0( "Target: ", opt$TF_label, ", step: ", 
                        opt$step, ", n = ", nrow(coefs) )
  
  
  cat("Compare to ChIP-seq\n")
  dir.create(paste0(opt$outdir, "/compare_to_ChIP"), showWarnings = TRUE)
  out_prefix_chip <- paste0(opt$outdir, "/compare_to_ChIP/target_", opt$TF_label,
                            "_step_", opt$step, "_")

  cat("Merge the GHT-SELEX coefficient table with the peak status table\n")
  # merge the GHT-SELEX coefficient table with the peak status table
  peak_stats <- fread(opt$chip_peaks,data.table=F)
  merged <- merge(stats,peak_stats,by="name")

  write.table(merged,paste0(out_prefix, "peaks_stats_merged_wChIP_seq.tab"),
              sep="\t",quote=F,row.names = F)
  
  # Mark the genomic bins that exactly overlap (zero distance) a ChIP-seq peak (MACS cutoff 50)
  merged$overlap_chipPeak <- merged$MACS2_score > 50 & merged$distance == 0
  exact_overlap <- paste0( "genomic bins that exactly overlap (zero distance) ",
                           "a ChIP-seq peak (MACS cutoff 50): ",
                           sum(merged$overlap_chipPeak) )

  cat(paste0("\n", exact_overlap, "\n"))

  cat("Write ROC and PR plots\n")
  # cat( array(batch_coeffs) )
  cat("\n\n\n")
  sapply( X = array(batch_coeffs),
          FUN = write_ROC_PR,
          merged = merged,
          out_prefix_chip = out_prefix_chip,
          this.title = this.title,
          macs50_num = opt$total_ChIP )
  
  cat("done write ROC and PR plots\n")
  ## if (opt$only_compare_to_chip) {
      
  ##   counts <- fread( opt$count_matrix, data.table = FALSE )
  ##   geneIDs <- counts[,4]
  ##   counts <- counts[,-(1:4)]
  ##   rownames(counts) <- geneIDs

  ##   design <- read.csv(file = paste0(out_prefix, "design.csv"), 
  ##                      header = TRUE, row.names = 1)
  ##   counts <- counts[,design$Sample]
  ##   }
    
  ## cat("\nVisualize log-normalized counts of bins that overlap ChIP-seq peaks\n")

  ## ### get the list of bins that overlap ChIP-seq peaks
  ## filter <- unique(merged$name[merged$overlap_chipPeak])

  ## ### This is simply the log-normalized read count (plus a pseudocount)
  ## logM_normalized <- as.matrix(log( t(t(counts[filter,]+0.01)/colSums(counts)) ))

  ## top_annotations <- HeatmapAnnotation( df = design[,c("Batch", "Target", "Cycle")] )

  ## mid <- min(logM_normalized) + abs(max(logM_normalized) - min(logM_normalized))/2
  ## col_fun <- colorRamp2( c( min(logM_normalized), mid, max(logM_normalized) ),
  ##                        c( "white", "orange" ,"red" ) )

  ## pdf( paste0(out_prefix, "logM_normalized_counts.pdf" ), width = 12, height = 8 )
  ## draw( ComplexHeatmap::Heatmap( as.matrix(logM_normalized),
  ##                                column_title = paste0(this.title,
  ##                                                      ", log-normalized read count\nShowing ",
  ##                                                      exact_overlap),
  ##                                cluster_rows = TRUE,
  ##                                cluster_columns = FALSE,
  ##                                show_row_names = FALSE,
  ##                                top_annotation = top_annotations,
  ##                                col = col_fun,
  ##                                use_raster = TRUE,
  ##                                name = "logM_norm",
  ##                                column_names_gp = gpar(fontsize = 6) ) )
  ## dev.off()


  ## mid <- min(Y[filter,]) + abs(max(Y[filter,]) - min(Y[filter,]))/2
  ## col_fun <- colorRamp2( c( min(Y[filter,]), mid, max(Y[filter,]) ),
  ##                        c( "white", "orange" ,"red" ) )

  ## pdf( paste0(out_prefix, "imputed_log_normalized_counts.pdf" ), width = 12, height = 8 )
  ## draw( ComplexHeatmap::Heatmap( as.matrix( Y[filter,] ),
  ##                                column_title = paste0(this.title,
  ##                                                      ", imputed log-normalized counts\nShowing ",
  ##                                                      exact_overlap),
  ##                                cluster_rows = TRUE,
  ##                                cluster_columns = FALSE,
  ##                                show_row_names = FALSE,
  ##                                top_annotation = top_annotations,
  ##                                col = col_fun,
  ##                                name = "imputed\nlogM_norm",
  ##                                use_raster = TRUE,
  ##                                column_names_gp = gpar(fontsize = 6) ) )
  ## dev.off()
    
  ## write.csv(Y[filter,], paste0(out_prefix, "imputed_log_normalized_counts_in_heatmap.csv"))
  ## write.csv(logM_normalized, paste0(out_prefix, "logM_normalized_counts_in_heatmap.csv"))
}


