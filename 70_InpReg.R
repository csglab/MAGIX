library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tictoc)
library(Matrix)
library(uwot)
library(pheatmap)
library(stringr)
library(data.table)
library(PRROC)
# Loading libraries and scripts that are needed for core functionality
# message("Loading required libraries and scripts...")
library(sads)
library(Rcpp)
library(optparse)
set.seed(1)

## ahcorcha
## Add metric logs from ./ref/compare_to_ChIP_ROC_PR.R
############################################################### Functions ######
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
  # this.coef <- "BatchYWE_B_AffSeq_H12_CTCF.Cycle"
  # this.coef <- "Experiment_IDYWI_B_AffSeq_C10_ZNF888_DBD.Cycle"
  
  ks_test <- ks.test( merged[ merged$overlap_chipPeak == TRUE , this.coef ],
                      merged[ merged$overlap_chipPeak == FALSE, this.coef ],
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
  
  ## pdf( paste0(out_prefix_chip, this.coef, "_ROC_curve.pdf" ), width = 7, height = 7 )
  ## plot( roc.res, main=paste0(this.title, "\nEnrichment in ", this.coef), color=F )
  ## abline(0,1)
  ## dev.off()
  
  ## pdf( paste0(out_prefix_chip, this.coef, "_PR_curve.pdf" ), width = 7, height = 7 )
  ## plot( pr.res, main=paste0(this.title, "\nEnrichment in ", this.coef), color=F )
  ## dev.off()
  
  ## pdf(file = paste0(out_prefix_chip, this.coef, "_jaccard.pdf"), width = 7, height = 7 )
  ## smoothScatter( unlist(merged[ordered_indices,this.coef]),
  ##                main = paste0(this.title, "\nEnrichment in ", this.coef),
  ##                jaccard,xlab="GHT-SELEX score cutoff",
  ##                nbin=1000,
  ##                bandwidth = c(0.01,0.001),
  ##                transformation = function(x)x^0.1)
  ## dev.off()
}



compare_to_ChIP <- function(){
  
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
  
  cat("\nCompare to ChIP-seq\n")
  dir.create(paste0(opt$outdir, "/compare_to_ChIP"), showWarnings = TRUE)
  out_prefix_chip <- paste0(opt$outdir, "/compare_to_ChIP/target_", opt$TF_label,
                            "_step_", opt$step, "_")

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

  
  sapply( X = array(batch_coeffs),
          FUN = write_ROC_PR,
          merged = merged,
          out_prefix_chip = out_prefix_chip,
          this.title = this.title,
          macs50_num = opt$total_ChIP )
  

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





################################################################################


########################################################   IN and load data ####
# /home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/CTCF_in/GHT_SELEX_CTCF_counts_aggregates_and_experiments_no_blacklisted_100k.bed"

option_list = list(
  make_option(c("-z", "--script_path"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/03_process_GHT_SELEX/",
              help=""),

  make_option(c("-o", "--outdir"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_CTCF",
              help=""),

  make_option(c("-s", "--step"), type="character", metavar="integer", default=4, help="" ),
  
  make_option(c("-m", "--count_matrix"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_CTCF/GHT_SELEX_CTCF_counts_aggregates_and_experiments_no_blacklisted_100k.bed",
              help=""),
  
  make_option(c("-d", "--design"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_CTCF/CTCF_design_matrix_per_TF.txt",
              help=""),
  
  make_option(c("-t", "--TF"), type="character", metavar="character", default="CTCF", help="Should be the same as in the design matrix"),
  
  make_option(c("-k", "--TF_label"), type="character", metavar="character", default="CTCF_FL", help="Used for the output file prefix"),  
  
  ## Step 2 exclusive options
  make_option(c("-b", "--step2_batches"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_CTCF/CTCF_good_experiments.txt",
              help=""),
  
  make_option(c("-x", "--ltr_sample_size"), type="character", metavar="integer",
              default=70761, help=""),
  
  ## Compare to ChIP-seq exclusive options
  make_option(c("-c", "--use_chip"), type="character", metavar="character",
              default="TRUE", help=""),
  
  make_option(c("-p", "--chip_peaks"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_PRDM10/whole_genome_200_bins_closest_summit_PRDM10_PRDM10ChIP2_IGO_10521_46_S110_small.tab",
              help=""),
  
  make_option(c("-l", "--only_compare_to_chip"), type="character", metavar="character",
              default="FALSE", help=""),
  ## Covariance account_covariance
  make_option(c("-a", "--account_covariance"), type="logical", default=TRUE, help=""),
  
  make_option(c("-n", "--previous_model"), type="character", metavar="character",
              default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/V13_Inp_Reg_results/InpReg_CTCF/target_CTCF_FL_step_2_model.RDS",
              help=""),

  make_option(c("-q", "--total_ChIP"), type="character", metavar="character",
              default="100",
              help=""),  
    
  make_option(c("-e", "--test_depletion"), type="logical", metavar="character",
              default=FALSE, help="")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
################################################################################

dir.create(opt$outdir, showWarnings = TRUE, recursive = TRUE )
out_prefix <- paste0(opt$outdir, "/target_", opt$TF_label, "_step_", opt$step, "_")

if(opt$only_compare_to_chip == "TRUE" ){
  
  opt$use_chip <- TRUE
  sink(file=paste0(out_prefix, "only_compare_to_ChIP_log.txt"))
  compare_to_ChIP()
  sink()
  q() }

sink(file=paste0(out_prefix, "log.txt"))

# Load the source codes for the old and the new GEDI
# sourceCpp(paste0(opt$script_path, "/src/Eigen.v7.cpp"))
sourceCpp(paste0(opt$script_path, "/src/Eigen.v8.cpp"))

# source(paste0(opt$script_path, "/src/logNormPois.Regression.v12.shrinkage_with_covariance.R"))
source(paste0(opt$script_path, "/src/logNormPois.Regression.v13.covariance_mask.R"))


################################################################# Load data ####
counts <- fread( opt$count_matrix, data.table = FALSE )
geneIDs <- counts[,4]
counts <- counts[,-(1:4)]
rownames(counts) <- geneIDs

design <- fread( opt$design, data.table = FALSE )
rownames(design) <- design$Sample

design <- design[ order(design$Batch),]
design <- design[ order(design$Cycle),]
design <- design[ order(design$Experiment_ID),]
design <- design[ order(design$Target),]

if ( opt$test_depletion ) { 
  cat("\nDepletion\n")
  design$Cycle <- -design$Cycle
  }


this.title <- paste0( "Target: ", opt$TF_label, ", step: ", opt$step, ", n = ", nrow(counts) )

cat(this.title)
################################################################################



###################################################### Setup and run lnpSVD ####
# In step 1, we will simply analyze all batches for this TF
if( opt$step == 1 ){ 
  batch <- unique( design$Batch[design$Target == opt$TF] ) 
  
  # Only keep (a) aggregate counts and (b) counts that correspond to the TF/batch of interest
  filter <- which( (design$Target=="Aggregate" & design$Batch %in% batch ) | 
                    ( design$Target != "Aggregate" & design$Batch %in% batch ) )
  }


# In step 2, we will limit to samples that are "good"
if( opt$step == 2 | opt$step == 3 | opt$step == 4 ){
  selected_samples <- read.csv(file = opt$step2_batches, header = FALSE)
  selected_samples <- as.vector(selected_samples$V1)
  selected_samples <- gsub("-", ".", selected_samples )
  
  batch <- unique(design[design$Experiment_ID %in% selected_samples, "Batch"])
  # design <- design[ design$Batch %in% exp_batches, ]
  # batch <- unique( design$Batch )
  
  # Only keep (a) aggregate counts and (b) counts that correspond to the TF/batch of interest
  filter <- which( ( design$Target == "Aggregate" & design$Batch %in% batch ) | 
                   ( design$Target != "Aggregate" & design$Experiment_ID %in% selected_samples ) )
  }

design <- design[filter,] # filter the design matrix

if( !all(design$Sample %in% colnames(counts))) {
  stop("Some samples in the design matrix do not have a match in the provided count table.")
  }

# now only keep columns in the count matrix that correspond to remaining samples in the design data frame
counts <- counts[,design$Sample]
# convert the count data frame to a sparse matrix
counts <- as(as.matrix(counts),"dgCMatrix")
num_bin <- as.integer(nrow(counts))

# For the aggregate cases, replace the Target with the Experiment, since in any
#  situation (even in step 2) we will consider each batch separately for the aggregates
design$Target[design$Target=="Aggregate"] <- design$Experiment_ID[design$Target=="Aggregate"]

write.csv( x = design, file = paste0(out_prefix, "design.csv") )

# create the model matrix that should be fitted to the data
### In step 1, each batch will get its own "slope"
if( opt$step == 1 ){

  if( length(unique(design$Batch)) >= 2 ){
    B <- model.matrix(~Batch+Experiment_ID:Cycle+0,design)
    } else{
       B <- model.matrix(~Experiment_ID:Cycle+1,design)
       }
  }

### In step 2, the good batches are combined, so that we are left with one
###  slope for the TF of interest (and one slope for the aggregate signal in each batch)
if( opt$step == 2 | opt$step == 3 | opt$step == 4 ){

  if( length(unique(design$Batch)) >= 2 ){
    B <- model.matrix(~Batch+Target:Cycle+0,design)
    } else{
       B <- model.matrix(~Target:Cycle+1,design)
       }
  }

rownames(B) <- design$Sample

write.csv(x = as.data.frame(B), file = paste0(out_prefix, "design_matrix.csv"))

draw_design_mat_heatmap(B = B, out_prefix = out_prefix, name = "design_matrix_heatmap.pdf" )

################################################################################


############################################################ Create mask    ####
# First, create a mask that shows which coefficients are in the same group
cat("\nCreate mask\n")
groups <- rep(0,ncol(B)) # 0 means intercept
groups[grep("aggregate",colnames(B))] <- 1 # the slope of enrichment for aggregates
groups[grep( opt$TF,colnames(B))] <- 2 # the slope of enrichment for the TF of interest
# for aesthetic purposes, sort the columns of B based on the groups
B <- B[,order(groups)]
groups <- groups[order(groups)]

draw_design_mat_heatmap(B = B, out_prefix = out_prefix,
                        name = "design_matrix_with_groups_heatmap.pdf" )

write.csv( x = data.frame( coefficient=cbind(colnames(B),group=groups) ), row.names = F,
           file = paste0( out_prefix, "coefficient_groups.csv" ) )

cov_mask <- sapply(groups,function(x)x==groups)*1
# prmatrix(cov_mask,rowlab=rep("",nrow(cov_mask)),collab=rep("",ncol(cov_mask)))


################################################################################


######################################################## Fitting the models ####
# First, we need to create an object for the model
cat("\nCreating an object for the model\n")
model <- new("lnpReg")

# Next, the model should be set up with the count matrix and the design matrix.
# Here, due to computational limitations of my laptop,
# I'm only using a random subset of the count matrix rows. In the actual application,
# either all rows should be used (if memory allows),
#  or the count matrix should be split into chunks
# (preferably after shuffling the rows to make sure each chunk contains
# genomic regions from various chromosomes), the code should be applied to
# each chunk separately, and then the results should be combined
# model$setup( counts[sample.int(nrow(counts),200000),], B )
cat("\nFit model parameters\n")

if( opt$account_covariance == FALSE ){
  model$setup( counts, B, cov_mask=NULL )
  } else if ( opt$account_covariance == TRUE ) {
  model$setup( counts, B, cov_mask )
} else{
  stop(paste0( "account_covariance parameter should be TRUE or FALSE is: ",
               opt$account_covariance )) }

cat("\n")

if (opt$step == 3 | opt$step == 4) {
  cat("\nSetting parameters from model in step 2\n")
  previous_model <- readRDS( opt$previous_model )
  model$params$s <- previous_model$params$s
  model$params$sigma2 <- previous_model$params$sigma2
  model$aux$fixed_global_params <- TRUE 
  model$hyperparams$Z_precision <- previous_model$hyperparams$Z_precision
  }


# The next step is to fit the model. Here, I'm using 100 iterations, but we
# may want to increase the iterations if the model does not converge
cat("\nOptimization\n")
tic("Optimization")
model$optimize(iterations = 100, track_internval=10, etol=1e-3)
toc()

# This function returns the "imputed" log-normalized counts.
# Imputation is done by calculating the *expected* log-normalized count of each
# genomic bin in each experiment, given the *observed* count of that bin in that
# experiment and the *fitted* count from the model.
cat("\nGet Y\n")
Y <- getY.lnpReg(model)
write.csv( x = Y, row.names = TRUE,
           file = paste0(out_prefix, "imputed_log_norm_counts_Y.csv") )

####### Visualize results and identify good batches
cat("\nVisualize estimated effective library size factor\n")
## Here, I am just plotting the estimated effective library size factors. The only
## pattern that I expect at this point is a decrease in library size as the cycles progress.
plot_estimated_library_sizes(model, out_prefix, this.title)


# Plot the pairwise correlation of the model coefficients. We expect the Intercept
# and the Aggregate coefficients to have low correlation with the other coefficients,
# while the others should form a cluster. If a particular batch is an outlier,
# it can be detected in this plot. If there’s other evidence that this batch has
# not worked, it can be removed in step 2.
#
# Note that, if genomic bins were analyzed in several chunks, the model coefficients
# of these chunks should be first combined prior to the downstream analyses shown in this notebook.


coefs <- getZ.lnpReg(model)
stats <- data.frame( name=rownames(coefs), coefs )


col_fun <- colorRamp2( c( -1, 0, 1 ), c( "blue", "white" ,"red" ) )

pdf( paste0(out_prefix, "coefs_heatmap.pdf" ), width = 10, height = 10 )
draw( ComplexHeatmap::Heatmap(as.matrix(cor(coefs)), col = col_fun, column_title = this.title,
                              column_names_gp = gpar(fontsize = 8),
                              row_names_gp = gpar(fontsize = 8) ) )
dev.off()

batch_coeffs <- colnames(coefs)
batch_coeffs <- batch_coeffs[ !grepl(pattern = "Aggregate", x = batch_coeffs) &
                              !grepl(pattern = "Intercept", x = batch_coeffs) &
                              !grepl(pattern = "aggregate", x = batch_coeffs) &
                              !grepl(pattern = "^Batch", x = batch_coeffs) ]
batch_coeffs <- gsub( ":",".", batch_coeffs )

if ( opt$step == 1 & length(batch_coeffs) >= 2 ){

    cat("\nQ-Q plots\n")
    dir.create(paste0(opt$outdir, "/qq_plots"), showWarnings = TRUE)
    out_prefix_qq <- paste0(opt$outdir, "/qq_plots/target_", opt$TF_label, "_step_", opt$step, "_")

    batch_coeffs_combinations <- as.data.frame(t(combn( batch_coeffs, m = 2 )))

    coefs_qq_plot <- apply( X=batch_coeffs_combinations,
                            MARGIN = 1,
                            FUN=coefs_qq_plot, coefs_df=coefs,
                            out_prefix_qq = out_prefix_qq )
    }

write.table(stats, paste0( out_prefix, "peak_stats.txt"),sep="\t",quote=F,row.names = F)
write.table(coefs,paste0( out_prefix, "coefs.txt"),sep="\t",quote=TRUE)


############################################################################# ##
### Create coefficient bed file
if( opt$step %in% c("2", "3", "4")) {
  
  coeff_bed <- str_split_fixed( stats$name, ":", 2 )
  coeff_bed <- cbind( coeff_bed[,1], str_split_fixed( coeff_bed[,2], "-", 2 ) )
  
  coeff_bed <- cbind( coeff_bed, 
                      stats[, paste0("Target", opt$TF, ".Cycle" ) ] )
  
   colnames(coeff_bed) <- c("chr", "start", "stop", 
                            paste0("Target", opt$TF, ".Cycle" ) )
  
  write.table(coeff_bed, paste0( out_prefix, "coefficient_br.bed"), 
              sep = "\t", quote = F, row.names = F)
}


 
####################################################################### LRT ####
if( opt$step == 4 ){
  cat("\nLRT\n")
  # model <- model2
  ## Find TargetCTCF:Cycle ( TargetCTCF:Cycle ) column number
  target_column_name <- paste0( "Target", opt$TF, ":Cycle" )
  rm_Index <- which( colnames(B) == target_column_name )
  
  
  # The above step only fits the model parameters (i.e. identifies maximum a-priori
  # values of the model coefficients). In order to perform statistical analysis,
  # we can use likelihood ratio test (LRT):

  cat("\nLikelihood ratio test\n")
  
  # if( opt$step == 3 ){
  #   lrt.res <- model$lrt(
  #   rmIndex = rm_Index, # This is the index of the variable (in the model matrix) that we want to test
  #   topEntries = as.numeric(opt$ltr_sample_size), # Here, I'm performing LRT for only 10,000 bins (out of 2,000,000) due to speed
  #   iterations = 10000, track_internval = 20, etol = 1e-3 )
  # }
  
  if( opt$step == 4 ){
    lrt.res <- model$lrt(
    rmIndex = rm_Index, # This is the index of the variable (in the model matrix) that we want to test
    topEntries = nrow(counts), # Here, I'm performing LRT for only 10,000 bins (out of 2,000,000) due to speed
    iterations = 10000, track_internval = 20, etol = 1e-3 )
  }
  
    cat("\n")
  
  lrt.title <- paste0( "Target: ", opt$TF_label, ", step: ", opt$step,
                       ", LTR n = ", opt$ltr_sample_size )

  
  # lrt.res[order(lrt.res$pvalue),]
  # The LRT function returns a data frame with multiple columns:
  # - coefficient.br: This is the coefficient of the variable of interest in the full model,
  #                   as fitted before calling the LRT function.
  # - coefficient.ar: The LRT function further "refines" the coefficients of the full model
  #                   to ensure all have converged. Furthermore, unlike the "optimize" function,
  #                   the LRT function does not perform coefficient shrinkage.
  # - full_LL: log-likelihood of the full model
  # - reduced_LL: log-likelihood of the reduced model (i.e. the model in which
  #               the variable of interest is removed)
  # - pvalue: LRT p-value
  # - fdr: FDR-adjusted p-value


  # Plot a volcano plot based on the LRT. Note that only the top 10,000 candidates
  # that were tested can be plotted, skewing the plot toward positive coefficients

  # Volcano plot using the before-refinement coefficient values
  pdf( paste0(out_prefix, "LRT_before_refinement_coefficients.pdf" ), width = 7, height = 7 )
    smoothScatter(lrt.res$coefficient.br,-log10(lrt.res$pvalue),
                  main=paste0(lrt.title, "\nVolcano plot of before-refinement coefficient values" ),
                 xlab = "before-refinement coefficients", ylab = "-log10(p-value)" )
  dev.off()

  # Volcano plot using the after-refinement coefficient values
  pdf( paste0(out_prefix, "LRT_after_refinement_coefficients.pdf" ), width = 7, height = 7 )
  smoothScatter(lrt.res$coefficient.ar,-log10(lrt.res$pvalue),
                main=paste0(lrt.title, "\nVolcano plot of after-refinement coefficient values"),
                xlab = "after-refinement coefficients", ylab = "-log10(p-value)" )
  dev.off()

  # Directly compare the before-refinement and after-refinement coefficient values
  # Note that after-refinement coefficients are generally expected to be larger, because there is no shrinkage
  pdf( paste0(out_prefix, "LRT_compare_before_and_after_refinement_coefficients.pdf" ), width = 7, height = 7 )
  smoothScatter(lrt.res$coefficient.br,lrt.res$coefficient.ar,
                main=paste0(lrt.title, "\nCompare before and after-refinement coefficient values"),
                xlab = "before-refinement", ylab = "after-refinement" )
  abline(0,1)
  dev.off()

  write.csv(lrt.res, paste0(out_prefix, "LTR_results.csv"), row.names = FALSE)
  }



######################################################  Compare to ChIP-seq ####
if ( opt$use_chip == "TRUE"){ 
  cat("\nCompare to ChIP-seq\n")
  compare_to_ChIP() }

######################################################################  End ####
# Save the model as RDS object. But first, I'm going to erase some of the internal
# objects in the model that are quite large but are not needed to be saved.
model$target$Y <- NULL
model$target$M <- NULL
model$aux$ZB <- NULL
saveRDS(model, paste0( out_prefix, "model.RDS") )
cat("\n\n\n\n")
sessionInfo()
sink()















