suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(ggplot2)
library(tictoc)
library(Matrix)
library(uwot)
# library(pheatmap)
library(stringr)
library(data.table)
library(PRROC)
# Loading libraries and scripts that are needed for core functionality
# message("Loading required libraries and scripts...")
library(sads)
library(Rcpp)
library(optparse)
set.seed(1)

########################################################   IN and load data ####
option_list = list(
  make_option(c("-a", "--script_path"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX",
              help=""),

  make_option(c("-b", "--outdir"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/OUT",
              help=""),

  make_option(c("-c", "--step"), type="character", metavar="integer", default=4, help="" ),
  
  make_option(c("-d", "--count_matrix"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/IN/GHT_SELEX_CTCF_counts_aggregates_and_experiments_no_blacklisted_100k.bed",
              help=""),
  
  make_option(c("-e", "--design"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/IN/CTCF_design_matrix_per_TF.txt",
              help=""),
  
  make_option(c("-f", "--TF"), type="character", metavar="character", default="CTCF", help="Should be the same as in the design matrix"),
  
  make_option(c("-g", "--TF_label"), type="character", metavar="character", default="CTCF_FL", help="Used for the output file prefix"),  
  
  ## Step 2 exclusive options
  make_option(c("-i", "--step2_batches"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/IN/CTCF_experiments_per_TF.txt",
              help="Either a file to a list with the selected experiment IDs or the string all"),
  
  make_option(c("-j", "--lrt_sample_size"), type="character", metavar="character",
              default="default_none", help=""),
  
  ## Compare to ChIP-seq exclusive options
  make_option(c("-k", "--use_chip"), type="character", metavar="character",
              default="FALSE", help=""),
  
  make_option(c("-l", "--chip_peaks"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/IN/whole_genome_200_bins_closest_summit_CTCF_CTCFChIP2_IGO_10521_25_S37_small_100k.tab",
              help=""),
  
  make_option(c("-m", "--only_compare_to_chip"), type="character", metavar="character",
              default="FALSE", help=""),
  ## Covariance account_covariance
  make_option(c("-n", "--account_covariance"), type="logical", default=TRUE, help=""),
  
  make_option(c("-o", "--previous_model"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/MAGIX/data/CTCF_demo/OUT/MAGIX_CTCF_FL/target_CTCF_FL_step_2_model.RDS",
              help=""),

  make_option(c("-p", "--total_ChIP"), type="character", metavar="character",
              default="100",
              help=""),  
    
  make_option(c("-q", "--test_depletion"), type="logical", metavar="character",
              default=FALSE, help=""),
  
  make_option(c("-r", "--remove_target_from_aggregate"), type="logical", metavar="character",
              default=FALSE, help="")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
################################################################################

opt$outdir <- paste0(opt$outdir, "/MAGIX_", opt$TF_label )
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
sourceCpp(paste0(opt$script_path, "/src/_C++/Eigen.v8.cpp"))
source(paste0(opt$script_path, "/src/_R/MAGIX_fn.R"))
source(paste0(opt$script_path, "/src/_R/logNormPois_Regression_covariance_mask_v13.R"))


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
  cat("\nTesting for depletion instead of enrichment.\n")
  design$Cycle <- -design$Cycle
  }


this.title <- paste0( "Target: ", opt$TF_label, ", step: ", opt$step, ", n = ", nrow(counts) )

cat(this.title)
################################################################################



############################# Substract target counts from aggregates counts ####
# opt$remove_target_from_aggregate <- TRUE
if ( as.logical(opt$remove_target_from_aggregate) ) {
  
  for (this_batch in unique( design$Batch ) ) {
    
    # this_batch <- "YWE_B"
    
    for ( this_cycle in unique( design$Cycle ) ) {
      
      # this_cycle <- 1
      this_samples <- design[ ( design$Batch == this_batch ) & 
                              ( design$Cycle == this_cycle ) , ]
      
      aggregate_names <- this_samples[ this_samples$Target == "Aggregate", "Sample" ]
      target_names <- this_samples[ this_samples$Target != "Aggregate", "Sample" ]
      
      if(length( aggregate_names ) > 1 ){
        cat("you should have one aggregate per batch/cycle.\n"); q()
      }
      
      if (length(target_names) > 1 ){
        
        tmp_col <- data.frame( rowSums( counts[, target_names] ) )
        counts[, aggregate_names] <- counts[, aggregate_names] - tmp_col
        rm(tmp_col)
      }
      else{
        counts[, aggregate_names] <- counts[, aggregate_names] - counts[, target_names]
      }
    }
  }
}

### Check if removing target from aggregate counts creates negative counts.
## This can happen when you didnt create the aggregate using the targets and
## specified --remove_target_from_aggregate = TRUE.
if ( sum( counts < 0 ) != 0) { cat("Negative values in count matrix.\n") }
# rm(aggregate_names, this_batch, this_samples, target_names )

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
if( opt$step == 2 | opt$step == 3 ){
    
    if( opt$step2_batches == "all"){
        selected_samples <- design$Experiment_ID
    }else{
        selected_samples <- read.csv(file = opt$step2_batches, header = FALSE)
        selected_samples <- as.vector(selected_samples$V1)
        selected_samples <- gsub("-", ".", selected_samples )
    }
    
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
if( opt$step == 2 | opt$step == 3 ){

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

if (opt$step == 3 ) {
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
if( opt$step %in% c("2", "3")) {
  
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
if( opt$step == 3 ){
  cat("\nLRT\n")
  # model <- model2
  ## Find TargetCTCF:Cycle ( TargetCTCF:Cycle ) column number
  target_column_name <- paste0( "Target", opt$TF, ":Cycle" )
  rm_Index <- which( colnames(B) == target_column_name )
  
  # The above step only fits the model parameters (i.e. identifies maximum a-priori
  # values of the model coefficients). In order to perform statistical analysis,
  # we can use likelihood ratio test (LRT):

  cat("\nLikelihood ratio test\n")
  if( (opt$lrt_sample_size == "default_none") | 
      ( as.numeric(opt$lrt_sample_size) > nrow(counts) ) ){ 
    
    opt$lrt_sample_size <- nrow(counts) 
    }
  
  cat( "LRT sample size: ", as.numeric(opt$lrt_sample_size), "\n" )
  
  lrt.res <- model$lrt(
    rmIndex = rm_Index, # This is the index of the variable (in the model matrix) that we want to test
    topEntries = as.numeric(opt$lrt_sample_size), # Here, I'm performing LRT for only 10,000 bins (out of 2,000,000) due to speed
    iterations = 10000, track_internval = 20, etol = 1e-3 )
  cat("\n")
  
  lrt.title <- paste0( "Target: ", opt$TF_label, ", step: ", opt$step,
                       ", LTR n = ", opt$lrt_sample_size )

  
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















