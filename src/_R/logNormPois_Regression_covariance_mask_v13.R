
setRefClass(
  
  "lnpReg",
  
  fields = list(
    params="list", # This is where the main parameters of the model are stored
    hyperparams="list", # This is where the hyper-parameters of the model are stored
    aux="list", # This is where the auxiliary variables that are useful for model fitting are stored
    target="list", # This is where the target matrices (Y and/or M) are stored
    convergence="list" # This is where the convergence information for model fitting progress is stored
  ),
  
  methods = list(
    
    #######################
    # Function to create the matrices and prepare other variables for a lnpReg model
    #######################
    
    setup = function(
      M, # The raw read count matrix
      B, # The design matrix
      cov_mask = NULL # This mask will be applied to the empirical prior for covariance
    ) {
      
      message("Setting up the lnpReg model...")
      
      ##### initialize the observation type

      # Note: M will be stored in the same format as the one provided (e.g. dgCMatrix)
      #  However, Y is always stored as `matrix`
      if( is.null(M) ) { stop("Error: M should be provided") }
      if( is.null(B) ) { stop("Error: B should be provided") }
      if( ncol(B) != ncol(M) ) {
        message("B and M have different number of columns. Checking if B should be transposed ...")
        if( nrow(B) == ncol(M) ) {
          B <- t(B)
        } else {
          stop("Error: B and M have incompatible dimensions")
        }
      }
      if( any( colnames(B) != colnames(M) ) ) {
        message("Warning: M and B column names are not the same. Column names of M will be used.")
      }
      
      Y <- as.matrix( log(0.01+M) ) # if M is provided as a single matrix, Y will be initialized (or over-written) as log of M
      aux$obs.type <<- "M"
      aux$geneIDs <<- rownames(M)
      aux$sampleIDs <<- colnames(M)
      aux$fixed_global_params <<- F # If T, then sigma2 and s will not be updated
      aux$optimize_all <<- T # If F, then only the Z rows that are marked in the next vector will be updated
      aux$optimize_indices <<- 1:nrow(Y) # Will be used only if aux$optimize_all is F
      aux$shrinkage <<- T # By default, shrinkage is on. It will be off only when performing LRT

      ##### Store the model dimensions
      
      aux$J <<- nrow(Y)
      aux$J_vec <<- rep(1,aux$J)
      aux$N <<- ncol(Y)
      aux$N_vec <<- rep(1,aux$N)
      aux$K <<- nrow(B)
      aux$diag_K <<- diag(aux$K)

      # if a mask is provided for the empirical prior covariance matrix, store it
      if( !is.null(cov_mask) ) {
        aux$cov_mask <<- cov_mask
      } else { # otherwise, use only the variances (and not the covariances)
        aux$cov_mask <<- aux$diag_K
      }
      
      # store additional auxiliary information
      hyperparams$Z_precision <<- matrix(0,nrow=aux$K,ncol=aux$K)
      
      aux$ite <<- 0 # the number of optimization iterations completed so far
      
      ##### initialize global parameters
      
      # initialize the sample-level scaling factors
      s_0 <- eigenVecMatProduct( aux$J_vec, as.matrix(M) ) / aux$J # s is initialized as the column-mean of M
      params$s <<- log(s_0) # convert s_0 to log-scale, which is compatible with Y

      # Initialize sigma2
      Yp <- Y - eigenVecVecProduct( aux$J_vec, params$s )
      params$sigma2 <<- rep(
        sum(c(Yp*Yp)) / (aux$N*aux$J+1),
        aux$N ) # Initially, sigma2 has the same value for all samples
      
      ##### Initialize other model parameters and matrices
      
      # Initialize Z and B
      params$Z <<- matrix(0,nrow=aux$J,ncol=aux$K)
      params$B <<- B

      # Initialize target matrices
      target$M <<- M
      target$Y <<- Y


      # Initialize convergence vectors
      convergence$dZ <<- NA # This vector tracks the change of Z in the last iteration; each element corresponds to one gene
      convergence$dY <<- NA # This vector tracks the change of Y in the last iteration; each element corresponds to one gene
      convergence$sigma2 <<- params$sigma2 # This vector tracks the value of sigma2 across iterations
      convergence$ds <<- NA # This vector tracks the change of s in the last iteration; each element corresponds to one sample
    },

    #######################
    # Functions for fitting individual components of the model (by block coordinate optimization)
    #######################
  

    ######
    # Solves for Z, the shared metagene matrix: Yi ~ o + oi + (Z+Qi) x Bi
    # This function enforces orthogonality of Z, by finding an orthonormal U and diagonal S so that Z=US
    # Note that this function requires QiDBi to be up-to-date
    solve.Z = function() {
      multipliers <- sqrt(1/params$sigma2)
      if( !aux$shrinkage ) { # no shrinkage is required
        diagLambda <- matrix(0,nrow=aux$K,ncol=aux$K)
      } else { # shrinkage using a prior distribution for Z~N(0,sigma2*precision^-1)
        diagLambda <- hyperparams$Z_precision
      }
      
      if( aux$optimize_all ) {
        params$Z <<- solveZ(
          eigenRowMult( target$Y, multipliers ),
          params$s * multipliers, rep(0,aux$J),
          eigenRowMult( params$B, multipliers ),
          aux$diag_K,
          diagLambda )
      } else {
        params$Z[aux$optimize_indices,] <<- solveZ(
          eigenRowMult( target$Y[aux$optimize_indices,,drop=F], multipliers ),
          params$s * multipliers, rep(0,length(aux$optimize_indices)),
          eigenRowMult( params$B, multipliers ),
          aux$diag_K,
          diagLambda )
      }
      
      if( !aux$fixed_global_params ) { # This is not a reduced model, and therefore s can be updated
        # If any column of Z has a non-zero mean, shift it to s
        for( i in 1:aux$K ) {
          meanZ <- mean(params$Z[,i])
          params$Z[,i] <<- params$Z[,i] - meanZ
          params$s <<- params$s + meanZ*params$B[i,]
        }
      }
      
    },

    ######
    # Update ZB
    update.ZB = function() {
      if( aux$optimize_all ) {
        aux$ZB <<- eigenMatProduct(
          params$Z,params$B )
      } else {
        aux$ZB[aux$optimize_indices,] <<- eigenMatProduct(
          params$Z[aux$optimize_indices,,drop=F],params$B )
      }
    },

    ######
    # Solves for si, the cell-specific library size
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZBi and QiDBi to be up-to-date
    solve.s = function()
    {
      if( !aux$fixed_global_params ) { # This is not a reduced model, and therefore s can be updated
        params$s <<- solveS(
          target$Y,
          aux$ZB,
          aux$J_vec,
          rep(0,aux$J),
          aux$J )
      }
    },
    
    
    ######### Up to here
    ######
    # Solve for sigma2, the mean squared error of the model, for the case where M is observed
    # Note that this function requires ZBi and QiDBi to be up-to-date for all samples
    solve.sigma2 = function()
    {
      if( !aux$fixed_global_params ) { # This is not a reduced model, and therefore sigma2 can be updated
        params$sigma2 <<- solveSigma2_M(
          target$Y,
          aux$ZB,
          params$s, rep(0,aux$J),
          aux$J,
          params$sigma2 )
        if( aux$shrinkage ) { # since shirnkage is requested, the prior for Z should also be updated
          # use an emprical Bayes approach
          # first, solve Z, but with no prior
          multipliers <- sqrt(1/params$sigma2)
          Z_noprior <- solveZ(
            eigenRowMult( target$Y, multipliers ),
            params$s * multipliers, rep(0,aux$J),
            eigenRowMult( params$B, multipliers ),
            aux$diag_K,
            matrix(0,nrow=aux$K,ncol=aux$K) )
          aux$Z_noprior <<- Z_noprior
          # Use the law of total covariance to calculate the covariance of Z
          Covar_E <- cov(Z_noprior)
          E_Covar <- solve(tcrossprod(eigenRowMult( params$B, multipliers )))
          hyperparams$Z_cov <<- (Covar_E+E_Covar)*aux$cov_mask
          hyperparams$Z_precision <<- solve(hyperparams$Z_cov)
          prmatrix(
            round(hyperparams$Z_cov,digits=2),
            rowlab=rep("",aux$K),collab=rep("",aux$K))
        }
      }

      convergence$sigma2 <<-
        rbind( convergence$sigma2, params$sigma2 )
    },
    
    ######
    # Solve Yi using raw counts Mi
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZBi and QiDBi to be up-to-date
    solve.Y = function( i )
    {
      if( aux$optimize_all ) {
        target$Y <<- solveY(
          target$M,
          target$Y,
          aux$ZB,
          params$s, rep(0,aux$J),
          params$sigma2 )
      } else {
        target$Y[aux$optimize_indices,] <<- solveY(
          target$M[aux$optimize_indices,,drop=F],
          target$Y[aux$optimize_indices,,drop=F],
          aux$ZB[aux$optimize_indices,,drop=F],
          params$s, rep(0,length(aux$optimize_indices)),
          params$sigma2 )
      }

      #return(1)
    },
    
    #######################
    # Function for calculating convergence metrics
    #######################
    
    # `ref`: the reference relative to which the parameter changes are calculated
    track_changes = function( ref.Z, ref.Y, ref.s ) {
      convergence$dZ <<- apply(abs(params$Z-ref.Z),1,max)
      convergence$dY <<- apply(abs(target$Y-ref.Y),1,max)
      convergence$ds <<- abs(params$s-ref.s)
    },
    
    # `return TRUE if the model has converged
    check_convergence = function( etol ) {
      # If the global parameters are fixed, they are assumed converged
      if( aux$fixed_global_params ) {
        max_ds <- 0
        max_dsigma2 <- 0
        sigma_converged <- T
        s_converged <- T
      } else { # otherwise, check their convergence explicitly
        if( any(is.na(convergence$ds)) ) { # ds is not calculated yet
          max_ds <- NA
          s_converged <- F
        } else {
          max_ds <- max(convergence$ds)
          if(  max_ds < etol ) {
            s_converged <- T
          } else {
            s_converged <- F
          }
        } 
        if( nrow(convergence$sigma2) <= 1 ) { # only one sigma2 calculation is available
          max_dsigma2 <- NA
          sigma_converged <- F
        } else {
          # dsigma2 is the log fold-change between the last two iterations
          # max_dsigma2 <-
          #   max( abs( log( convergence$sigma2[nrow(convergence$sigma2),] /
          #                    convergence$sigma2[nrow(convergence$sigma2)-1,] ) ) )
          
          # dsigma2 is the difference between standard deviations in the last two iterations
          max_dsigma2 <-
            max( abs( sqrt( convergence$sigma2[nrow(convergence$sigma2),] ) -
                        sqrt( convergence$sigma2[nrow(convergence$sigma2)-1,] ) ) )
          if( max_dsigma2 < etol ) {
            sigma_converged <- T
          } else {
            sigma_converged <- F
          }
        } 
      }
      
      max_dZ <- max( convergence$dZ )
      if( max_dZ < etol ) {
        Z_converged <- T
      } else {
        Z_converged <- F
      }

      max_dY <- max( convergence$dY )
      if( max_dY < etol ) {
        Y_converged <- T
      } else {
        Y_converged <- F
      }
      message("    dSigma2:",max_dsigma2," dS:",max_ds," dZ:",max_dZ," dY:",max_dY)
      
      if( !aux$optimize_all ) {
        aux$optimize_indices <<-
          which( convergence$dZ > etol | convergence$dY > etol )
        message("    N:",length(aux$optimize_indices))
      }
    
      return( sigma_converged & s_converged & Y_converged & Z_converged )
    },
    
    
    #######################
    # Function for iterative optimization of the model
    #######################
    
    # `iterations`: the number of iterations for this round of optimization
    # `track_interval`: the interval for convergence statistics to be calculated and stored
    optimize = function( iterations=1000, track_internval=5, etol=1e-4 ) {
      
      message("Performing block coordinate descent optimization...")
      
      # iteratively optimize the model through block coordinate descent
      for( ite in 1:iterations ) {
        aux$ite <<- aux$ite+1
        # if convergence info needs to be updated, store the current model parameters
        if( aux$ite %% track_internval == 1 | track_internval==1 ) {
          prev_Z <- params$Z
          prev_Y <- target$Y
          prev_s <- params$s
          tic( paste0(
            "  Iteration ",ite,"/",iterations," (total ",aux$ite,")") )
        }
        # Solve model parameters
        solve.Z()
        update.ZB() # update ZB
        solve.s()
        solve.Y()
        solve.sigma2()
        # Calculate and store the convergence statistics
        if( aux$ite %% track_internval == 1 | track_internval==1 ) {
          toc()
          track_changes( prev_Z, prev_Y, prev_s )
          if( check_convergence( etol ) ) { # use the convergence metrics to check for convergence
            break
          }
        }
      }
      if( !check_convergence( etol ) ) {
        message("Warning: model did not reach convergence criteria by the end of iterations.")
      }
    },
    
    
    
    
    #######################
    # Function for fitting a reduced model to data, followed by LRT
    #######################
    lrt = function(
      rmIndex, # The index of the variable to be removed from the design matrix
      topEntries=100000, # The LRT will be performed for these many top genes only
      coeffSign=c("positive","negative","both"),
      iterations=1000, track_internval=5, etol=1e-4 ) {
  
      ##### store the fitted full model
      aux$B.full <<- params$B
      aux$Z.full <<- params$Z
      aux$Y.full <<- target$Y
      aux$M.full <<- target$M
      
      # Create a template data frame for the results      
      res <- data.frame(
        name=aux$geneIDs,
        coefficient.br = params$Z[,rmIndex], # before refinement
        coefficient.ar = NA, # after refinement
        full_LL = NA,
        reduced_LL = NA,
        pvalue = NA,
        fdr = NA )
      
      
      # Identify the top entries that should be kept for lrt
      if( coeffSign[1]=="positive" ) {
        topIndices <- order( params$Z[,rmIndex], decreasing=T )[1:topEntries]
      } else if( coeffSign[1]=="negative" ) {
        topIndices <- order( params$Z[,rmIndex], decreasing=F )[1:topEntries]
      } else if( coeffSign[1]=="both" ) {
        topIndices <- order( abs(params$Z[,rmIndex]), decreasing=T )[1:topEntries]
      } else {
        stop("Unrecognized coeffSign.")
      }

      # Fix the global parameters, and optimize only the genes that have not converged
      aux$fixed_global_params <<- T
      aux$optimize_all <<- F
      aux$optimize_indices <<- 1:topEntries
      aux$shrinkage <<- F
      
      message("Refining the full model for top candidates, with no shrinkage...")
      ##### Keep only the top candidates and optimize until convergence
      params$Z <<- params$Z[topIndices,,drop=F]
      target$Y <<- target$Y[topIndices,,drop=F]
      target$M <<- target$M[topIndices,,drop=F]
      optimize(iterations=iterations,track_internval=track_internval,etol=etol)
      # Store the results
      params$Z.full.noShrinkage <<- params$Z
      target$Y.full.noShrinkage <<- target$Y
      res$coefficient.ar[topIndices] <- params$Z[,rmIndex]
      
      
      message("Performing block coordinate descent optimization for the reduced model...")
      ##### modify B and Z by removing the requested columns
      params$Z <<- params$Z[,-rmIndex,drop=F]
      params$B <<- params$B[-rmIndex,,drop=F]
      aux$K <<- nrow(params$B) # recalculate the dimension of the model
      aux$diag_K <<- diag(aux$K)
      # Repopulate the vector of genes that need to be optimized
      aux$optimize_indices <<- 1:topEntries
      optimize(iterations=iterations,track_internval=track_internval,etol=etol)
      
      ##### store the fitted reduced model, and then restore the full model
      params$B.reduced <<- params$B
      params$Z.reduced <<- params$Z
      target$Y.reduced <<- target$Y
      aux$reduced.indices <<- topIndices
      # restore the full model
      params$B <<- aux$B.full
      params$Z <<- aux$Z.full
      target$Y <<- aux$Y.full
      target$M <<- aux$M.full
      aux$fixed_global_params <<- F
      aux$optimize_all <<- T
      aux$optimize_indices <<- 1:aux$J
      aux$K <<- nrow(params$B)
      aux$diag_K <<- diag(aux$K)
      aux$shrinkage <<- T
      aux$B.full <<- NULL
      aux$Z.full <<- NULL
      aux$Y.full <<- NULL
      aux$M.full <<- NULL
      gc()
      
      # Yhat is the estimated fragment abundance based on the fitted model
      Yhat.full <-
        t( t( params$Z.full.noShrinkage %*% params$B )
           + params$s )
      Yhat.reduced <-
        t( t( params$Z.reduced %*% params$B.reduced ) +
             params$s )

      # Function to calculate log-likelihood for a genomic region
      lrt.local <- function(x,sigma2) {
        m <- x[1:aux$N]
        mu <- x[(aux$N+1):(2*aux$N)]
        sum( sapply(1:aux$N, function(x) { dpoilog(m[x],mu[x],sigma2[x],T) } ) )
      }
      
      message("Calculating the log-likelihoods for the full model...")
      full_LL <- apply(
        cbind(target$M[topIndices,,drop=F],Yhat.full), 1,
        lrt.local, sigma2=params$sigma2 )
      message("Calculating the log-likelihoods for the reduced model...")
      reduced_LL <- apply(
        cbind(target$M[topIndices,,drop=F],Yhat.reduced), 1,
        lrt.local, sigma2=params$sigma2 )
      res$full_LL[topIndices] <- full_LL
      res$reduced_LL[topIndices] <- reduced_LL
      
      res$pvalue[topIndices] <- pchisq(2*(full_LL-reduced_LL),df = 1,lower.tail = F)
      lrt_p <- res$pvalue
      lrt_p[is.na(lrt_p)] <- 1
      lrt_fdr <- p.adjust(lrt_p,method="fdr")
      res$fdr[topIndices] <- lrt_fdr[topIndices]
      
      return(res)
    }
  )
)

################ Auxiliary functions

# Loading libraries and scripts that are needed for auxiliary functions
#library(jointseg)

#######################
# Function for imputing Y given raw counts, after model fitting is complete
#
# The effect of sample-specific distortions and cell-specific library sizes
#######################

# `logScale`: whether the imputed values should be on the logarithmic scale
# `rowCentre`: whether the global gene-specific offsets should be removed from imputed values
getY.lnpReg <- function( object, logScale=T ) {
  # Remove the effects of sample-specific factors from Yi, and concatenate
  Y_res <- scale( object$target$Y, center=object$params$s, scale=F)
  rownames(Y_res) <- object$aux$geneIDs
  colnames(Y_res) <- object$aux$sampleIDs
  # Return the imputed values
  if( logScale ) {
    return( Y_res ) # Y_res is itself on the log-scale
  } else { # If observation is simply the counts, return exp of Y
    return( exp(Y_res) )
  }
}

#######################
# Function to return the lnpReg B, Z or ZB projection
#######################

getB.lnpReg <- function( object ) {
  return( object$params$B )
}

getZ.lnpReg <- function( object ) {
  Z <- object$params$Z
  # set the row and column names, and return the result
  rownames(Z) <- object$aux$geneIDs
  colnames(Z) <- rownames(object$params$B)
  return( Z )
}

getZB.lnpReg <- function( object ) {
  # Concatenate the ZxDxBi matrices
  ZB <- object$params$Z %*% object$params$B
  # set the row and column names, and return the result
  rownames(ZB) <- object$aux$geneIDs
  colnames(ZB) <- object$aux$sampleIDs
  return( ZB )
}
