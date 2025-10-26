#############################################################################################
# This script corresponds to real data 10 fold ECV.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

# The mouse brain real dataset was downloaded from:
# https://www.kaggle.com/datasets/aayush9753/singlecell-rnaseq-data-from-mouse-brain

# loading packages
library(readr)
library(POET)
library(irlba)
library(bcv)
library(ICtest)
library(GFM)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)


brain_metadata <- read_csv("~/glmcv/realV2/gene/brain_metadata.csv")
brain_counts <- read_csv("~/glmcv/realV2/gene/brain_counts.csv")
colnames(brain_metadata)[1] <- "ID"
colnames(brain_counts)[1] <- "ID"

merged_data <- merge(brain_metadata, brain_counts, by = "ID", all = TRUE)
filtered_data <- merged_data[merged_data$mouse.sex == "F", ]
brain_count_f <- filtered_data[, -c(1:6)]

if (any(is.na(brain_count_f))) {
  #cat("move NA...\n")
  brain_count_f <- brain_count_f[, colSums(is.na(brain_count_f)) == 0]  
}

if (any(colSums(brain_count_f) == 0)) {
  #cat("move zero gene...\n")
  brain_count_f <- brain_count_f[, colSums(brain_count_f) > 0] 
}

if (any(rowSums(brain_count_f) == 0)) {
  cat("move zero cells...\n")
  brain_count_f <- brain_count_f[rowSums(brain_count_f) > 0, ]  
}
# print("dim(brain_count_f)")
# print(dim(brain_count_f))
# 707 18439


library(Seurat)
################### data processing
seurat_obj <- CreateSeuratObject(counts = t(brain_count_f))

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 1000)

highly_variable_genes <- VariableFeatures(seurat_obj)

brain_count_f1000 <- brain_count_f[, highly_variable_genes]
quantile_limit <- quantile(brain_count_f1000, 0.97, na.rm = TRUE)
# print("quantile_limit")
# print(quantile_limit)
brain_count_f1000[brain_count_f1000 > quantile_limit] <- quantile_limit
print(dim(brain_count_f1000))




for (n in c(707)) {
  for (p in c(1000)) {
    # 
    # cat("\n", paste(rep("=", 50), collapse=""), "\n")
    # cat("Real data start sample =", n, "feature =", p, "\n")
    # cat(paste(rep("=", 50), collapse=""), "\n")
    # 
    # Initialize results matrix for this (n,p) combination
    results <- matrix(0, nrow = 1, ncol = 63)
    
    Rcpp::sourceCpp("~/glmcv/simuV2/confirm_jmle_omp_poisson_with_intercept_missing.cpp")
    source("~/glmcv/simuV2/glmcv_functions.R")
    
    # Set seed for reproducibility
    set.seed(666)
    
    resp = as.matrix(brain_count_f1000)
    est_c=estimate_C(log(resp+1),qmax = 9)

    init=GFM::gfm(list(resp),types = 'poisson',q=8,verbose = F)
    theta0=cbind(rep(1,n),init$hH)
    A0 <- cbind(init$hmu,init$hB)
    
    
    C_ada=(est_c$C_est)
    # cat("  C_ada = ",C_ada,"\n")
    
    # Initialize result matrix for cross-validation methods
    mat_es <- matrix(0, nrow = 5, ncol = 10)
    
    # ===============================
    # Cross-Validation Methods
    # ===============================
    
    # ECV with K = 10 only
    for (K in c(10)) { 
      set.seed(K+2024)
      split_matrices <- split_matrix(n, p, K)
      
      # Initialize loss storage
      loss_values <- matrix(0, nrow = K, ncol = 8)
      
      # Perform K-fold cross-validation
      for (k in 1:K) {
        #cat("Then start ", k, " fold  \n")
        # Create train/test masks
        test_mask <- split_matrices[[k]]
        train_mask <- Reduce("+", split_matrices[-k]) > 0
        
        # Create training set
        train_set <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
        train_set[train_mask] <- resp[train_mask]
        
        num_train <- sum(train_mask == 1)
        tol_val <- 0.01 / num_train
        
        # Fit models with different numbers of factors
        for (q_est in 1:8) {
          jml_res <- confirm_CJMLE_poisson_cpp(
            train_set, train_mask, 
            theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)], 
            matrix(TRUE, p, q_est + 1), 
            C = C_ada, tol = tol_val
          )
          
          jml_result <- JL(resp, jml_res$theta, jml_res$A, test_mask)
          loss_values[k, q_est] <- -jml_result$lik
          cat(" q = ",q_est, "\n")
        }
        
        cat_es_loss <- loss_values[k, ]  # ES (first fold)
        pen1 <- max(n, p) * log(min(n, p))
        pen2 <- (n + p) * log((n * p) / (n + p))
        
        # ES methods
        # 
        # cat("cat_es_loss = ",(cat_es_loss),"\n")
        # cat("cat_es_loss+p1 = ",(cat_es_loss + (1 / (K - 1)) * (1:8) * pen1),"\n")
        # cat("cat_es_loss+p2 = ",(cat_es_loss + (1 / (K - 1)) * (1:8) * pen2),"\n")
        # cat("cat_es_loss+p3 = ",(cat_es_loss + (K / (K - 1)) * (1:8) * pen1),"\n")
        # cat("cat_es_loss+p4 = ",(cat_es_loss + (K / (K - 1)) * (1:8) * pen2),"\n")
        # 
        # 
        # cat("cat_es_loss = ",which.min(cat_es_loss),"\n")
        # cat("cat_es_loss+p1 = ",which.min(cat_es_loss + (1 / (K - 1)) * (1:8) * pen1),"\n")
        # cat("cat_es_loss+p2 = ",which.min(cat_es_loss + (1 / (K - 1)) * (1:8) * pen2),"\n")
        # cat("cat_es_loss+p3 = ",which.min(cat_es_loss + (K / (K - 1)) * (1:8) * pen1),"\n")
        # cat("cat_es_loss+p4 = ",which.min(cat_es_loss + (K / (K - 1)) * (1:8) * pen2),"\n")
        # 
        
      }
      
      # Determine row index for results
      # if (K == 2) r <- 1  # COMMENTED OUT
      # else
      if (K == 5) r <- 2  # COMMENTED OUT
      else if (K == 10)  r <- 3  # COMMENTED OUT
      # r <- 3 # MODIFIED: Directly set r = 2 for K = 5
      
      # Calculate different penalty methods
      es_loss <- loss_values[1, ]  # ES (first fold)
      pen1 <- max(n, p) * log(min(n, p))
      pen2 <- (n + p) * log((n * p) / (n + p))
      
      # ES methods
      mat_es[r, 1] <- which.min(es_loss)
      mat_es[r, 2] <- which.min(es_loss + (1 / (K - 1)) * (1:8) * pen1)
      mat_es[r, 3] <- which.min(es_loss + (1 / (K - 1)) * (1:8) * pen2)
      mat_es[r, 4] <- which.min(es_loss + (K / (K - 1)) * (1:8) * pen1)
      mat_es[r, 5] <- which.min(es_loss + (K / (K - 1)) * (1:8) * pen2)
      
      # ECV methods
      cv_loss <- colSums(loss_values)
      
      # cat("cv_loss = ",(cv_loss/(n*p)),"\n")
      # cat("cv_loss+p1 = ",((cv_loss + (K / (K - 1)) * (1:8) * pen1)/(n*p)),"\n")
      # cat("cv_loss+p2 = ",((cv_loss + (K / (K - 1)) * (1:8) * pen2)/(n*p)),"\n")
      # cat("cv_loss+p3 = ",((cv_loss + ((K * K) / (K - 1)) * (1:8) * pen1)/(n*p)),"\n")
      # cat("cv_loss+p4 = ",((cv_loss + ((K * K) / (K - 1)) * (1:8) * pen2)/(n*p)),"\n")
      # 
      mat_es[r, 6] <- which.min(cv_loss)
      mat_es[r, 7] <- which.min(cv_loss + (K / (K - 1)) * (1:8) * pen1)
      mat_es[r, 8] <- which.min(cv_loss + (K / (K - 1)) * (1:8) * pen2)
      mat_es[r, 9] <- which.min(cv_loss + ((K * K) / (K - 1)) * (1:8) * pen1)
      mat_es[r, 10] <- which.min(cv_loss + ((K * K) / (K - 1)) * (1:8) * pen2)
    }
    
    
    print("Cross-validation results:")
    print(mat_es)
    
    # ===============================
    # Information Criterion Methods
    # ===============================
    
    # Fit full models for information criteria
    true_matrix <- matrix(TRUE, nrow = n, ncol = p)
    tol_full <- 0.01 / (n * p)
    
    # Fit models with different numbers of factors
    jic_results <- list()
    lnlik <- numeric(8)
    
    for (q_est in 1:8) {
      jic_results[[q_est]] <- confirm_CJMLE_poisson_cpp(
        resp, true_matrix, 
        theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)], 
        matrix(TRUE, p, q_est + 1), 
        C = C_ada, tol = tol_full
      )
      
      jic_result <- JL(resp, jic_results[[q_est]]$theta, jic_results[[q_est]]$A, true_matrix)
      lnlik[q_est] <- jic_result$lik
    }
    
    # Calculate penalties
    pen1 <- max(n, p) * log(min(n, p))
    pen2 <- (n + p) * log((n * p) / (n + p))
    
    # JIC methods
    jic1 <- -2 * lnlik + (1:8) * pen1
    jic2 <- -2 * lnlik + (1:8) * pen2
    re_jic1 <- which.min(jic1)
    re_jic2 <- which.min(jic2)
    
    cat("jic1 = ",re_jic1,"\n")
    cat("jic2 = ",re_jic2,"\n")
    # AIC methods
    aic1 <- -2 * lnlik + (1:8) * (2 * p)
    aic2 <- -2 * lnlik + (1:8) * (2 * n)
    re_aic1 <- which.min(aic1)
    re_aic2 <- which.min(aic2)
    
    # BIC methods
    bic1 <- -2 * lnlik + (1:8) * (log(n) * p)
    bic2 <- -2 * lnlik + (1:8) * (log(p) * n)
    bic3 <- numeric(8)
    for (qq in 1:8) {
      bic3[qq] <- -2 * lnlik[qq] + qq * (n + p - qq) * log(n * p)
    }
    re_bic1 <- which.min(bic1)
    re_bic2 <- which.min(bic2)
    re_bic3 <- which.min(bic3)
    
    cat("bic3 = ",re_bic3,"\n")
    # ===============================
    # Other Methods
    # ===============================
    
    # Log-transformed response for some methods
    logresp <- log(resp + 1)
    
    # Bi-cross validation methods
    bcv1 <- cv.svd.gabriel(logresp, krow = 2, kcol = 2, maxrank = 8)
    re_bcv1 <- which.min(as.vector(colMeans(bcv1$msep))[-1])
    
    # bcv2 <- cv.svd.wold(logresp, k = 5, maxrank = 8)
    re_bcv2 <- 0
    #which.min(as.vector(colMeans(bcv2$msep))[-1])
    
    # Ladle method
    #ladle <- PCAladle(logresp, ncomp = 8)
    re_ladle <- NA
    #ladle$k

    re_ratio <- 0

    # Sample CV methods
    samplecv2_result <- factor_selection_cv(
      resp = resp,
      theta0 = theta0,
      A0 = A0,
      K = 5,
      max_factors = 8,
      C = C_ada,
      # also use
      tol_factor = 0.01,
      see = 1
    )
    re_samplecv2 <- samplecv2_result$optimal_factors
    
    cv_results_bayes <- factor_selection_cv_bayes(
      logresp, 
      qmax = 8, 
      K_folds = 5, 
      seed = 1
    )
    re_samplecv_bayes <- cv_results_bayes$facnum_samplecv
    
    # # Print results for this simulation
    # cat("Real completed\n")
    # cat("JIC1:", re_jic1, "JIC2:", re_jic2, "\n")
    # cat("AIC1:", re_aic1, "AIC2:", re_aic2, "\n")
    # cat("BIC1:", re_bic1, "BIC2:", re_bic2, "BIC3:", re_bic3, "\n")
    # cat("BCV1:", re_bcv1, "BCV2:", re_bcv2, "\n")
    # cat("Ladle:", re_ladle, "Ratio:", re_ratio, "\n")
    # cat("sample CV:", re_samplecv2, "sample CV bayes:", re_samplecv_bayes, "\n")
    # # Return results vector
    re_all= c(
      as.vector(t(mat_es)),  # CV results (50 values)
      re_jic1, re_jic2,     # JIC results (2 values)
      re_aic1, re_aic2,     # AIC results (2 values)
      re_bic1, re_bic2, re_bic3,  # BIC results (3 values)
      re_bcv1, re_bcv2,     # BCV results (1 values)
      re_ladle,             
      re_ratio,             
      re_samplecv2,         # Sample CV2 result (1 value)
      re_samplecv_bayes     # Sample CV Bayes result (1 value)
    )
    
    # Store results
    results[1, ] <- re_all
  }
  
  # ===============================
  # Process and Save Results for this (n,p) combination
  # ===============================
  
  # Define column names corresponding to estimation methods
  column_names <- c(
    # Cross-validation methods (K=2)
    "ES_K2", "ES_K2_pen1", "ES_K2_pen2", "ES_K2_pen3", "ES_K2_pen4",
    "GCV_K2", "GCV_K2_pen1", "GCV_K2_pen2", "GCV_K2_pen3", "GCV_K2_pen4",
    
    # Cross-validation methods (K=5) 
    "ES_K5", "ES_K5_pen1", "ES_K5_pen2", "ES_K5_pen3", "ES_K5_pen4",
    "GCV_K5", "GCV_K5_pen1", "GCV_K5_pen2", "GCV_K5_pen3", "GCV_K5_pen4",
    
    # Cross-validation methods (K=10)
    "ES_K10", "ES_K10_pen1", "ES_K10_pen2", "ES_K10_pen3", "ES_K10_pen4",
    "GCV_K10", "GCV_K10_pen1", "GCV_K10_pen2", "GCV_K10_pen3", "GCV_K10_pen4",
    
    # Reverse cross-validation methods (K=5)
    "ES_rev_K5", "ES_rev_K5_pen1", "ES_rev_K5_pen2", "ES_rev_K5_pen3", "ES_rev_K5_pen4",
    "GCV_rev_K5", "GCV_rev_K5_pen1", "GCV_rev_K5_pen2", "GCV_rev_K5_pen3", "GCV_rev_K5_pen4",
    
    # Reverse cross-validation methods (K=10)
    "ES_rev_K10", "ES_rev_K10_pen1", "ES_rev_K10_pen2", "ES_rev_K10_pen3", "ES_rev_K10_pen4",
    "GCV_rev_K10", "GCV_rev_K10_pen1", "GCV_rev_K10_pen2", "GCV_rev_K10_pen3", "GCV_rev_K10_pen4",
    
    # Information criteria
    "JIC1", "JIC2",
    "AIC1", "AIC2",
    "BIC1", "BIC2", "BIC3",
    
    # Other methods
    "BCV_Gabriel", "BCV_Wold",
    "Ladle", "Ratio",
    "SampleCV2", "SampleCV_Bayes"
  )
  
  # Assign column names
  colnames(results) <- column_names
  
  # Calculate frequency matrix
  re_frequ <- sapply(1:ncol(results), function(j) {
    tab <- table(factor(results[, j], levels = 1:8))
    as.numeric(tab)
  })
  colnames(re_frequ) <- column_names
  rownames(re_frequ) <- paste("Factor", 1:8)
  
  # Print results for this (n,p) combination
  cat("\nResults Summary for n =", n, ", p =", p, ":\n")
  print(results)
  cat("\nFrequency Table for n =", n, ", p =", p, ":\n")
  print(re_frequ)
}



