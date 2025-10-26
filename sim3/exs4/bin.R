#############################################################################################
# This script corresponds to Example S4 in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

# Load required libraries
library(bcv)
library(ICtest)
library(GFM)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(POET)
library(openxlsx)
# we use "mirtjml" to deal with binary factor model
library(mirtjml)
# Set base directory for saving results
base_dir <- "~/glmcv/simuV2/result_bin/V2bin"

# Create directory if it doesn't exist
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# Load required functions
source("~/glmcv/simuV2/glmcv_functions.R")

# Simulation parameters
q <- 2
n_simulations <- 100
# Define (n,p) combinations to test
n_values <- c(50,100,300)
p_values <- c(50,100,300)



# Main loop for (n,p) combinations
for (n in n_values) {
  for (p in p_values) {
    
    # cat("\n", paste(rep("=", 50), collapse=""), "\n")
    # cat("Starting simulation with n =", n, "p =", p, "q =", q, "\n")
    # cat(paste(rep("=", 50), collapse=""), "\n")
    # 
    # Initialize results matrix for this (n,p) combination
    results <- matrix(0, nrow = n_simulations, ncol = 63)
    
    # Simulation loop
    for (s in 1:n_simulations) {
      # cat("Processing simulation:", s, "out of", n_simulations, "\n")
      # 
      # Set seed for reproducibility
      set.seed(s)

      re_all <- tryCatch({
        
        # Generate simulation data
        data <- sim1_bin(n, p, q)
        resp <- data$resp
        
        C_ada=max(estimate_C_binary(resp,qmax = 9)$C_est,15)
        # cat(" For ",s," repeat, C_ada = ",C_ada,"\n")
        # 
        
        # Initialize factor loading matrices
        theta0 <- cbind(data$F.matr, matrix(runif(n * (8 - q), -0.2, 0.2), n, (8 - q)))
        A0 <- cbind(data$A.matr, matrix(runif(p * (8 - q), -0.2, 0.2), p, (8 - q)))
        d0 = data$d.vec
        # Initialize result matrix for cross-validation methods
        mat_es <- matrix(0, nrow = 5, ncol = 10)
        
        # ===============================
        # Cross-Validation Methods
        # ===============================
     
        for (K in c(5)) {
          split_matrices <- split_matrix(n, p, K)
          
          # Initialize loss storage
          loss_values <- matrix(0, nrow = K, ncol = 8)
          
          # Perform K-fold entrywise splitting cross-validation
          for (k in 1:K) {
            # Create train/test masks
            test_mask <- split_matrices[[k]]
            train_mask <- Reduce("+", split_matrices[-k]) > 0
            
            # Create training set
            train_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
            train_set[train_mask] <- resp[train_mask]
            test_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
            test_set[test_mask] <- resp[test_mask]
            num_train <- sum(train_mask == 1)
            
            # Fit models with different numbers of factors
            for (q_est in 1:8) {
              tryCatch({
                jml_res <- mirtjml_expr(
                  train_set, q_est, 
                  theta0 = theta0[, 1:(q_est)], 
                  A0 = A0[, 1:(q_est)], d0 = d0,
                  tol = 0.01*num_train, cc = C_ada,
                  print_proc = F
                )
                #
                jml_result <- JL_bin(test_set, jml_res$A_hat, jml_res$theta_hat, jml_res$d_hat)
                loss_values[k, q_est] <- -jml_result$lik
              }, error = function(e) {
                loss_values[k, q_est] <- NA  
              })
            }
          }
          
          # Determine row index for results
          r <- 2
          
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
          cv_loss <- apply(loss_values, 2, function(x) {
            non_na <- x[!is.na(x)]
            if (length(non_na) == 0) Inf else mean(non_na)
          })

          mat_es[r, 6] <- which.min(cv_loss)
          mat_es[r, 7] <- which.min(cv_loss + (1 / (K - 1)) * (1:8) * pen1)
          mat_es[r, 8] <- which.min(cv_loss + (1 / (K - 1)) * (1:8) * pen2)
          mat_es[r, 9] <- which.min(cv_loss + (( K) / (K - 1)) * (1:8) * pen1)
          mat_es[r, 10] <- which.min(cv_loss + (( K) / (K - 1)) * (1:8) * pen2)
        }


        
        # ===============================
        # Information Criterion Methods
        # ===============================
        
        # Fit full models for information criteria
        true_matrix <- matrix(TRUE, nrow = n, ncol = p)

        # Fit models with different numbers of factors
        jic_results <- list()
        lnlik <- numeric(8)
        
        
        for (q_est in 1:8) {
          jic_results[[q_est]] <- mirtjml_expr(
            resp, q_est, 
            theta0[, 1:(q_est)], A0[, 1:(q_est)], 
            d0=d0,tol = 0.01*n*p,
            cc = C_ada, print_proc = F
          )
          
          jic_result <- JL_bin(resp, jic_results[[q_est]]$A_hat, jic_results[[q_est]]$theta_hat, jic_results[[q_est]]$d_hat)
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
        
        # cat("IC success \n")
        # ===============================
        # Other Methods
        # ===============================

        # Bi-cross validation methods
        bcv1 <- cv.svd.gabriel(resp, krow = 2, kcol = 2, maxrank = 8)
        re_bcv1 <- which.min(as.vector(colMeans(bcv1$msep))[-1])
        
        #bcv2 <- cv.svd.wold(resp, k = 5, maxrank = 8)
        re_bcv2 <- 0
        #which.min(as.vector(colMeans(bcv2$msep))[-1])
        
        # cat(" CV success \n")
        # 
        # Ladle method
        #ladle <- PCAladle(resp, ncomp = 8)
        re_ladle <- 0
        #ladle$k
        re_ratio <- 0

        # Sample CV methods
        samplecv2_result <- factor_selection_cv_binary(
          resp = resp,
          theta0 = theta0,
          A0 = A0,
          d0 = d0,
          K = 5,
          max_factors = 8,
          C = C_ada,
          tol_factor = 0.01,
          see = s
        )
        re_samplecv2 <- samplecv2_result$optimal_factors
        
        # cat("sample cv2 success \n")
        cv_results_bayes <- factor_selection_cv_bayes(
          resp, 
          qmax = 8, 
          K_folds = 5, 
          seed = s
        )
        re_samplecv_bayes <- cv_results_bayes$facnum_samplecv
        
        # cat("bayes CV success \n")
        # # Print results for this simulation
        # cat("Simulation", s, "completed\n")
        # cat("JIC1:", re_jic1, "JIC2:", re_jic2, "\n")
        # cat("AIC1:", re_aic1, "AIC2:", re_aic2, "\n")
        # cat("BIC1:", re_bic1, "BIC2:", re_bic2, "BIC3:", re_bic3, "\n")
        # cat("BCV1:", re_bcv1, "BCV2:", re_bcv2, "\n")
        # cat("Ladle:", re_ladle, "Ratio:", re_ratio, "\n")
        # cat("sample CV:", re_samplecv2, "sample CV bayes:", re_samplecv_bayes, "\n")
        # # Return results vector
        c(
          as.vector(t(mat_es)),  # CV results (50 values)
          re_jic1, re_jic2,     # JIC results (2 values)
          re_aic1, re_aic2,     # AIC results (2 values)
          re_bic1, re_bic2, re_bic3,  # BIC results (3 values)
          re_bcv1, re_bcv2,     # BCV results (2 values)
          re_ladle,            
          re_ratio,            
          re_samplecv2,         # Sample CV2 result (1 value)
          re_samplecv_bayes     # Sample CV Bayes result (1 value)
        )
        
      }, error = function(e) {
        # cat("Error in simulation", s, ":", e$message, "\n")
        # cat("Error in simulation", s, ":", e$message, "\n", 
        #     file = file.path(base_dir, "error_log.txt"), append = TRUE)
        # # Return zero vector if error occurs
        rep(0, 63)
      })
      
      # Store results
      results[s, ] <- re_all
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
    
    # Save results for this (n,p) combination
    filename_facnum <- sprintf("%s/n%d_p%d_q%d_facnum-bin-cada-u8q2-tol001.rds", base_dir, n, p, q)
    filename_frequency <- sprintf("%s/n%d_p%d_q%d_frequency-bin-cada-u8q2-tol001.rds", base_dir, n, p, q)
    
    saveRDS(results, file = filename_facnum)
    saveRDS(re_frequ, file = filename_frequency)
    
    
    
    # Save results for this (n,p) combination
    filename_facnum <- sprintf("%s/n%d_p%d_q%d_facnum-bin-cada-u8q2-tol001.xlsx", base_dir, n, p, q)
    filename_frequency <- sprintf("%s/n%d_p%d_q%d_frequency-bin-cada-u8q2-tol001.xlsx", base_dir, n, p, q)
    
    # excel
    write.xlsx(results, file = filename_facnum)
    write.xlsx(re_frequ, file = filename_frequency)
    
    cat("\nResults for n =", n, ", p =", p, "saved to:\n")
    cat("Factor numbers:", filename_facnum, "\n")
    cat("Frequencies:", filename_frequency, "\n")
    cat("Completed (n,p) = (", n, ",", p, ")\n")
  }
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("All simulations completed successfully!\n")
cat("Results saved for", length(n_values) * length(p_values), "(n,p) combinations\n")
cat(paste(rep("=", 60), collapse=""), "\n")




