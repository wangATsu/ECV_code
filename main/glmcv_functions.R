
#####################################################################
########## This script contains the functions used in the analysis.
#####################################################################


# load library
library(Rcpp)
library(mvtnorm)
library(irlba)




### function for Poisson factor model

sim1 <- function(N, J, K){
  F.matr <- matrix(runif(N*K,-2,2), N, K);
  F.matr <- cbind(rep(1,N), F.matr)
  A.matr <- matrix(runif(J*(K+1),-2,2), J, K+1);                 
  
  temp = F.matr %*% t(A.matr);
  resp = matrix(0, N, J);
  resp[] = rpois(N*J, lambda = exp(temp))
  
  list(resp = resp, F.matr = F.matr, A.matr = A.matr,  M = temp)
}


split_matrix <- function(n, p, K) {
  original_matrix <- matrix(1, nrow = n, ncol = p)
  split_matrices <- replicate(K, matrix(0, nrow = n, ncol = p), simplify = FALSE)
  for (i in 1:n) {
    for (j in 1:p) {
      counts <- rmultinom(1, size = original_matrix[i, j], prob = rep(1/K, K))
      for (k in 1:K) {
        split_matrices[[k]][i, j] <- counts[k]
      }
    }
  }
  split_matrices <- lapply(split_matrices, function(x) x == 1)
  
  return(split_matrices)
}


JL <- function(resp, theta, A, nonmis_ind){
  temp = theta%*% t(A)
  temp1 = resp * temp - exp(temp)  
  lik = sum(temp1[nonmis_ind])
  list(lik = lik, M = temp)
}


### function for binary factor model




JL_bin <- function(data, A, Theta, d){
  N = nrow(data)
  temp = Theta %*% t(A) + rep(1, N) %*% t(d)
  M = temp
  prob = 1/(1+exp(-temp))
  temp = data * log(prob) + (1-data) * log(1-prob);
  temp[is.na(temp)] = 0;
  list(lik = sum(temp), M = M)
  
}



### function for gaussian factor model

sim1_norm <- function(N, J, K){
  F.matr <- matrix(runif(N*K,-4,4), N, K);
  F.matr <- cbind(rep(1,N), F.matr)
  A.matr <- matrix(runif(J*(K+1),-4,4), J, K+1);                 
  
  temp = F.matr %*% t(A.matr);
  resp = matrix(0, N, J);
  resp[] = temp + rnorm(N*J)
  
  list(resp = resp, F.matr = F.matr, A.matr = A.matr,  M = temp)
}



JL_norm <- function(resp, theta, A, nonmis_ind){
  temp = theta%*% t(A)
  
  lik =sum(((log(dnorm(resp, mean = temp, sd = 1)))[nonmis_ind]))
  
  list(lik = lik, M = temp)
}




################## function for sample cv

factor_selection_cv <- function(resp, theta0, A0, K = 5, max_factors = 8, 
                                C = 5, tol_factor = 0.01,see=123) {

  if (ncol(theta0) < max_factors + 1) {
    stop("theta0  at least max_factors + 1")
  }
  if (ncol(A0) < max_factors + 1) {
    stop("A0 at least max_factors + 1")
  }
  
  n <- nrow(resp)
  p <- ncol(resp)
  
  tol <- tol_factor / (n * p)
  
  set.seed(see)  
  fold_indices <- sample(rep(1:K, length.out = n))
  
  cv_losses <- matrix(0, nrow = K, ncol = max_factors)
  
  for (k in 1:K) {
    #cat("Processing fold", k, "of", K, "\n")
    test_idx <- which(fold_indices == k)
    train_idx <- which(fold_indices != k)
    
    X_train <- resp[train_idx, , drop = FALSE]
    X_test <- resp[test_idx, , drop = FALSE]
    theta0_train <- theta0[train_idx, , drop = FALSE]

    true_matrix <- matrix(TRUE, nrow = nrow(X_train), ncol = ncol(X_train))
    
    n_test <- nrow(X_test)
    
    for (q in 1:max_factors) {
      #cat("  Testing q =", q, "\n")
      
      jic_res <- confirm_CJMLE_poisson_cpp(
        X_train, 
        true_matrix, 
        theta0_train[, 1:(q+1), drop = FALSE], 
        A0[, 1:(q+1), drop = FALSE], 
        matrix(TRUE, p, q+1), 
        C = C, 
        tol = tol
      )
      
      A_estimated <- jic_res$A
      alpha <- A_estimated[, 1]  
      B <- A_estimated[, 2:(q+1), drop = FALSE]  

      F_test <- matrix(0, nrow = n_test, ncol = q)

      for (i in 1:n_test) {
        y_i <- X_test[i, ]

        glm_data <- data.frame(
          y = y_i,
          alpha = alpha
        )
        
        for (j in 1:q) {
          glm_data[[paste0("B", j)]] <- B[, j]
        }
        

        formula_str <- paste("y ~ offset(alpha) +", 
                             paste(paste0("B", 1:q), collapse = " + "), 
                             "- 1")
        
        tryCatch({
          glm_result <- glm(as.formula(formula_str), 
                            data = glm_data, 
                            family = poisson(link = "log"))

          coefficients <- coef(glm_result)

          if (any(is.na(coefficients)) || any(is.nan(coefficients)) || any(is.infinite(coefficients))) {
            F_test[i, ] <- rep(0, q)
            #warning(paste("GLM coefficients contain NA/NaN/Inf for sample", i, ": using zero vector"))
          } else {
            F_test[i, ] <- coefficients
          }
          
        }, error = function(e) {
          F_test[i, ] <- rep(0, q)
          #warning(paste("GLM failed for sample", i, ": using zero vector"))
        })
      }
      
      
      eta <- F_test %*% t(B) + matrix(1, nrow = n_test, ncol = 1) %*% t(alpha)
      
      likelihood_loss <- -sum(X_test * eta - exp(eta))
      
      cv_losses[k, q] <- likelihood_loss
    }
  }
  
  mean_cv_losses <- colMeans(cv_losses)
  
  optimal_q <- which.min(mean_cv_losses)

  result <- list(
    optimal_factors = optimal_q,
    cv_losses = cv_losses,
    mean_cv_losses = mean_cv_losses,
    fold_indices = fold_indices
  )
  
  #cat("Optimal number of factors:", optimal_q, "\n")
  
  return(result)
}

factor_selection_cv_normal <- function(resp, theta0, A0, K = 5, max_factors = 8, 
                                       C = 5, tol_factor = 0.01, see = 123) {

  if (ncol(theta0) < max_factors + 1) {
    stop("theta0 at least max_factors + 1")
  }
  if (ncol(A0) < max_factors + 1) {
    stop("A0的列数必须至少为max_factors + 1")
  }
  
  n <- nrow(resp)
  p <- ncol(resp)
  
  tol <- tol_factor / (n * p)
  
  set.seed(see)  
  fold_indices <- sample(rep(1:K, length.out = n))
 
  cv_losses <- matrix(0, nrow = K, ncol = max_factors)
  
  for (k in 1:K) {
    #cat("Processing fold", k, "of", K, "\n")
   
    test_idx <- which(fold_indices == k)
    train_idx <- which(fold_indices != k)
    
    X_train <- resp[train_idx, , drop = FALSE]
    X_test <- resp[test_idx, , drop = FALSE]
    
    theta0_train <- theta0[train_idx, , drop = FALSE]
    
    true_matrix <- matrix(TRUE, nrow = nrow(X_train), ncol = ncol(X_train))
    
    n_test <- nrow(X_test)
    
    for (q in 1:max_factors) {

      jic_res <- CJMLE_linear(
        X_train, 
        true_matrix, 
        theta0_train[, 1:(q+1), drop = FALSE], 
        A0[, 1:(q+1), drop = FALSE], 
        matrix(TRUE, p, q+1), 1,
        C = C, 
        tol = tol,F
      )
      
      A_estimated <- jic_res$A
      alpha <- A_estimated[, 1] 
      B <- A_estimated[, 2:(q+1), drop = FALSE]
     
      F_test <- matrix(0, nrow = n_test, ncol = q)
      for (i in 1:n_test) {
        y_i <- X_test[i, ]
        
        glm_data <- data.frame(
          y = y_i,
          alpha = alpha
        )
       
        for (j in 1:q) {
          glm_data[[paste0("B", j)]] <- B[, j]
        }
       
        formula_str <- paste("y ~ offset(alpha) +", 
                             paste(paste0("B", 1:q), collapse = " + "), 
                             "- 1")
        tryCatch({
          glm_result <- glm(as.formula(formula_str), 
                            data = glm_data, 
                            family = gaussian(link = "identity"))
         
          coefficients <- coef(glm_result)
          if (any(is.na(coefficients)) || any(is.nan(coefficients)) || any(is.infinite(coefficients))) {
            F_test[i, ] <- rep(0, q)
           # warning(paste("GLM coefficients contain NA/NaN/Inf for sample", i, ": using zero vector"))
          } else {
            F_test[i, ] <- coefficients
          }
          
        }, error = function(e) {
          F_test[i, ] <- rep(0, q)
          #warning(paste("GLM failed for sample", i, ": using zero vector"))
        })
      }
      mu <- F_test %*% t(B) + matrix(1, nrow = n_test, ncol = 1) %*% t(alpha)
      
      likelihood_loss <- sum((X_test - mu)^2) / 2
      
      cv_losses[k, q] <- likelihood_loss
    }
  }
  mean_cv_losses <- colMeans(cv_losses)
  
  optimal_q <- which.min(mean_cv_losses)
  
  result <- list(
    optimal_factors = optimal_q,
    cv_losses = cv_losses,
    mean_cv_losses = mean_cv_losses,
    fold_indices = fold_indices
  )
  
  #cat("Optimal number of factors:", optimal_q, "\n")
  
  return(result)
}



factor_selection_cv_bayes <- function(resp, qmax = 8, K_folds = 5, C=0.5,seed = 123) {

  set.seed(seed)
  

  n <- nrow(resp)
  p <- ncol(resp)
  resp_centered=resp

  fold_indices <- sample(rep(1:K_folds, length.out = n))
  

  loss_matrix <- matrix(0, nrow = K_folds, ncol = qmax)
  

 
  for (fold in 1:K_folds) {
    

    test_indices <- which(fold_indices == fold)
    train_indices <- which(fold_indices != fold)
    
    train_data <- resp_centered[train_indices, ]
    test_data <- resp_centered[test_indices, ]

    for (q in 1:qmax) {
      
      tryCatch({

        train_data_center=(apply(train_data, 2, function(x) x - mean(x)))
        re_poet <- POET(t(train_data_center), K = q)


        SigmaY_est <- re_poet$SigmaY
        

        if (nrow(SigmaY_est) != p || ncol(SigmaY_est) != p) {
          stop(sprintf("covariance dim: expectation %d×%d, practice %d×%d", 
                       p, p, nrow(SigmaY_est), ncol(SigmaY_est)))
        }
        
 
        eigen_values <- eigen(SigmaY_est, only.values = TRUE)$values
        if (any(eigen_values <= 1e-10)) {
 
         # cat(sprintf(" warning not positive\n"))
          SigmaY_est <- SigmaY_est + diag(1e-6, p)
        }
        
        n_test <- nrow(test_data)
        
        test_data_center=(apply(test_data, 2, function(x) x - mean(x)))
        log_likelihood <- 0
        for (i in 1:n_test) {
          test_sample <- test_data_center[i, , drop = FALSE]  
          log_lik_i <- dmvnorm(test_sample, 
                               mean = rep(0, p), 
                               sigma = SigmaY_est, 
                               log = TRUE)
          log_likelihood <- log_likelihood + log_lik_i
        }

        loss_matrix[fold, q] <- -log_likelihood / n_test  
        
      }, error = function(e) {

       # cat(sprintf("the %d fold, factor number %d estimates error: %s\n", fold, q, e$message))
        loss_matrix[fold, q] <- Inf
      })
    }
  }

  average_loss <- colMeans(loss_matrix, na.rm = TRUE)
  

  facnum_samplecv <- which.min(average_loss)
  
  #cat("\n=== CV results ===\n")
  #cat(sprintf("estimate number %d\n", facnum_samplecv))
  #cat(sprintf("loss: %.6f\n", average_loss[facnum_samplecv]))
  
  results <- list(
    facnum_samplecv = facnum_samplecv,
    average_loss = average_loss,
    loss_matrix = loss_matrix,
    fold_indices = fold_indices
  )
  
  return(results)
}




### code for choosing constraint constant C

estimate_C <- function(X, qmax = 8, safety = 1.2) {


  svd_res <- irlba::irlba(X, nv = qmax, nu = qmax)
  U <- svd_res$u           # n×qmax
  V <- svd_res$v           # p×qmax
  d <- svd_res$d[1:qmax]   

  D_sqrt <- diag(sqrt(d), nrow = qmax, ncol = qmax)
  A_hat  <- U %*% D_sqrt   # n×qmax
  B_hat  <- V %*% D_sqrt   # p×qmax

  a_norms <- sqrt(rowSums(A_hat^2))  
  b_norms <- sqrt(rowSums(B_hat^2))  

  C_norm_hat <- max(c(a_norms, b_norms))
  C_est      <- C_norm_hat * safety
  
  list(
    qmax        = qmax,
    safety      = safety,
    C_norm_hat  = C_norm_hat,  
    C_est       = C_est,       
    a_norms     = a_norms,
    b_norms     = b_norms
  )
}

estimate_C_binary <- function(X,
                              qmax   = 8,     # truncated 
                              safety = 1.5,   
                              eps    = 1e-12,  # smooth
                              radius = 1     
) {
  
  n <- nrow(X); p <- ncol(X)
  
  P_hat    <- matrix(NA_real_, n, p)
  N_counts <- matrix(0L,      n, p)
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      window_values <- c(X[i, j])

      for (d in seq_len(radius)) {
        if (i - d >= 1)       window_values <- c(window_values, X[i - d, j])
        if (i + d <= n)       window_values <- c(window_values, X[i + d, j])
      }

      for (d in seq_len(radius)) {
        if (j - d >= 1)       window_values <- c(window_values, X[i, j - d])
        if (j + d <= p)       window_values <- c(window_values, X[i, j + d])
      }

      vals <- window_values[!is.na(window_values)]
      P_hat[i, j]    <- mean(vals)
      N_counts[i, j] <- length(vals)
    }
  }
  

  P_adj <- (P_hat * N_counts + eps) / (N_counts + 5*eps)
  Mhat  <- log(P_adj / (1 - P_adj))
  

  svd_res <- irlba::irlba(Mhat, nv = qmax, nu = qmax)
  U <- svd_res$u     # n×qmax
  V <- svd_res$v     # p×qmax
  d <- svd_res$d[1:qmax]
  Dsqrt <- diag(sqrt(d), qmax, qmax)
  A_hat <- U %*% Dsqrt  # n×qmax
  B_hat <- V %*% Dsqrt  # p×qmax
  

  a_norms <- sqrt(rowSums(A_hat^2))
  b_norms <- sqrt(rowSums(B_hat^2))
  

  C0    <- max(c(a_norms, b_norms))
  C_est <- C0 * safety
  
  list(
    radius   = radius,
    qmax     = qmax,
    safety   = safety,
    C0       = C0,
    C_est    = C_est,
    a_norms  = a_norms,
    b_norms  = b_norms,
    Mhat     = Mhat,
    P_smooth = P_hat,
    N_counts = N_counts
  )
}


sim1_bin <- function(N, J, K){
  F.matr <- matrix(runif(N*K,-8,8), N, K);
  F.matr <- cbind(rep(1,N), F.matr)
  A.matr <- matrix(runif(J*(K+1),-8,8), J, K+1);                 
  
  temp = F.matr %*% t(A.matr);
  prob = 1/(1+exp(-temp));
  resp = matrix(0, N, J);
  resp[] = rbinom(N*J, 1, prob)
  
  list(resp = resp, F.matr = F.matr[,-1], A.matr = A.matr[,-1], d.vec = A.matr[,1], M = temp)
}



factor_selection_cv_binary <- function(resp, theta0, A0, d0, K = 5, max_factors = 8, 
                                       C = 5, tol_factor = 0.01, see = 123) {  
  if (ncol(theta0) < max_factors) {
    stop("theta0 at least max_factors")
  }
  if (ncol(A0) < max_factors) {
    stop("A0 at least max_factors")
  }
  
  n <- nrow(resp)
  p <- ncol(resp)
  

  set.seed(see)
  fold_indices <- sample(rep(1:K, length.out = n))
  
  cv_losses <- matrix(0, nrow = K, ncol = max_factors)

  for (k in 1:K) {
    #cat("Processing fold", k, "of", K, "\n")
    
    test_idx <- which(fold_indices == k)
    train_idx <- which(fold_indices != k)
    
    X_train <- resp[train_idx, , drop = FALSE]
    X_test <- resp[test_idx, , drop = FALSE]
    num_train <- nrow(X_train) * p

    theta0_train <- theta0[train_idx, , drop = FALSE]

    true_matrix <- matrix(TRUE, nrow = nrow(X_train), ncol = ncol(X_train))
    
    n_test <- nrow(X_test)
    

    for (q in 1:max_factors) {
     # cat("  Testing q =", q, "\n")

      jml_result <- tryCatch({
        mirtjml_expr(
          X_train, q, 
          theta0 = theta0_train[, 1:q, drop = FALSE], 
          A0 = A0[, 1:q, drop = FALSE], 
          d0 = d0,
          tol = 0.01 * num_train, 
          cc = C,
          print_proc = FALSE
        )
      }, error = function(e) {
        #cat("    Error in mirtjml_expr for q =", q, ":", e$message, "\n")
        return(NULL)
      })

      if (is.null(jml_result)) {
        cv_losses[k, q] <- Inf
        next
      }
      
      alpha <- jml_result$d_hat 
      B <- jml_result$A_hat

      B_cond <- kappa(B)
      if (B_cond > 1e10) {
        #cat("    Warning: Loading matrix is ill-conditioned (kappa =", B_cond, ")\n")
      }

      F_test <- matrix(0, nrow = n_test, ncol = q)
      failed_samples <- 0

      for (i in 1:n_test) {
        y_i <- X_test[i, ]

        if (all(y_i == 0) || all(y_i == 1)) {
          F_test[i, ] <- rep(0, q)
          failed_samples <- failed_samples + 1
          next
        }
        

        glm_success <- FALSE

        glm_data <- data.frame(
          y = y_i,
          alpha = alpha
        )
        
        for (j in 1:q) {
          glm_data[[paste0("B", j)]] <- B[, j]
        }
        

        formula_str <- paste("y ~ offset(alpha) +", 
                             paste(paste0("B", 1:q), collapse = " + "), 
                             "- 1")

        glm_result <- tryCatch({
          glm(as.formula(formula_str), 
              data = glm_data, 
              family = binomial(link = "logit"),
              control = list(epsilon = 1e-6, maxit = 50))
        }, warning = function(w) {

          suppressWarnings(glm(as.formula(formula_str), 
                               data = glm_data, 
                               family = binomial(link = "logit"),
                               control = list(epsilon = 1e-6, maxit = 50)))
        }, error = function(e) {
          NULL
        })

        if (!is.null(glm_result) && glm_result$converged) {
          coefficients <- coef(glm_result)
          if (!any(is.na(coefficients)) && !any(is.infinite(coefficients)) && 
              all(abs(coefficients) < 1e10)) {
            F_test[i, ] <- coefficients
            glm_success <- TRUE
          }
        }

        if (!glm_success) {
          F_test[i, ]=rep(0,q)
        }
      }
      
      if (failed_samples > 0) {
       # cat("    Warning:", failed_samples, "samples failed in factor score estimation\n")
      }
      

      eta <- F_test %*% t(B) + matrix(1, nrow = n_test, ncol = 1) %*% t(alpha)

      eta <- pmax(eta, -20)
      eta <- pmin(eta, 20)

      prob_matrix <- plogis(eta)
      
      likelihood_loss <- -sum(dbinom(X_test, size = 1, prob = prob_matrix, log = TRUE))

      if (is.na(likelihood_loss) || is.infinite(likelihood_loss)) {
        #cat("    Warning: Invalid likelihood loss for q =", q, "\n")
        cv_losses[k, q] <- Inf
      } else {
        cv_losses[k, q] <- likelihood_loss
      }
    }
  }

  mean_cv_losses <- colMeans(cv_losses)

  valid_losses <- which(is.finite(mean_cv_losses))
  
  if (length(valid_losses) == 0) {
    #warning("No valid models found. Check your data and parameters.")
    optimal_q <- 1
  } else {
    optimal_q <- valid_losses[which.min(mean_cv_losses[valid_losses])]
  }
  

  result <- list(
    optimal_factors = optimal_q,
    cv_losses = cv_losses,
    mean_cv_losses = mean_cv_losses,
    fold_indices = fold_indices,
    n_valid_models = length(valid_losses)
  )
  
  #cat("\nOptimal number of factors:", optimal_q, "\n")
  #cat("Valid models found for q =", paste(valid_losses, collapse = ", "), "\n")
  
  return(result)
}




sim1_depen <- function(N, J, K){

  F.matr <- matrix(runif(N*K, -2, 2), N, K)
  F.matr <- cbind(rep(1, N), F.matr)

  A.matr <- matrix(runif(J*(K+1), -2, 2), J, K+1)
  temp = F.matr %*% t(A.matr)

  Sigma <- matrix(0, J, J)
  

  diag(Sigma) <- 0.05

  for(i in 1:(J-1)){
    Sigma[i, i+1] <- 0.01
    Sigma[i+1, i] <- 0.01
  }
  

  noise <- mvrnorm(n = N, mu = rep(0, J), Sigma = Sigma)

  temp <- temp + noise
  

  resp <- matrix(0, N, J)
  resp[] <- rpois(N*J, lambda = exp(temp))
  
  list(resp = resp, F.matr = F.matr, A.matr = A.matr, M = temp, 
       noise = noise, Sigma = Sigma)
}


sim1_overdis <- function(N, J, K){

  F.matr <- matrix(runif(N*K, -2, 2), N, K)
  F.matr <- cbind(rep(1, N), F.matr)

  A.matr <- matrix(runif(J*(K+1), -2, 2), J, K+1)

  temp = F.matr %*% t(A.matr)

  lambda <- exp(temp)

  resp <- matrix(0, N, J)

  random_matrix <- matrix(runif(N*J), N, J)

  selection_threshold1 <- quantile(random_matrix, 0.1) 
  selection_threshold2 <- quantile(random_matrix, 0.2) 

  for (i in 1:N) {
    for (j in 1:J) {
      if (random_matrix[i, j] <= selection_threshold1) {

        theta <- 1
        size <- theta
        prob <- size / (size + lambda[i, j])
        resp[i, j] <- rnbinom(1, size, prob)
      } else if (random_matrix[i, j] <= selection_threshold2) {

        theta <- 10
        size <- theta
        prob <- size / (size + lambda[i, j])
        resp[i, j] <- rnbinom(1, size, prob)
      } 
      else {

        resp[i, j] <- rpois(1, lambda[i, j])
      }
    }
  }
  
  list(resp = resp, F.matr = F.matr, A.matr = A.matr, M = temp)
}


