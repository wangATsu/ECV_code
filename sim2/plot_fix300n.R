
#############################################################################################
# This script corresponds to reproduce Figure S2 (right plot) in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

# loading packages
library(readr)
library(dplyr)
library(tidyr)
library(writexl)


base_dir <- "~/glmcv/simuV2/result_poi/V2fix300n"

actual_column_names <- c("ES_K2", "ES_K2_pen1", "ES_K2_pen2", "ES_K2_pen3", "ES_K2_pen4",
                         "GCV_K2", "GCV_K2_pen1", "GCV_K2_pen2", "GCV_K2_pen3", "GCV_K2_pen4",
                         "ES_K5", "ES_K5_pen1", "ES_K5_pen2", "ES_K5_pen3", "ES_K5_pen4",
                         "GCV_K5", "GCV_K5_pen1", "GCV_K5_pen2", "GCV_K5_pen3", "GCV_K5_pen4",
                         "ES_K10", "ES_K10_pen1", "ES_K10_pen2", "ES_K10_pen3", "ES_K10_pen4",
                         "GCV_K10", "GCV_K10_pen1", "GCV_K10_pen2", "GCV_K10_pen3", "GCV_K10_pen4",
                         "ES_rev_K5", "ES_rev_K5_pen1", "ES_rev_K5_pen2", "ES_rev_K5_pen3", "ES_rev_K5_pen4",
                         "GCV_rev_K5", "GCV_rev_K5_pen1", "GCV_rev_K5_pen2", "GCV_rev_K5_pen3", "GCV_rev_K5_pen4",
                         "ES_rev_K10", "ES_rev_K10_pen1", "ES_rev_K10_pen2", "ES_rev_K10_pen3", "ES_rev_K10_pen4",
                         "GCV_rev_K10", "GCV_rev_K10_pen1", "GCV_rev_K10_pen2", "GCV_rev_K10_pen3", "GCV_rev_K10_pen4",
                         "JIC1", "JIC2", "AIC1", "AIC2", "BIC1", "BIC2", "BIC3",
                         "BCV_Gabriel", "BCV_Wold", "Ladle", "Ratio", "SampleCV2", "SampleCV_Bayes")


process_file <- function(n, p, file_path) {

  matrix_data <- readRDS(file_path)
  
  cat("Processing file:", basename(file_path), "\n")
  cat("Matrix dimensions:", dim(matrix_data), "\n")
  if(!is.matrix(matrix_data)) {
    matrix_data <- as.matrix(matrix_data)
  }
  
  rightselection_counts <- colSums(matrix_data == 4) / nrow(matrix_data)

  if(length(rightselection_counts) != length(actual_column_names)) {
    warning(paste("Column count mismatch! Expected:", length(actual_column_names), 
                  "Got:", length(rightselection_counts)))
  }

  rightselection_df <- data.frame(
    selection = "rightselection", 
    n = n, 
    p = p,
    stringsAsFactors = FALSE
  )

  for(i in 1:length(actual_column_names)) {
    if(i <= length(rightselection_counts)) {
      rightselection_df[[actual_column_names[i]]] <- rightselection_counts[i]
    } else {
      rightselection_df[[actual_column_names[i]]] <- NA
    }
  }
  
  return(rightselection_df)
}

process_file_all_stats <- function(n, p, file_path) {

  matrix_data <- readRDS(file_path)
  
  # cat("Processing file:", basename(file_path), "\n")
  # cat("Matrix dimensions:", dim(matrix_data), "\n")
  # cat("Unique values in matrix:", sort(unique(as.vector(matrix_data))), "\n\n")
  # 

  if(!is.matrix(matrix_data)) {
    matrix_data <- as.matrix(matrix_data)
  }
  
  underselection_counts <- colSums(matrix_data >= 1 & matrix_data <= 3) / nrow(matrix_data)
  rightselection_counts <- colSums(matrix_data == 4) / nrow(matrix_data)
  overselection_counts <- colSums(matrix_data >= 5 & matrix_data <= 8) / nrow(matrix_data)

  results_list <- list()
  
  underselection_df <- data.frame(
    selection = "underselection", 
    n = n, 
    p = p,
    stringsAsFactors = FALSE
  )
  for(i in 1:length(actual_column_names)) {
    if(i <= length(underselection_counts)) {
      underselection_df[[actual_column_names[i]]] <- underselection_counts[i]
    }
  }
  results_list[[1]] <- underselection_df
  
  rightselection_df <- data.frame(
    selection = "rightselection", 
    n = n, 
    p = p,
    stringsAsFactors = FALSE
  )
  for(i in 1:length(actual_column_names)) {
    if(i <= length(rightselection_counts)) {
      rightselection_df[[actual_column_names[i]]] <- rightselection_counts[i]
    }
  }
  results_list[[2]] <- rightselection_df
  
  overselection_df <- data.frame(
    selection = "overselection", 
    n = n, 
    p = p,
    stringsAsFactors = FALSE
  )
  for(i in 1:length(actual_column_names)) {
    if(i <= length(overselection_counts)) {
      overselection_df[[actual_column_names[i]]] <- overselection_counts[i]
    }
  }
  results_list[[3]] <- overselection_df
  
  return(bind_rows(results_list))
}

combinations <- expand.grid(n = c(300), p = c(30, 50, 75, 100, 150, 200, 300,400,500,600,700, 800,900, 1000))

results <- list()

for(i in 1:nrow(combinations)) {
  n <- combinations$n[i]
  p <- combinations$p[i]
  
  file_name <- sprintf("n%d_p%d_q4_facnum-fixn-cada.rds", n, p)
  file_path <- file.path(base_dir, file_name)
  
  if(file.exists(file_path)) {
    tryCatch({
      result <- process_file(n, p, file_path)
      
      results[[length(results) + 1]] <- result
    }, error = function(e) {
      # cat("Error processing file:", file_path, "\n")
      # cat("Error message:", e$message, "\n\n")
    })
  } else {
    cat("File not found:", file_path, "\n")
  }
}


final_results <- bind_rows(results)


final_results_sorted <- final_results %>%
  arrange(desc(selection), n, p)



write_xlsx(final_results_sorted, path = file.path(base_dir, "result_rate_300n.xlsx"))
result_rate_300n=final_results_sorted[, c("p", "GCV_rev_K10","GCV_rev_K5","GCV_K2","GCV_K5","GCV_K10")]
View(result_rate_300n)
write_xlsx(result_rate_300n, path = file.path(base_dir, "result_rate_300n_plot.xlsx"))


################# plot 

library(ggplot2)
data_long <- tidyr::pivot_longer(
  result_rate_300n,
  cols = -p,  
  names_to = "Method",
  values_to = "Rate"
)





p <-ggplot(data_long, aes(x = p, y = Rate, color = Method, shape = Method)) +
  geom_line() +                    
  geom_point(size = 2) +            
  scale_color_manual(
    breaks = c("GCV_rev_K10", "GCV_rev_K5", "GCV_K2", "GCV_K5", "GCV_K10"),  # 添加这一行
    values = c("red", "blue", "green", "purple", "orange"),
    labels = c(expression(pi == 0.1), expression(pi == 0.2), 
               expression(pi == 0.5), expression(pi == 0.8), 
               expression(pi == 0.9))
  ) +
  scale_shape_manual(
    breaks = c("GCV_rev_K10", "GCV_rev_K5", "GCV_K2", "GCV_K5", "GCV_K10"),  # 添加这一行
    values = c(16, 17, 15, 18, 19),
    labels = c(expression(pi == 0.1), expression(pi == 0.2), 
               expression(pi == 0.5), expression(pi == 0.8), 
               expression(pi == 0.9))
  ) +
  scale_x_continuous(
    breaks = unique(data_long$p)  
  ) +  
  labs(x = "p",                    
       y = "Correct selection rate") + 
  theme_classic() +                
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),        
    legend.key.size = unit(1.2, "cm")             
  )



ggsave(
  filename = file.path(base_dir, "correct_selection_rate_plot.pdf"),
  plot = p,
  device = "pdf",
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)
