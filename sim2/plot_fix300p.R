#############################################################################################
# This script corresponds to reproduce Figure S2 (left plot) in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

# load packages
library(readr)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)


base_dir <- "~/glmcv/simuV2/result_poi/V2fix300p"

# Here, the GCV formulation corresponds to the ECV method in the paper.
# Note: pen3 and pen4 correspond to f1 and f2 in the main text,
# while pen1 and pen2 correspond to f3 and f4.
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

  # cat("Processing file:", basename(file_path), "\n")
  # cat("Matrix dimensions:", dim(matrix_data), "\n")
  # 
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


combinations <- expand.grid(
  n = c(30, 50, 75, 100, 150, 200, 300,400,500,600,700, 800,900, 1000), 
  p = c(300)
)


results <- list()


for(i in 1:nrow(combinations)) {
  n <- combinations$n[i]
  p <- combinations$p[i]
  
  file_name <- sprintf("n%d_p%d_q4_facnum-fixp-cada.rds", n, p)
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


write_xlsx(final_results_sorted, path = file.path(base_dir, "result_rate_300p.xlsx"))

result_rate_300p <- final_results_sorted[, c("n", "GCV_rev_K10", "GCV_rev_K5", "GCV_K2", "GCV_K5", "GCV_K10")]
View(result_rate_300p)
write_xlsx(result_rate_300p, path = file.path(base_dir, "result_rate_300p_plot.xlsx"))




########## code for plot 



data_long <- tidyr::pivot_longer(
  result_rate_300p,
  cols = -n,  
  names_to = "Method",
  values_to = "Rate"
)


p <- ggplot(data_long, aes(x = n, y = Rate, color = Method, shape = Method)) +
  geom_line() +                     
  geom_point(size = 2) +            
  scale_color_manual(
    breaks = c("GCV_rev_K10", "GCV_rev_K5", "GCV_K2", "GCV_K5", "GCV_K10"),
    values = c("red", "blue", "green", "purple", "orange"),
    labels = c(expression(pi == 0.1), expression(pi == 0.2), 
               expression(pi == 0.5), expression(pi == 0.8), 
               expression(pi == 0.9))
  ) +
  scale_shape_manual(
    breaks = c("GCV_rev_K10", "GCV_rev_K5", "GCV_K2", "GCV_K5", "GCV_K10"),
    values = c(16, 17, 15, 18, 19),
    labels = c(expression(pi == 0.1), expression(pi == 0.2), 
               expression(pi == 0.5), expression(pi == 0.8), 
               expression(pi == 0.9))
  ) +
  scale_x_continuous(
    breaks = unique(data_long$n)  
  ) +
  labs(x = "n",                   
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

# save pdf
ggsave(
  filename = file.path(base_dir, "correct_selection_rate_plot_300p.pdf"),
  plot = p,
  device = "pdf",
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)
