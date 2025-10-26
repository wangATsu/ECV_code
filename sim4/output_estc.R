#############################################################################################
# This script corresponds to reproduce Figure S1 (left plot) in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

data_folder <- "~/glmcv/simuV2/estimate_c/"

n_val <- 50
p_val <- 50
q_val <- 4
C_vals <- 1:10  

file_pattern <- "n{n}_p{p}_q{q}_C{C}_frequency.rds"

selected_columns <- c(
   "ES_K5", "ES_K5_pen1", "ES_K5_pen2",
  "ES_K5_pen3", "ES_K5_pen4", "GCV_K5", "GCV_K5_pen1",
  "GCV_K5_pen2", "GCV_K5_pen3", "GCV_K5_pen4"
)

output_rds_name <- "result_facnum_mat_np50k5.rds"
output_xlsx_name <- "result_facnum_mat_np50k5.xlsx"

# ======================================
library(openxlsx)

result_list <- list()
column_sums_info <- list()

for(C in C_vals) {

  filename <- glue::glue(file_pattern, n = n_val, p = p_val, q = q_val, C = C)
  filepath <- file.path(data_folder, filename)
  
 # cat("file ", filename, "\n")

  if(!file.exists(filepath)) {
    cat("no file", filepath, "\n")
    next
  }
  
  data_matrix <- readRDS(filepath)

  if(!is.matrix(data_matrix)) {
    cat("not matrix", filename, "\n")
    next
  }
  
  if(!all(selected_columns %in% colnames(data_matrix))) {
    missing_cols <- selected_columns[!selected_columns %in% colnames(data_matrix)]
    cat("file ", filename, "lack ", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  selected_data <- data_matrix[, selected_columns, drop = FALSE]

  result_vector <- c()
  col_sums <- c()
  
  for(col_name in selected_columns) {
    col_data <- selected_data[, col_name]

    col_sum <- sum(col_data, na.rm = TRUE)
    col_sums <- c(col_sums, col_sum)

    if(all(is.na(col_data)) || all(col_data == 0, na.rm = TRUE)) {
      result_vector <- c(result_vector, "/")
    } else {

      under_count <- sum(col_data[1:(q_val-1)], na.rm = TRUE)

      over_count <- sum(col_data[(q_val+1):8], na.rm = TRUE)
      
      result_string <- paste0(under_count, "/", over_count)
      result_vector <- c(result_vector, result_string)
    }
  }

  final_vector <- c(C, n_val, p_val, result_vector)

  result_list[[length(result_list) + 1]] <- final_vector

  names(col_sums) <- selected_columns
  column_sums_info[[paste0("C", C)]] <- col_sums
}

  result_facnum_mat <- do.call(rbind, result_list)

  colnames(result_facnum_mat) <- c("C", "n", "p", selected_columns)

  result_df <- as.data.frame(result_facnum_mat, stringsAsFactors = FALSE)
  result_df$C <- as.numeric(result_df$C)
  result_df$n <- as.numeric(result_df$n)
  result_df$p <- as.numeric(result_df$p)

  result_df <- result_df[order(result_df$C, result_df$n, result_df$p), ]

  result_facnum_mat <- as.matrix(result_df)

  output_rds_path <- file.path(data_folder, output_rds_name)
  output_xlsx_path <- file.path(data_folder, output_xlsx_name)
  
  saveRDS(result_facnum_mat, output_rds_path)
  write.xlsx(result_facnum_mat, output_xlsx_path, rowNames = FALSE)
  


