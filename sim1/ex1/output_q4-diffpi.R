
#############################################################################################
# This script corresponds to reproduce Table S1 in the supplement materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

# -----------------------------------------------------------------------------
base_folder <- "~/glmcv/simuV2/result_poi/V2q4/"  
output_folder <- "~/glmcv/simuV2/result_poi/V2q4/" 
file_prefix <- "facnum-cada"  


n_values <- c(50, 100, 300)  
p_values <- c(50, 100, 300)  
q_value <- 4  



target_columns <- c(
  "ES_K2" ,    "ES_K2_pen3" , "ES_K2_pen4" , "ES_K2_pen1", "ES_K2_pen2" ,  "GCV_K2",              
  "GCV_K2_pen3" ,  "GCV_K2_pen4" ,  "GCV_K2_pen1" ,"GCV_K2_pen2", "ES_K5"  ,            
  "ES_K5_pen3" ,  "ES_K5_pen4",  "ES_K5_pen1",   "ES_K5_pen2",  "GCV_K5" ,         
  "GCV_K5_pen3",  "GCV_K5_pen4",   "GCV_K5_pen1",  "GCV_K5_pen2",    
  "ES_K10",   "ES_K10_pen3", "ES_K10_pen4",     "ES_K10_pen1",      "ES_K10_pen2",   "GCV_K10",             
  "GCV_K10_pen3",   "GCV_K10_pen4" ,   "GCV_K10_pen1",     "GCV_K10_pen2",  "ES_rev_K5",      
  "ES_rev_K5_pen3",   "ES_rev_K5_pen4" ,  "ES_rev_K5_pen1",  "ES_rev_K5_pen2", "GCV_rev_K5" ,     
  "GCV_rev_K5_pen3",  "GCV_rev_K5_pen4",  "GCV_rev_K5_pen1",  "GCV_rev_K5_pen2",
  "ES_rev_K10",       "ES_rev_K10_pen3", "ES_rev_K10_pen4",   "ES_rev_K10_pen1",  
  "ES_rev_K10_pen2", "GCV_rev_K10", "GCV_rev_K10_pen3","GCV_rev_K10_pen4", "GCV_rev_K10_pen1", "GCV_rev_K10_pen2"
)

# Note: pen3 and pen4 correspond to f1 and f2 in the main text, while pen1 and pen2 correspond to f3 and f4.


# -----------------------------------------------------------------------------
file_info <- expand.grid(n = n_values, p = p_values, stringsAsFactors = FALSE)
file_info$q <- q_value


file_info$filename <- paste0("n", file_info$n, "_p", file_info$p, "_q", file_info$q, "_", file_prefix, ".rds")
file_info$filepath <- file.path(base_folder, file_info$filename)

print("files:")
print(file_info[, c("n", "p", "q", "filename")])


# -----------------------------------------------------------------------------
missing_files <- file_info$filepath[!file.exists(file_info$filepath)]
if(length(missing_files) > 0) {
  cat("no files \n")
  cat(paste(missing_files, collapse = "\n"))
  cat("\n")
}

existing_files <- file_info[file.exists(file_info$filepath), ]
cat(paste("find ", nrow(existing_files), "files\n"))


# -----------------------------------------------------------------------------
process_single_file <- function(filepath, n_val, p_val, q_val, target_cols) {
  cat("For ", basename(filepath), "\n")
  
  data <- readRDS(filepath)
  
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  } else if(is.matrix(data)) {
  } else {
    cat("unknown type", class(data), "\n")
    data <- as.matrix(data)
  }
  

  if(is.null(colnames(data))) {
    stop("No column")
  }
  

  missing_cols <- target_cols[!target_cols %in% colnames(data)]
  if(length(missing_cols) > 0) {
    cat("file", basename(filepath), "lacks", paste(missing_cols, collapse = ", "), "\n")
  }
  

  available_cols <- target_cols[target_cols %in% colnames(data)]
  extracted_data <- data[, available_cols, drop = FALSE]
  

  col_sums <- apply(extracted_data, 2, function(x) {
    if(all(is.na(x) | x == 0)) return(NA)
    else return(sum(x, na.rm = TRUE))
  })
  

  result_vector <- sapply(available_cols, function(col_name) {
    col_data <- extracted_data[, col_name]  
    
    if(all(is.na(col_data) | col_data == 0)) {
      return("/")
    }
    
    valid_data <- col_data[!is.na(col_data) & col_data != 0]
    
    if(length(valid_data) == 0) {
      return("/")
    }
    
    under_count <- sum(valid_data < q_val)
    over_count <- sum(valid_data > q_val)
    total_count <- length(valid_data)
    
    under_freq <- under_count 
    over_freq <- over_count
    
    return(paste0(under_freq, "/", over_freq))
  })

  final_result <- c(q_val, n_val, p_val, result_vector)
  names(final_result) <- c("q", "n", "p", available_cols)
  
  return(list(result = final_result, col_sums = col_sums, filename = basename(filepath)))
}


# -----------------------------------------------------------------------------
cat("\n Start...\n")
all_results <- list()
all_col_sums <- list()

for(i in 1:nrow(existing_files)) {
  file_result <- process_single_file(
    existing_files$filepath[i],
    existing_files$n[i],
    existing_files$p[i],
    existing_files$q[i],
    target_columns
  )
  
  all_results[[i]] <- file_result$result
  all_col_sums[[i]] <- file_result$col_sums
  names(all_col_sums)[i] <- file_result$filename
}


# -----------------------------------------------------------------------------

final_colnames <- c("q", "n", "p", target_columns)


result_facnum_mat <- matrix("", nrow = length(all_results), ncol = length(final_colnames))
colnames(result_facnum_mat) <- final_colnames


for(i in 1:length(all_results)) {
  result_row <- all_results[[i]]

  for(col in final_colnames) {
    if(col %in% names(result_row)) {
      result_facnum_mat[i, col] <- as.character(result_row[col])
    } else {
      result_facnum_mat[i, col] <- "/"
    }
  }
}


result_df <- as.data.frame(result_facnum_mat, stringsAsFactors = FALSE)
result_df$q <- as.numeric(result_df$q)
result_df$n <- as.numeric(result_df$n)
result_df$p <- as.numeric(result_df$p)


result_df <- result_df[order(result_df$q, result_df$n, result_df$p), ]


result_facnum_mat <- as.matrix(result_df)

# cat("\n dim", dim(result_facnum_mat), "\n")
# cat("column name: \n")
# print(colnames(result_facnum_mat))



if(!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


rds_filename <- file.path(output_folder, "result_facnum_mat-diffpi.rds")
saveRDS(result_facnum_mat, rds_filename)
cat("RDS save as", rds_filename, "\n")

library(openxlsx)
xlsx_filename <- file.path(output_folder, "result_facnum_mat-diffpi.xlsx")
write.xlsx(result_facnum_mat, xlsx_filename, rowNames = FALSE)
cat("Excel save as", xlsx_filename, "\n")
