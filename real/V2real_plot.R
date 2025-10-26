#############################################################################################
# This script corresponds to real data Figure 1 in main text and Figure S3 in the supplementary material.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################


# The mouse brain real dataset was downloaded from:
# https://www.kaggle.com/datasets/aayush9753/singlecell-rnaseq-data-from-mouse-brain

# loading packages
library(pheatmap)
library(Seurat)
library(irlba)
library(readr)
library(GFM)
library(Rcpp)
library(RcppArmadillo)


### Function for identifiability
Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}
add_identifiability <- function(H, B, mu){
  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  signB1 <- sign(svdHB$v[1,])
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)
  
  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)
  
  return(list(H=H, B=B, mu=mu))
}

s=666
set.seed(s)  

save_path <- paste0("~/glmcv/realV2/gene/V2real_plot", s, "/")

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

### loading data
brain_metadata <- read_csv("~/glmcv/realV2/gene/brain_metadata.csv")
brain_counts <- read_csv("~/glmcv/realV2/gene/brain_counts.csv")
colnames(brain_metadata)[1] <- "ID"
colnames(brain_counts)[1] <- "ID"

merged_data <- merge(brain_metadata, brain_counts, by = "ID", all = TRUE)
filtered_data <- merged_data[merged_data$mouse.sex == "F", ]
brain_count_f <- filtered_data[, -c(1:6)]

if (any(is.na(brain_count_f))) {
  cat("move NA ...\n")
  brain_count_f <- brain_count_f[, colSums(is.na(brain_count_f)) == 0] 
}

if (any(colSums(brain_count_f) == 0)) {
  cat("move zero gene...\n")
  brain_count_f <- brain_count_f[, colSums(brain_count_f) > 0]  
}

if (any(rowSums(brain_count_f) == 0)) {
  cat("move zero cells...\n")
  brain_count_f <- brain_count_f[rowSums(brain_count_f) > 0, ]  
}
# print("dim(brain_count_f)")
# print(dim(brain_count_f))
# 707 18439

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
# print(dim(brain_count_f1000))  


n=dim(brain_count_f1000)[1]
p=dim(brain_count_f1000)[2]

mat_es=matrix(0, nrow = 5, ncol = 10)


source("~/glmcv/simuV2/glmcv_functions.R")
Rcpp::sourceCpp("~/glmcv/simuV2/confirm_jmle_omp_poisson_with_intercept_missing.cpp")

resp = as.matrix(brain_count_f1000)

true_matrix <- matrix(TRUE, nrow = n, ncol = p)


est_c=estimate_C(log(resp+1),qmax = 6)

init=GFM::gfm(list(resp),types = 'poisson',q=15,verbose = F)
theta0=cbind(rep(1,n),init$hH)
A0 <- cbind(init$hmu,init$hB)

C_ada=(est_c$C_est)


jic.res15 = confirm_CJMLE_poisson_cpp(as.matrix(resp),  true_matrix, theta0, A0, matrix(TRUE, p, 16), C = C_ada, tol = 0.01/n/p)

saveRDS(jic.res15, file = file.path(save_path, "jicres15.rds"))





######################################### plot
##### loading genes Figure S3 
all.genes <- rownames(seurat_obj)
pcadata=ScaleData(FindVariableFeatures(NormalizeData(seurat_obj), selection.method = "vst", nfeatures = 1000), features = all.genes)
seurat_objpca <- RunPCA(pcadata, features = VariableFeatures(object = seurat_obj),  npcs = 15)

jicres15 <- readRDS(file.path(save_path, "jicres15.rds"))
loading=jicres15$A[,-1]
embedding=jicres15$theta[,-1]


para_add=add_identifiability(embedding,loading,jicres15$A[,1])

loadings <- para_add$B
rownames(loadings)=rownames(seurat_objpca@reductions$pca@feature.loadings)
colnames(loadings)=colnames(seurat_objpca@reductions$pca@feature.loadings)
seurat_objpca@reductions$pca@feature.loadings=loadings
embeddings <- para_add$H
rownames(embeddings)=rownames(seurat_objpca@reductions$pca@cell.embeddings)
colnames(embeddings)=colnames(seurat_objpca@reductions$pca@cell.embeddings)
seurat_objpca@reductions$pca@cell.embeddings=embeddings


variances <- apply(loadings, 2, var)
stdev <- sqrt(variances)
seurat_objpca@reductions$pca@stdev=stdev

seurat_objpca@reductions$pca@key <- "Factor"


for (i in 1:15) {
  filename <- paste0(save_path, "load", i, ".pdf")
  
  tryCatch({
    pdf(filename, width = 8, height = 8)
    print(VizDimLoadings(seurat_objpca, dims = i, reduction = "pca"))

    dev.off()
    
    cat("Save successfully ", filename, "\n")
    
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    cat("error", filename, "save", e$message, "\n")
  })
}
## give the important feature
print(seurat_objpca[["pca"]], dims = 1:15, nfeatures = 10)


########### heatmap plot Figure 1
log_transformed_data <- log1p(t(brain_count_f1000))
pcadata=ScaleData(FindVariableFeatures(NormalizeData(seurat_obj), selection.method = "vst", nfeatures = 1000), features = all.genes)

variable_genes <- VariableFeatures(pcadata)
final_expression_data <- GetAssayData(pcadata, slot = "data")[variable_genes, ]

for (i in 1:15) {
  load <- loadings[, i]
  fac <- embeddings[, i]

  top20_names <- rownames(log_transformed_data)[order(load, decreasing = TRUE)[1:20]]
  bottom20_names <- rownames(log_transformed_data)[order(load, decreasing = FALSE)[1:20]]
  selected_gene_names <- unique(c(top20_names, bottom20_names))

  top50_names <- colnames(log_transformed_data)[order(fac, decreasing = TRUE)[1:150]]
  bottom50_names <- colnames(log_transformed_data)[order(fac, decreasing = FALSE)[1:150]]
  selected_cell_names <- unique(c(top50_names, bottom50_names))

  selected_datafinal <- (as.matrix(final_expression_data))[selected_gene_names, selected_cell_names]

  filename <- paste0(save_path, "cell", i, ".pdf")

  pdf(filename, width = 8, height = 8)

  pheatmap(
    selected_datafinal,
    main = paste0("Factor ", i), 
    legend = FALSE,
    show_colnames = FALSE,
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "azure", "red"))(1000),
    na_col = "azure",
    fontsize_row = 7.5,   
    fontsize = 15,      
    fontsize_main = 60  
  )

  dev.off()
  cat("Save", filename, " (dim=", dim(selected_datafinal)[1], "×", dim(selected_datafinal)[2], ")\n")
}






