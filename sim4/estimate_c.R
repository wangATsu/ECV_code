#############################################################################################
# This script corresponds to reproduce Figure S1 (right plot) in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################

library(openxlsx)
source("~/glmcv/simuV2/glmcv_functions.R")
dir.create("~/glmcv/simuV2/estimate_c", recursive = TRUE, showWarnings = FALSE)

#### get the estimators
n_values <- c(300)
p_values <- c(50, 100, 300)
S=100
for (n in n_values) {
  estc_matrix <- matrix(NA, nrow = S, ncol = length(p_values))
  colnames(estc_matrix) <- paste0("p=", p_values)
  
  for (p_idx in 1:length(p_values)) {
    p <- p_values[p_idx]
    
    for (s in 1:S) {
      set.seed(s)
      resppoi <- sim1(n, p, 4)$resp
      estc_matrix[s, p_idx] <- estimate_C(log(resppoi + 1),qmax=9)$C_est
      
      cat("n =", n, " p =", p, " s =", s, " C est =", estc_matrix[s, p_idx], "\n")
    }
  }

  filename <- paste0("~/glmcv/simuV2/estimate_c/estc_poi_n", n)
  saveRDS(estc_matrix, paste0(filename, ".rds"))
  write.xlsx(as.data.frame(estc_matrix), paste0(filename, ".xlsx"), rowNames = TRUE)
}


library(ggplot2)
library(tidyr)
library(dplyr)


######## code for plot

estc_matrix <- readRDS("~/glmcv/simuV2/estimate_c/estc_poi_n300.rds")

estc_df <- as.data.frame(estc_matrix)
estc_df$s <- 1:nrow(estc_df)

estc_long <- estc_df %>%
  pivot_longer(cols = -s, names_to = "p_value", values_to = "cad")

estc_long$p_value <- factor(estc_long$p_value, 
                            levels = c("p=50", "p=100", "p=300"))

colors <- c("p=50" = "blue", "p=100" = "green", "p=300" = "red")

p <- ggplot(estc_long, aes(x = p_value, y = cad, color = p_value)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  geom_hline(yintercept = sqrt(20), linetype = "dashed", color = "black") +
  scale_color_manual(values = colors) +
  scale_y_continuous(limits = c(3, 10)) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "none", 
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 12, color = "black")  
  ) +
  labs(y = "Data-driven tuning of C")

print(p)

ggsave("~/glmcv/simuV2/estimate_c/estc_poi_n300_boxplot.pdf", 
       plot = p, width = 9, height = 7, dpi = 300)

