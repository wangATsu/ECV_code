#############################################################################################
# This script corresponds to reproduce Figure S1 (left plot) in the supplementary materials.
# Before running, change the relevant paths to the corresponding locations on your computer,
# and place the related functions in the appropriate directory for easy sourcing.
#############################################################################################


# loading packages
library(ggplot2)
library(tidyr)
library(dplyr)

# Note: pen3 and pen4 correspond to f1 and f2 in the main text,
# while pen1 and pen2 correspond to f3 and f4.
# Here, the GCV formulation corresponds to the ECV method in the paper.

facnum_rate_estc_n50 <- matrix(NA, nrow = 10, ncol = 6)
colnames(facnum_rate_estc_n50) <- c("ES_K5", "ES_K5_pen3", "ES_K5_pen4", 
                                    "GCV_K5", "GCV_K5_pen3", "GCV_K5_pen4")
rownames(facnum_rate_estc_n50) <- paste0("C", 1:10)

for(i in 1:10) {
  file_path <- paste0("~/glmcv/simuV2/estimate_c/n50_p50_q4_C", i, "_frequency.rds")
  data <- readRDS(file_path)

  columns_needed <- c("ES_K5", "ES_K5_pen3", "ES_K5_pen4", 
                      "GCV_K5", "GCV_K5_pen3", "GCV_K5_pen4")
  

  for(j in 1:length(columns_needed)) {

    col_data <- data[, columns_needed[j]]
    correct_rate <- col_data[4] / sum(col_data)
    facnum_rate_estc_n50[i, j] <- correct_rate
  }
}


########### plot 
plot_data <- as.data.frame(facnum_rate_estc_n50)
plot_data$C <- 1:10

data_long <- plot_data %>%
  pivot_longer(cols = -C, names_to = "Method", values_to = "Rate")

data_long <- data_long %>%
  mutate(
    Method_group = case_when(
      Method %in% c("ES_K5", "ES_K5_pen3", "ES_K5_pen4") ~ "ES",
      Method %in% c("GCV_K5", "GCV_K5_pen3", "GCV_K5_pen4") ~ "GCV"
    ),
    Method_type = case_when(
      Method %in% c("ES_K5", "GCV_K5") ~ "base",
      Method %in% c("ES_K5_pen3", "GCV_K5_pen3") ~ "pen1",
      Method %in% c("ES_K5_pen4", "GCV_K5_pen4") ~ "pen2"
    )
  )


data_long$Method_factor <- factor(data_long$Method,
                                  levels = c("ES_K5", "ES_K5_pen3", "ES_K5_pen4",
                                             "GCV_K5", "GCV_K5_pen3", "GCV_K5_pen4"))
p2 <- ggplot(data_long, aes(x = C, y = Rate, color = Method_factor, 
                            shape = Method_factor, group = Method_factor)) +
  geom_line() +                     
  geom_point(size = 3) +           
  geom_vline(xintercept = sqrt(20), linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("ES_K5" = "red", "ES_K5_pen3" = "red", "ES_K5_pen4" = "red",
               "GCV_K5" = "blue", "GCV_K5_pen3" = "blue", "GCV_K5_pen4" = "blue"),
    labels = c("ES_K5" = expression(ES[5]), 
               "ES_K5_pen3" = expression(p[1]*ES[5]), 
               "ES_K5_pen4" = expression(p[2]*ES[5]),
               "GCV_K5" = expression(ECV[5]), 
               "GCV_K5_pen3" = expression(p[1]*ECV[5]), 
               "GCV_K5_pen4" = expression(p[2]*ECV[5]))
  ) +  
  scale_shape_manual(
    values = c("ES_K5" = 0, "ES_K5_pen3" = 2, "ES_K5_pen4" = 1,
               "GCV_K5" = 0, "GCV_K5_pen3" = 2, "GCV_K5_pen4" = 1),
    labels = c("ES_K5" = expression(ES[5]), 
               "ES_K5_pen3" = expression(p[1]*ES[5]), 
               "ES_K5_pen4" = expression(p[2]*ES[5]),
               "GCV_K5" = expression(ECV[5]), 
               "GCV_K5_pen3" = expression(p[1]*ECV[5]), 
               "GCV_K5_pen4" = expression(p[2]*ECV[5]))
  ) + 
  scale_x_continuous(
    breaks = 1:10  
  ) + 
  labs(x = "C",                     
       y = "Correct selection rate") +  
  theme_classic() +                 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = c(1, 0),      
    legend.justification = c(1.005, 0),
    legend.margin = margin(5, 5, 5, 5), 
    legend.box.background = element_rect(colour = "black", fill = "white", linewidth = 0.5),  # 使用box.background
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.title = element_blank(),
    legend.text = element_text(size = 14),      
    legend.key.size = unit(1, "cm"),            
    legend.spacing.x = unit(0.2, "cm"),         
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#print(p2)
ggsave("~/glmcv/simuV2/estimate_c/factor_selection_rate.pdf", p2, width = 9, height = 7, units = "in")

