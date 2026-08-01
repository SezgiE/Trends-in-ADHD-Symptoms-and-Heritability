## Sezgi Ercan - 20.06.2025
## Script 4 for item level analysis

##----------------------- OpenMx Model - Result Plotting ----------------------##

library(tidyverse)
library(readxl)
library(gridExtra)
library(grid)


rm(list = ls(all=TRUE)) # clean memory

# This script takes OpenMx result files from third script as input.
# Please set the working directory to where the above-mentioned files are.
setwd("./item_level_analyses/model_results_IL")

# Please indicate your desired output directory. Plots will be automatically saved there.
output_dir <- "./item_level_analyses/model_results_IL"

data <- data.frame(read_excel("openmx_results.xlsx"))
trend_df <- data.frame(read_excel("trend_plot_df.xlsx"))


#
#----------------------------------Code Starts---------------------------------#
#


alpha_threshold <- round(0.05/36, 5)

# selecting necessary variables
data <- data  %>%
  filter(comparison == "no_trend") %>%
  mutate(sig = ifelse(p < alpha_threshold, "Significant", "Non-significant")) %>%
  dplyr::select(q_name, rater, age, p, sig)


# Data and label adjustments 
data$age <- factor(data$age, levels=c("7", "10", "12"))
data$respondent <- gsub('m\\b', 'Mother', data$rater)
data$respondent <- gsub('f\\b', 'Father', data$respondent)
data$respondent <- gsub("t\\b", 'Teacher', data$respondent)
data$question_number <- as.numeric(str_extract(data$q_name, "\\d+"))


# ranking function
rank_p_centered <- function(p_values) {
  ranked <- rank(p_values)
  
  middle_rank <- median(ranked)
  
  centered_ranks <- ranked - middle_rank
  
  return(centered_ranks)
}

trend_direct <- trend_df %>%
  mutate(item_key = paste0(q_name, rater, rated_age),
         trend_direction = ifelse(beta_bc > 0, "UP", "DOWN")) %>%
  dplyr::select(item_key, trend_direction)

data <- data %>%
  group_by(rater, age) %>%
  mutate(p_rank = rank_p_centered(p),
         q_rank = rank_p_centered(question_number),
         item_name = gsub("q", "", q_name),
         item_key = paste0(q_name, rater, age)) %>%
  ungroup() %>%
  left_join(distinct(trend_direct, item_key, .keep_all = TRUE), by = "item_key") %>%
  # Create a unified trend variable based on significance and direction
  mutate(trend_status = case_when(
    sig == "Non-significant" ~ "No Trend",
    trend_direction == "UP" ~ "Upward",
    trend_direction == "DOWN" ~ "Downward",
    TRUE ~ "No Trend"
  )) %>%
  mutate(trend_status = factor(trend_status, levels = c("Upward", "Downward", "No Trend")))

# Color Palette Hex Codes
npg_colors <- c("Upward" = "#E64B35", "Downward" = "#4DBBD5", "No Trend" = "#CCCCCC")

my_plot <- data %>%
  ggplot(aes(age, q_rank)) +
  # Map color to the new trend_status and reduce point size for print
  geom_point(aes(color = trend_status), size = 8.5) + 
  geom_text(aes(label = ifelse(sig == "Significant", 
                               paste0(item_name, "\n", ifelse(trend_direction == "UP", "▲", "▼")), 
                               item_name)), 
            size = 2.5, lineheight = 0.8, color = "black") + # Text smaller and strictly black
  facet_wrap(~respondent, strip.position = "bottom") +
  scale_x_discrete(expand = c(.5, .5)) +
  theme_minimal()  +
  scale_color_manual(values = npg_colors) + 
  theme(
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(-0.5, "cm"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -30, r = 0, b = 0, l = 0),
    axis.line.x = element_line(color = "black", linewidth = 0.5), # Thinner lines
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(0.1, "cm"),
    strip.text = element_text(size = 10, face = "bold"), 
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", margin = margin(b = 15, t = 15))) +
  labs(color = "Trend Direction", 
       x = "", 
       y = "CBCL/6-18 Items") +
  guides(color = guide_legend(override.aes = list(size = 3))) # Smaller legend points


my_plot_path <- file.path(output_dir, "model_res_plot.pdf")
ggsave(my_plot_path,
       plot = my_plot, 
       width = 8, height = 4.5, units = "in",
       bg = "white", device = cairo_pdf)


# Item level trend predicted by the model
trend_df <- trend_df %>%
  mutate(instensity = ifelse(significance, 1, 0.75)) %>%
  filter(significance == TRUE)

items <- unique(trend_df$q_name)

p_list <- list()

for (item in items){
  

  plot_df <- trend_df %>%
    filter(q_name == item) %>%
    mutate(rated_age = as.numeric(as.character(rated_age)))
  
  ## plotting the trend
  p <- ggplot(plot_df, aes(x = cohort, y = mean_liability, color = rater)) +
    geom_line(aes(linetype = as.factor(rated_age), group = interaction(rater, rated_age)), size = 1) +
    labs(
      title = toupper(item),
      x = "Birth Cohort",
      y = "Mean Liability",
      color = "Respondent",
      linetype = "Age"
    ) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"),
      plot.title = element_text(hjust = 0.5, size = 25, family = "Arial"),
      axis.title = element_text(size = 22, family = "Arial"),
      axis.text = element_text(size = 20, family = "Arial"),
      legend.title = element_text(size = 20, family = "Arial"),
      legend.text = element_text(size = 20, family = "Arial"),
      legend.position = "bottom",  # legend to the bottom
      legend.box = "horizontal",  # legend items horizontally
      panel.grid.major = element_line(color = "gray90", size = 0.2),  # light grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      strip.background = element_blank(),  # Remove background for faceted plots
      strip.text = element_text(size = 15, family = "Arial")  # Customize facet labels
    ) +
    scale_color_manual(values = c("m" = "#1f77b4", "f" = "#ff7f0e", "t" = "#2ca02c")) + 
    scale_linetype_manual(values = c("7" = "solid", "10" = "dashed", "12" = "dotted")) +
    scale_x_continuous(breaks = c(1990, 2000, 2010))
    guides(linetype = guide_legend(override.aes = list(color = "black")))

  
  p
  
  p_list[[item]] <- p

}

plot_names <- names(p_list)

p_list_sorted <- p_list[order(as.numeric(gsub("q", "", plot_names)))]

# Trends in Mean Liability Over 30 Years by Birth Cohort, Rater and Age
trend_plot <- grid.arrange(grobs = p_list_sorted, ncol = 4, nrow = 3)


plot_path <- file.path(output_dir, "model_trend_plot.pdf")

ggsave(plot_path,
       plot = trend_plot, width = 32, height = 18, units = "in", bg = "white", device = cairo_pdf)










