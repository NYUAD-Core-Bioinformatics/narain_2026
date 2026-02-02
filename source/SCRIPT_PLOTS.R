#Scripts 
###Script for pathway analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

#upload data 
df <- read_csv("file_path.csv")
df$Condition <- factor(df$Condition, levels = c("Resilient", "SI", "SPT", "Susceptible", "TST")) #factor the condition

#set up matrix
all_pathways <- unique(df$Pathways)
all_conditions <- levels(df$Condition)
full_grid <- expand.grid(Pathways = all_pathways, Condition = all_conditions)
df_complete <- left_join(full_grid, df, by = c("Pathways", "Condition"))

#ranking pathways in PValueLog
pathway_order <- df_complete %>%
  group_by(Pathways) %>%
  summarise(avg_logp = mean(PValueLog, na.rm = TRUE)) %>%
  arrange(desc(avg_logp)) %>%
  pull(Pathways)

df_complete$Pathways <- factor(df_complete$Pathways, levels = pathway_order)

# plots
bubble_plot <- ggplot(df_complete, aes(x = Condition, y = Pathways)) +
  geom_point(aes(size = PValueLog, color = Zscore), shape = 16, na.rm = FALSE) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, na.value = "grey60", name = "Z-score"
  ) +
  scale_size_continuous(name = "-Log10(p-value)", range = c(3, 10)) +
  theme_bw(base_size = 12) + 
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    axis.ticks = element_line(size = 0.3),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# saving 
ggsave(
  filename = "file_path.png",
  plot = bubble_plot,
  width = 8, height = 10, dpi = 300, bg = "white"
)

###Script for Bar GRAPHS
library(dplyr)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)

#load data
sample_data <- read.csv("file_path.csv")


sample_data$phenotype <- factor(
  sample_data$phenotype,
  levels = c("Control","Resilient","SI","SPT","TST","Susceptible")
)

#statistics_summary
score_data <- sample_data %>%
  group_by(phenotype) %>%
  summarise(
    si_mean = mean(si_score), si_se = sd(si_score)/sqrt(n()),
    spt_mean = mean(spt_score), spt_se = sd(spt_score)/sqrt(n()),
    tst_mean = mean(tst_score), tst_se = sd(tst_score)/sqrt(n()),
    n = n()
  )

#tukey sig
get_tukey_sig <- function(model, group1, group2)
{
  tuk <- TukeyHSD(model)[[1]]
  row1 <- paste0(group1, "-", group2)
  row2 <- paste0(group2, "-", group1)
  
  if (row1 %in% rownames(tuk)) p <- tuk[row1, "p adj"]
  else if (row2 %in% rownames(tuk)) p <- tuk[row2, "p adj"]
  else return("")
  
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}

#comparisions
comparisons <- list(
  si = list(
    c("Resilient","SI"),
    c("Resilient","Susceptible"),
    c("SI","Susceptible"),
    c("Control","SI"),
    c("Control","Susceptible")
  ),
  spt = list(
    c("Resilient","SPT"),
    c("Resilient","Susceptible"),
    c("SPT","Susceptible"),
    c("SPT","Control"),
    c("Susceptible","Control")
  ),
  tst = list(
    c("Resilient","TST"),
    c("Resilient","Susceptible"),
    c("Control","Susceptible"),
    c("Control","TST"),
    c("Control","SI"),
    c("Control","SPT")
  )
)

#anova+tukey
significance_labels <- list()

for(test in names(comparisons)){
  
  score_col <- paste0(test, "_score")
  
  model <- aov(sample_data[[score_col]] ~ sample_data$phenotype)
  
  stars <- sapply(
    comparisons[[test]],
    function(x) get_tukey_sig(model, x[1], x[2])
  )
  
  significance_labels[[test]] <- stars
}

#plot
ylabs <- c(
  si  = "SI score",
  spt = "Sucrose preference (%)",
  tst = "Total immobility time (s)"
)

tests <- c("si","spt","tst")

for(test in tests){
  
  mean_col <- paste0(test, "_mean")
  se_col   <- paste0(test, "_se")
  
  comp_list <- comparisons[[test]]
  lab_list  <- significance_labels[[test]]
  
  sig_idx <- which(lab_list != "")
  comp_list <- comp_list[sig_idx]
  lab_list  <- lab_list[sig_idx]
  
  n_comp <- length(comp_list)
  
  y_max <- max(score_data[[mean_col]] + score_data[[se_col]], na.rm = TRUE)
  extra <- 20 * n_comp
  y_limit <- y_max + extra + 30
  
  y_positions <- seq(
    from = y_max + 15,
    by = 20,
    length.out = n_comp
  )
  
  p <- ggplot(score_data, aes(x = phenotype, y = .data[[mean_col]], fill = phenotype)) +
    geom_bar(stat="identity", width=0.55, color="black", size=0.7) +
    geom_errorbar(
      aes(
        ymin = .data[[mean_col]] - .data[[se_col]],
        ymax = .data[[mean_col]] + .data[[se_col]]
      ),
      width = 0.22,
      size = 0.6
    ) +
    scale_fill_manual(values = c(
      "Control"     = brewer.pal(9,"Greys")[5],
      "Resilient"   = brewer.pal(9,"Blues")[7],
      "SI"          = brewer.pal(9,"Reds")[7],
      "SPT"         = brewer.pal(9,"Oranges")[5],
      "TST"         = brewer.pal(9,"Greens")[6],
      "Susceptible" = brewer.pal(9,"Purples")[5]
    )) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0,0)) +
    labs(x = "Phenotype", y = ylabs[[test]]) +
    theme(
      panel.background = element_rect(fill="white"),
      plot.background  = element_rect(fill="white"),
      panel.grid = element_blank(),
      axis.line  = element_line(color="black", size=0.6),
      axis.ticks = element_line(color="black", size=0.5),
      axis.text  = element_text(size=12),
      axis.title = element_text(size=14, face="bold"),
      legend.position = "none"
    )
  
  if(n_comp > 0){
    p <- p + geom_signif(
      comparisons = comp_list,
      annotations = lab_list,
      y_position = y_positions,
      tip_length = 0.01,
      size = 0.6,
      textsize = 5
    )
  }
  #saving 
  ggsave(
    paste0("file_path", test, "_tukey.png"),
    p, width = 7, height = 5, dpi = 600
  )
}

###Script for ROC (I am not sure if we used this as we had a mistake in TST but i will attach it here)

#  libraries
library(ComplexHeatmap)
library(ComplexUpset)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(pROC)
library(RColorBrewer)
library(readxl)

#initialize loop variables
lhb_celltype <- c("Oval-Medial", "Marginal", "Lateral", "HbX")
phenotype <- c("resilient", "susceptible", "si", "spt", "tst")
test <- c("si", "spt", "tst")

#plot roc curves
sample_data <- read.csv("file_path.csv")
for (i in 3:length(phenotype)) {
  roc_test <- paste0(phenotype[i], "_score")
  roc_list <- roc(sample_data[, "condition"], sample_data[, roc_test], direction = ifelse(phenotype[i] == "tst", "<", ">"))
  roc_data <- data.frame(tpr = roc_list$sensitivities,
                         fpr = 1 - roc_list$specificities,
                         golden_ratio = roc_list$sensitivities + roc_list$specificities,
                         threshold = roc_list$thresholds)
  roc_cutoffs <- roc_data[order(roc_data$golden_ratio, decreasing = T), c("tpr", "fpr")][1, ]
  roc_plot <- ggroc(roc_list, legacy.axes = T, size = 2, color = brewer.pal(9, "Blues")[5]) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                 linetype = "dashed", color = "black",
                 size = 1) +
    geom_segment(aes(x = roc_cutoffs$fpr, xend = roc_cutoffs$fpr, y = -0.02, yend = roc_cutoffs$tpr),
                 linetype = "solid", color = brewer.pal(9, "Reds")[7],
                 size = 3) +
    geom_segment(aes(x = -0.025, xend = roc_cutoffs$fpr + 0.00525, y = roc_cutoffs$tpr, yend = roc_cutoffs$tpr),
                 linetype = "solid", color = brewer.pal(9, "Reds")[7],
                 size = 3) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.02, 1.01), labels = scales::percent) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.025, 1.01), labels = scales::percent) +
    ylab("True positive rate") +
    xlab("False positive rate") +
    theme(plot.margin = margin(0.4, 1.8, 0.6, 0.2, unit = "cm"),
          panel.background = element_blank(),
          axis.line = element_line(size = 2, color = "black"),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(0.4, "cm"),
          axis.title.y = element_text(face = "bold", size = 34),
          axis.title.x = element_text(face = "bold", size = 34, vjust = -1),
          axis.text.y = element_text(size = 30, color = "black"),
          axis.text.x = element_text(size = 30, color = "black", vjust = -0.2, hjust = 0.3))
  roc_plot
  ggsave(paste0("~/Downloads/", phenotype[i], "_roc.png"), roc_plot, width = 10, height = 8) #saving 
  
}


