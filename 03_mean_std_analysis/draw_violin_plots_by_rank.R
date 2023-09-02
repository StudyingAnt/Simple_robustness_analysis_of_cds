# import package
library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(egg)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read file
input_file <- paste(base_path, "all_gencode_rums_profile.csv", sep="")

data <- read.csv(input_file)
n_transcript <- length(data$Transcript)

signatures = c(
  "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
  "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
  "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
  "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
  "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
  "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
  "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
  "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
  "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
  "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
  "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
  "SBS95"
)

signatures_noseqerr = c(
  "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", 
  "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d", 
  
  "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
  "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", 
  
  "SBS15", "SBS16", "SBS17a", "SBS17b", "SBS18", 
  "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", 
  
  "SBS24", "SBS25", "SBS26",  "SBS28", "SBS29", 
  "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", 
  
  "SBS35", "SBS36", "SBS37", "SBS38", "SBS39", 
  "SBS40", "SBS41", "SBS42", "SBS44", "SBS84", 
  
  "SBS85", "SBS86", "SBS87", "SBS88", "SBS89", 
  "SBS90", "SBS91", "SBS92", "SBS93", "SBS94"
)

signatures_aetiology = c(
  "SBS1", "SBS2", "SBS13", "SBS84", "SBS85",     # Deamination activity
  "SBS4", "SBS29", "SBS92",                      # Tobacco exposure
  "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38",   # UV exposure
  "SBS18",                                       # ROS exposure
  "SBS22", "SBS24", "SBS42", "SBS88", "SBS90",   # Chemical exposure
  "SBS11", "SBS25", "SBS31", "SBS32", "SBS35",   # Anti-cancer treatment
  "SBS86", "SBS87",
  "SBS3", "SBS6", "SBS9", "SBS10a", "SBS10b",    # Defective DNA repair and
  "SBS10c", "SBS10d", "SBS14", "SBS15", "SBS20", # hypermutations 
  "SBS21", "SBS26", "SBS30", "SBS36", "SBS44"
)

# change data format
sbs <- rep(signatures, each=n_transcript)

data_flatten <- c()
for (i in 1:79) {
  data_flatten <- c(data_flatten, data[, signatures[i]])
  }

data_plot <- data.frame(
  signature = sbs,
  robustness = data_flatten
  )

# FOR ALL SIGNATURES
data_mean_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(avg=mean(robustness)) %>% 
  arrange(desc(avg))

data_std_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(std=sd(robustness)) %>% 
  arrange(desc(std))

# all signatures
idx_starts = c(1, 11, 21, 31, 41, 51, 61, 71)
idx_ends = c(10, 20, 30, 40, 50, 60, 70, 79)

for (i in 1:length(idx_starts)){
  sig_mean <- data_mean_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  sig_std <- data_std_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  
  data_mean <- data_plot %>% filter(signature %in% sig_mean)
  data_std <- data_plot %>% filter(signature %in% sig_std)
  
  violin.plot.mean <- ggplot(data_mean, 
                             aes(x = fct_rev(factor(signature, levels=sig_mean)), 
                                 y = robustness, 
                                 fill = fct_rev(factor(signature, levels=sig_mean)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by average") +
    xlab("Signature") +
    ylab("Robustness") +
    #ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    ) 
  
  violin.plot.mean.fixed <- set_panel_size(violin.plot.mean,
                                           width = unit(40, "mm"),
                                           height = unit(30, "mm"))
  
  out_file_mean <- paste(base_path, "all_signature_violin_rank_by_mean_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_mean, plot = violin.plot.mean.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
  violin.plot.std <- ggplot(data_std, 
                         aes(x = fct_rev(factor(signature, levels=sig_std)), 
                             y = robustness, 
                             fill = fct_rev(factor(signature, levels=sig_std)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by standard deviation") +
    xlab("Signature") +
    ylab("Robustness") +
    #ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    )
  
  violin.plot.std.fixed <- set_panel_size(violin.plot.std,
                                       width = unit(40, "mm"),
                                       height = unit(30, "mm"))
  
  out_file_std <- paste(base_path, "all_signature_violin_rank_by_std_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_std, plot = violin.plot.std.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
}


# REMOVE SEQUENCING ERROR
data_mean_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(avg=mean(robustness)) %>% 
  arrange(desc(avg)) %>% 
  filter(signature %in% signatures_noseqerr)

data_std_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(std=sd(robustness)) %>% 
  arrange(desc(std)) %>% 
  filter(signature %in% signatures_noseqerr)

# removed sequencing errors
idx_starts = c(1, 11, 21, 31, 41, 51)
idx_ends = c(10, 20, 30, 40, 50, 60)

for (i in 1:length(idx_starts)){
  sig_mean <- data_mean_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  sig_std <- data_std_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  
  data_mean <- data_plot %>% filter(signature %in% sig_mean)
  data_std <- data_plot %>% filter(signature %in% sig_std)
  
  violin.plot.mean <- ggplot(data_mean, 
                             aes(x = fct_rev(factor(signature, levels=sig_mean)), 
                                 y = robustness, 
                                 fill = fct_rev(factor(signature, levels=sig_mean)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by average") +
    xlab("Signature") +
    ylab("Robustness") +
    ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    ) 
  
  violin.plot.mean.fixed <- set_panel_size(violin.plot.mean,
                                           width = unit(40, "mm"),
                                           height = unit(30, "mm"))
  
  out_file_mean <- paste(base_path, "noseqerr_signature_violin_rank_by_mean_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_mean, plot = violin.plot.mean.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
  violin.plot.std <- ggplot(data_std, 
                            aes(x = fct_rev(factor(signature, levels=sig_std)), 
                                y = robustness, 
                                fill = fct_rev(factor(signature, levels=sig_std)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by standard deviation") +
    xlab("Signature") +
    ylab("Robustness") +
    ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    )
  
  violin.plot.std.fixed <- set_panel_size(violin.plot.std,
                                          width = unit(40, "mm"),
                                          height = unit(30, "mm"))
  
  out_file_std <- paste(base_path, "noseqerr_signature_violin_rank_by_std_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_std, plot = violin.plot.std.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
}


# SIGNATURES WITH AETIOLOGY
data_mean_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(avg=mean(robustness)) %>% 
  arrange(desc(avg)) %>% 
  filter(signature %in% signatures_aetiology)

data_std_sorted_desc <- data_plot %>% 
  group_by(signature) %>% 
  summarise(std=sd(robustness)) %>% 
  arrange(desc(std)) %>% 
  filter(signature %in% signatures_aetiology)

# signatures with aetiology
idx_starts = c(1, 11, 21, 31, 32)
idx_ends = c(10, 20, 30, 40, 41)

for (i in 1:length(idx_starts)){
  sig_mean <- data_mean_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  sig_std <- data_std_sorted_desc$signature[idx_starts[i]:idx_ends[i]]
  
  data_mean <- data_plot %>% filter(signature %in% sig_mean)
  data_std <- data_plot %>% filter(signature %in% sig_std)
  
  violin.plot.mean <- ggplot(data_mean, 
                             aes(x = fct_rev(factor(signature, levels=sig_mean)), 
                                 y = robustness, 
                                 fill = fct_rev(factor(signature, levels=sig_mean)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by average") +
    xlab("Signature") +
    ylab("Robustness") +
    ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    ) 
  
  violin.plot.mean.fixed <- set_panel_size(violin.plot.mean,
                                           width = unit(40, "mm"),
                                           height = unit(30, "mm"))
  
  out_file_mean <- paste(base_path, "aetiology_signature_violin_rank_by_mean_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_mean, plot = violin.plot.mean.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
  violin.plot.std <- ggplot(data_std, 
                            aes(x = fct_rev(factor(signature, levels=sig_std)), 
                                y = robustness, 
                                fill = fct_rev(factor(signature, levels=sig_std)))) +
    geom_violin(scale = "width", lwd = 0.1, alpha = 0.5) + 
    geom_jitter(height = 0, width = 0.01, size = 1, shape = ".", stroke = 0.01) +
    ggtitle("Rank by standard deviation") +
    xlab("Signature") +
    ylab("Robustness") +
    ylim(0, 6.5) +
    coord_flip() + 
    theme_light() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 7, family = "Arial"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(
      values = c("#FF61AB", "#FF6176", "#FF8161", "#FFB561", "#FFEA62",
                 "#DFFF61", "#ABFF61", "#76FF61", "#61FF81", "#61FFB5")
    )
  
  violin.plot.std.fixed <- set_panel_size(violin.plot.std,
                                          width = unit(40, "mm"),
                                          height = unit(30, "mm"))
  
  out_file_std <- paste(base_path, "aetiology_signature_violin_rank_by_std_", as.character(idx_starts[i]), "_", as.character(idx_ends[i]), ".png", sep="")
  ggsave(out_file_std, plot = violin.plot.std.fixed, dpi = 1200, width = 60, height = 42.5, units = "mm")
  
}
