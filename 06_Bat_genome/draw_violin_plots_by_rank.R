# import package
library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(egg)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read all data
all_sig_data <- list()
for (species in c("all_gencode", "molMol2", "myoMyo6")) {
    rum_file <- paste(base_path, species, "_rums_profile.csv", sep="")
    all_sig_data[[species]] <- read.csv(rum_file)
}

#
signatures = c(
  "SBS1", "SBS2", "SBS3", "SBS4", "SBS6", "SBS9", "SBS13", "SBS18"
)

sbss <- c()
speciess <- c()
data_flatten <- c()
for (species in c("all_gencode", "molMol2", "myoMyo6")) {
  # all_sig_data[[species]]
  n <- length(all_sig_data[[species]]$Transcript)
  
  sbss <- c(sbss, rep(signatures, each=n))
  
  for (i in 1:length(signatures)) {
    data_flatten <- c(data_flatten, all_sig_data[[species]][, signatures[i]])
  }
  
  if (species == "all_gencode") {
    speciess <- c(speciess, rep("human",n*length(signatures)) )
  } else {
    speciess <- c(speciess, rep(species,n*length(signatures)) )
  }
}

all_sbs_data_plot <- data.frame(
  Species = speciess,
  Signature = sbss,
  Robustness = data_flatten
)

# M. molossus, M. myotis, H. sapiens
all_sbs_data_plot$Species <- factor(all_sbs_data_plot$Species, 
                                    levels = c("human", "molmol2", "myomyo2"), 
                                    labels = c("H. sapiens", "M. molossus", "M. myotis"))


for (sig in c("SBS1", "SBS4", "SBS6", "SBS18")) {
  big.plot <- ggplot(all_sbs_data_plot %>% filter(Signature == sig), aes(x=Species, y=Robustness, fill = Species)) + 
    geom_violin(width = 0.8, linewidth = 0.2) +
    geom_boxplot(width = 0.3, size = 0.2, outlier.size = 0.3, outlier.stroke = 0) +
    ylim(0,3) +
    scale_x_discrete(labels = c(expression(italic("H. sapiens")),
                                expression(italic("M. molossus")),
                                expression(italic("M. myotis"))
    )) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )
  
  ggsave(paste(base_path, sig, ".pdf", sep = ""), plot = big.plot, width = 30, height = 60, units = "mm")
  
}

for (sig in c("SBS2", "SBS13")) {
  big.plot <- ggplot(all_sbs_data_plot %>% filter(Signature == sig), aes(x=Species, y=Robustness, fill = Species)) + 
    geom_violin(width = 0.8, linewidth = 0.2) +
    geom_boxplot(width = 0.3, size = 0.2, outlier.size = 0.3, outlier.stroke = 0) +
    ylim(0,35) +
    scale_x_discrete(labels = c(expression(italic("H. sapiens")),
                                expression(italic("M. molossus")),
                                expression(italic("M. myotis"))
    )) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2.5, color = "red", fill = NA, size = 0.5)
  
  small.plot <- ggplot(all_sbs_data_plot %>% filter(Signature == sig), aes(x=Species, y=Robustness, fill = Species)) + 
    geom_violin(width = 0.8, linewidth = 0.2) +
    geom_boxplot(width = 0.3, size = 0.2, outlier.size = 0.3, outlier.stroke = 0) +
    ylim(0,2.5) +
    #ggtitle(bquote(.(sigs[i]) ~ "("*italic(.(spes[i]))*")")) +
    scale_x_discrete(labels = c(expression(italic("H. sapiens")),
                                expression(italic("M. molossus")),
                                expression(italic("M. myotis"))
    )) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),  # Remove x-axis title
      axis.text.x = element_blank(),   # Remove x-axis text
      axis.text.y = element_text(size = 5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "none",
      plot.background = element_rect(colour = "red", fill = NA, size = 0.5)
    )
  ggsave(paste(base_path, sig, ".pdf", sep = ""), plot = big.plot, width = 30, height = 60, units = "mm")
  ggsave(paste(base_path, sig, "_smallRange_v2.pdf", sep = ""), plot = small.plot, width = 25, height = 47, units = "mm")
}









sbs1.plot <- ggplot(data_plot %>% filter(signature == "SBS1"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS1 Spontaneous deamination") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs1.comp.png", sbs1.plot, dpi = 600, width = 80, height = 35, units = "mm")

sbs2.plot <- ggplot(data_plot %>% filter(signature == "SBS2"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS2 APOBEC activity") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs2.comp.png", sbs2.plot, dpi = 600, width = 80, height = 35, units = "mm")


# Calculate quartiles
result <- data_plot %>%
  group_by(species, signature) %>%
  summarise(
    Q1_robustness = quantile(robustness, 0.25),  # Calculating the third quartile
    Q3_robustness = quantile(robustness, 0.75),  # Calculating the third quartile
    IQR_robustness = IQR(robustness),
    .groups = 'drop'  # This drops the grouping structure after summarizing
  ) %>% 
  mutate(lower_bound = Q1_robustness - 1.5*IQR_robustness, upper_bound = Q3_robustness +1.5*IQR_robustness)

result




sbs2.no.plot <- ggplot(data_plot %>% filter(signature == "SBS2"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS2 APOBEC activity") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs2.comp.noout.png", sbs2.no.plot, dpi = 600, width = 80, height = 35, units = "mm")


sbs4.plot <- ggplot(data_plot %>% filter(signature == "SBS4"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS4 Tobacco smoking") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs4.comp.png", sbs4.plot, dpi = 600, width = 80, height = 35, units = "mm")

sbs6.plot <- ggplot(data_plot %>% filter(signature == "SBS6"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS6 Defective DNA mismatch repair") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs6.comp.png", sbs6.plot, dpi = 600, width = 80, height = 35, units = "mm")

sbs13.plot <- ggplot(data_plot %>% filter(signature == "SBS13"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS13 APOBEC activity") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs13.comp.png", sbs13.plot, dpi = 600, width = 80, height = 35, units = "mm")

sbs18.plot <- ggplot(data_plot %>% filter(signature == "SBS18"), aes(x = species, y = robustness, fill = species)) + 
  geom_violin(width = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.3, size = 0.5, outlier.size = 0.2) +
  ggtitle("SBS18 Damage by ROS") + 
  coord_flip() +
  ylab("Robustness") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),              
    axis.title = element_text(size = 7),              
    axis.text = element_text(size = 6),
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none", 
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
  )

ggsave("sbs18.comp.png", sbs18.plot, dpi = 600, width = 80, height = 35, units = "mm")



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



# signatures = c(
#   "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
#   "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
#   "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
#   "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
#   "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
#   "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
#   "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
#   "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
#   "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
#   "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
#   "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
#   "SBS95"
# )
# 
# signatures_noseqerr = c(
#   "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", 
#   "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d", 
#   
#   "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
#   "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", 
#   
#   "SBS15", "SBS16", "SBS17a", "SBS17b", "SBS18", 
#   "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", 
#   
#   "SBS24", "SBS25", "SBS26",  "SBS28", "SBS29", 
#   "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", 
#   
#   "SBS35", "SBS36", "SBS37", "SBS38", "SBS39", 
#   "SBS40", "SBS41", "SBS42", "SBS44", "SBS84", 
#   
#   "SBS85", "SBS86", "SBS87", "SBS88", "SBS89", 
#   "SBS90", "SBS91", "SBS92", "SBS93", "SBS94"
# )
# 
# signatures_aetiology = c(
#   "SBS1", "SBS2", "SBS13", "SBS84", "SBS85",     # Deamination activity
#   "SBS4", "SBS29", "SBS92",                      # Tobacco exposure
#   "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38",   # UV exposure
#   "SBS18",                                       # ROS exposure
#   "SBS22", "SBS24", "SBS42", "SBS88", "SBS90",   # Chemical exposure
#   "SBS11", "SBS25", "SBS31", "SBS32", "SBS35",   # Anti-cancer treatment
#   "SBS86", "SBS87",
#   "SBS3", "SBS6", "SBS9", "SBS10a", "SBS10b",    # Defective DNA repair and
#   "SBS10c", "SBS10d", "SBS14", "SBS15", "SBS20", # hypermutations 
#   "SBS21", "SBS26", "SBS30", "SBS36", "SBS44"
# )
