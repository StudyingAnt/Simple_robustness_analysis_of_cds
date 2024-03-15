library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(egg)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

gc_bias_file <- paste(base_path, "SBS_roughness_gc_bias.csv", sep="")
gc_bias_data <- read.csv(gc_bias_file)

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


means_data <- data_mean_sorted_desc %>% 
  arrange(factor(signature, levels=signatures_noseqerr))
stds_data <- data_std_sorted_desc %>% 
  arrange(factor(signature, levels=signatures_noseqerr))
entropies_data <- gc_bias_data %>% 
  filter(Signature %in% signatures_noseqerr) %>% 
  arrange(factor(Signature, levels=signatures_noseqerr))

mean_vs_entropy <- data.frame(
  Signature = factor(signatures_noseqerr, levels=signatures_noseqerr),
  entropy = entropies_data$SBS_entropy,
  gc_bias_tr = entropies_data$SBS_gc_bias_tr,
  mean_robustness = means_data$avg,
  std_robustness = stds_data$std
)

rr_plot <- ggplot(mean_vs_entropy, aes(x=entropy, y=mean_robustness, color=gc_bias_tr)) +
  geom_point(size = 0.35) +
  scale_color_gradient2(midpoint=0, 
                        low="blue", 
                        mid="#f2cc8f",
                        high="red") +
  xlab("Entropy") +
  ylab("Robustness") +
  labs(color="GC bias") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 7, family = "Arial"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "entropy_vs_robustness.png", sep="")
ggsave(out_file, plot = rr_plot, dpi=1200, width = 60, height = 50, units = "mm")

rr_plot <- ggplot(mean_vs_entropy, aes(x=entropy, y=std_robustness, color=gc_bias_tr)) +
  geom_point(size = 0.35) +
  scale_color_gradient2(midpoint=0, 
                        low="blue", 
                        mid="#f2cc8f",
                        high="red") +
  xlab("Entropy") +
  ylab("SD of robustness") +
  labs(color="GC bias") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 7, family = "Arial"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "entropy_vs_robustness_std.png", sep="")
ggsave(out_file, plot = rr_plot, dpi=1200, width = 60, height = 50, units = "mm")


