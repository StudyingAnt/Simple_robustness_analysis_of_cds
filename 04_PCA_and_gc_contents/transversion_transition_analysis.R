library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(egg)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

trv_trs_file <- paste(base_path, "SBS_transversion_transition.csv", sep="")
trv_trs_data <- read.csv(trv_trs_file)

gc_bias_file <- paste(base_path, "SBS_roughness_gc_bias.csv", sep="")
gc_bias_data <- read.csv(gc_bias_file)

entropies_data <- gc_bias_data %>% 
  filter(Signature %in% signatures_noseqerr) %>% 
  arrange(factor(Signature, levels=signatures_noseqerr))

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
trv_trs_sorted_data <- trv_trs_data %>% 
  filter(Signature %in% signatures_noseqerr) %>% 
  arrange(factor(Signature, levels=signatures_noseqerr))

mean_vs_trv_trs <- data.frame(
  Signature = factor(signatures_noseqerr, levels=signatures_noseqerr),
  trv = trv_trs_sorted_data$SBS_transversion,
  trv_tr = trv_trs_sorted_data$SBS_trv_tr,
  trs = trv_trs_sorted_data$SBS_transition,
  trs_tr = trv_trs_sorted_data$SBS_trs_tr,
  mean_robustness = means_data$avg,
  std_robustness = stds_data$std,
  entropy = entropies_data$SBS_entropy
)

## 1
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=entropy, y=mean_robustness, color=trs)) +
  geom_point(size = 0.35) +
  scale_color_gradient2(midpoint=0.4, 
                        low="#4b0092", 
                        mid="#f2cc8f",
                        high="#117733") +
  xlab("Entropy") +
  ylab("Robustness") +
  labs(color="Transition ratio") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "entropy_vs_robustness_trs.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 68.6, height = 50, units = "mm")

## 2
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=entropy, y=std_robustness, color=trs)) +
  geom_point(size = 0.35) +
  scale_color_gradient2(midpoint=0.4, 
                        low="#4b0092", 
                        mid="#f2cc8f",
                        high="#117733") +
  xlab("Entropy") +
  ylab("SD of robustness") +
  labs(color="Transition\n ratio") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "entropy_vs_robustness_trs_std.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 62.1, height = 50, units = "mm")

## 3
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs, y=mean_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transition") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 4
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs, y=std_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transition") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_std.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

################################################################################
############ with line #########################################################
################################################################################
## 1
m <- lm(mean_robustness ~ trv, mean_vs_trv_trs)

eq <- substitute(R^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
lb <- as.character(as.expression(eq))

rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv, y=mean_robustness)) +
  geom_smooth(method = lm, linewidth=0.5, color="red") +
  geom_point(size = 0.35) +
  #geom_text(x=0.83, y=0.78, label= lb, parse = TRUE, size =7/.pt) +
  xlab("Transversion ratio") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 2
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv, y=std_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transversion") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_std_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 3
m <- lm(mean_robustness ~ trs, mean_vs_trv_trs)

eq <- substitute(R^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
lb <- as.character(as.expression(eq))

rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs, y=mean_robustness)) +
  geom_smooth(method = lm, linewidth=0.5, color="red") +
  geom_point(size = 0.35) +
  #geom_text(x=0.17, y=0.78, label= lb, parse = TRUE, size =7/.pt) +
  xlab("Transition ratio") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 4
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs, y=std_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transition") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_std_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

################################################################################
############ transformed #######################################################
################################################################################
## 1
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv_tr, y=mean_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transversion") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_tr.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 2
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv_tr, y=std_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transversion") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_std_tr.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 3
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs_tr, y=mean_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transition") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_tr.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 4
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs_tr, y=std_robustness)) +
  geom_point(size = 0.35) +
  xlab("Transition") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_std_tr.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

################################################################################
############ transformed with line #############################################
################################################################################
## 1
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv_tr, y=mean_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transversion") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_tr_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 2
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trv_tr, y=std_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transversion") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transversion_vs_robustness_std_tr_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 3
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs_tr, y=mean_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transition") +
  ylab("Robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size =5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_tr_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

## 4
rr_plot <- ggplot(mean_vs_trv_trs, aes(x=trs_tr, y=std_robustness)) +
  geom_point(size = 0.35) +
  geom_smooth(method = lm) +
  xlab("Transition") +
  ylab("SD of robustness") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size =5),
    axis.title = element_text(size =7),
    legend.title = element_text(size =7),
    legend.text = element_text(size =5),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

out_file <- paste(base_path, "transition_vs_robustness_std_tr_line.pdf", sep="")
ggsave(out_file, plot = rr_plot, width = 60, height = 50, units = "mm")

