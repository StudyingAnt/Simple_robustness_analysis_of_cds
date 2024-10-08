library(tidyverse)
library(dplyr)
library(extrafont)
library(ggplot2)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read file
pca_shoulder_file <- paste(base_path, "all_gencode_noseqerr_signatures_PCA_shoulder.csv", sep="")
pca_shoulder <- read.csv(pca_shoulder_file)

pca_shoulder

shoulder_plot <- ggplot(pca_shoulder, aes(x=N_components, y=Cumulative_explained_variance*100)) +
  geom_point(size = 0.35) + #, stroke = 0.01) +
  geom_line(linewidth = 0.25) +
  scale_x_continuous(breaks = seq(10, 60, by = 10)) +
  xlab("Number of principal components") +
  ylab("Cumulative explained variance (%)") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7)
  )


shoulder_plot

out_file <- paste(base_path, "all_gencode_noseqerr_signatures_PCA_shoulder_v2.png", sep="")
ggsave(out_file, plot = shoulder_plot, dpi=600, width = 85, height = 60, units = "mm")



