# import package
library(tidyverse)
library(dplyr)
library(ggplot2)
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

# 
length_file <- paste(base_path, "gencode.v40.pc_transcripts.nopary.cdsplus.length.csv", sep = "")
length_data <- read.csv(length_file)

data2 <- data %>% arrange(Transcript)
length_data2 <- length_data %>% arrange(transcript_name)

new.df <- data.frame(
  rums = data2$SBS2,
  transcript_length = length_data2$length
)

plot.lst <- list()
for (sig in signatures) {
  new.df <- data.frame(
    rums = data2[[sig]],
    transcript_length = length_data2$length
  )
  
  lenrums.plot <- ggplot(new.df, aes(x=log10(transcript_length), y=rums)) + 
    geom_point(size = 0.001)+  # Adding a title for each plot
    theme(
      axis.text = element_blank(),  # Remove axis text/labels
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.title = element_blank()   # Remove axis titles
    )
  
  plot.lst[[sig]] <- lenrums.plot
  
  #ggsave(paste(base_path, "rums_length_", sig, ".png", sep=""), plot = lenrums.plot, dpi = 600, width = 100, height = 100, unit = "mm")
  
}
library(gridExtra)

all.plot.grid <- grid.arrange(grobs = plot.lst, ncol = 8, left = "Robustness", bottom = "Log10(length)")

ggsave(paste(base_path, "rumslength.all.png", sep=""), plot = all.plot.grid, dpi = 600, width = 180, height = 200, unit = "mm")


test.rlt <- cor.test(data2$SBS1, length_data2$length, method = "spearman")

rhos <- c()
p.vals <- c()
for (sig in signatures) {
  test.rlt <- cor.test(data2[[sig]], length_data2$length, method = "spearman")
  
  rhos <- c(rhos, test.rlt$estimate)
  p.vals <- c(p.vals, test.rlt$p.value)
  
}

test.df <- data.frame(
  signature = factor(signatures, levels = signatures),
  rho = rhos,
  p.val = -log10(p.vals)
)


lenrums.rho.plot <- ggplot(test.df, aes(x=signature, y=rho, size = p.val, color = p.val)) +
  geom_point() +
  ylim(-1, 1) + xlab("Signature") + ylab("Spearman's rho") +
  scale_size(name = "-Log10(P)",range = c(0.5, 3)) +
  scale_color_gradient(name = "-Log10(P)", low = "black", high = "red") +  # Name the color scale
  guides(
    color = guide_legend(),  # Use default legend for color
    size = guide_legend()  # Use default legend for size
  )+  # Name and colors for gradient
  theme(
    axis.title = element_text(size = 8),  # Font size for axis titles (x and y)
    axis.text = element_text(size = 6),   # Font size for axis text (x and y)
    legend.title = element_text(size = 8),  # Font size for legend title
    legend.text = element_text(size = 6),  # Font size for legend items
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)  # Specific for x-axis text with rotation
  )

ggsave(paste(base_path, "rumslength.rho.png", sep=""), plot = lenrums.rho.plot, dpi = 600, width = 180, height = 200, unit = "mm")
