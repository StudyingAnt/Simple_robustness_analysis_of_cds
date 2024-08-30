library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(egg)

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

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read file
pc_matrix_file <- paste(base_path, "all_gencode_noseqerr_signatures_PCA_matrix.csv", sep="")
pc_matrix <- read.csv(pc_matrix_file)
gc_bias_file <- paste(base_path, "SBS_roughness_gc_bias.csv", sep="")
gc_bias <- read.csv(gc_bias_file)

gc_bias_noseqerr <- gc_bias %>% arrange(SBS_gc_bias) %>% filter(Signature %in% signatures_noseqerr)

pcs <- c()
for (i in 1:length(pc_matrix$Principal_component)){pcs <- c(pcs, paste("PC", as.character(i), sep=""))}

pcm_heatmap <- ggplot(pc_matrix, aes(x = factor(Principal_component, levels=pcs), y = factor(Signature, levels = gc_bias_noseqerr$Signature),fill=Coefficient)) +
  geom_tile() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                       high="red") +
  xlab("Principal Component") +
  ylab("Signature sorted by GC targeting prefernce") +
  labs(fill="Weight")  +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
  )

pcm_heatmap.fixed <- set_panel_size(pcm_heatmap,
                                    width = unit(140, "mm"),
                                    height = unit(140, "mm"))

out_file <- paste(base_path, "pc_matrix_noseqerr_v2.png", sep="")
ggsave(out_file, plot = pcm_heatmap.fixed, dpi = 600, width = 170, height = 155, units = "mm")


