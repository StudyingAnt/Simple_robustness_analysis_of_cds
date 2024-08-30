library(tidyverse)
library(dplyr)
library(ggplot2)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read file
pca_transformed_file <- paste(base_path, "all_gencode_noseqerr_signatures_PCA_transformed.csv", sep="")
pca_transfomred <- read.csv(pca_transformed_file)
gc_contents_file <- paste(base_path, "gencode.v40.pc_transcripts.nopary.cdsplus.gc_contents.csv", sep="")
gc_contents <- read.csv(gc_contents_file)

pca_transfomred_sorted <- pca_transfomred %>% arrange(Transcript) 
gc_contents_sorted <- gc_contents %>% arrange(Transcript_name)

data <- cbind(pca_transfomred_sorted, gc_contents_sorted)

#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

for (i in 1:4){
  for (j in (i+1):5){
    if (i != j) {
      pc_x <- paste("PC", as.character(i), sep = "")
      pc_y <- paste("PC", as.character(j), sep = "")
      
      pca_plot <- ggplot(data, aes(x=get(pc_x), y=get(pc_y), color=GC_wobble)) +
        geom_point(size = 0.0001, shape=20) +
        scale_color_gradient(low = "#e3f2fd", high = "#072ac8", limits = c(0, 1)) +
        xlab(pc_x) +
        ylab(pc_y) +
        labs(color="GC Wobble")  +
        theme_light() +
        theme(
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.height = unit(3, "mm"),
          legend.key.width = unit(3, "mm"),
          legend.box.margin = margin(-10,-10,-10,-10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
        )
      
      out_file <- paste(base_path, "all_gencode_noseqerr_signatures_PCA_", pc_x, "_", pc_y, "_GC_wobble_v2.png", sep="")
      ggsave(out_file, plot = pca_plot, dpi=600, width = 85, height = 60, units = "mm")
    }
  }
}

