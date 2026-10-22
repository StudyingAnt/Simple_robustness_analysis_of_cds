library(tidyverse)
library(ggplot2)
library(lsa)

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

sim_dir <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/validation_out/run_1/"
ori_dir <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/validation_out/"
out_dir <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

cos_sims <- c()
for (signature in signatures){
  sim_file <- paste(sim_dir, tolower(signature),"/", signature, "_length_2200_gc_0.4_sim_1.csv", sep="")
  ori_file <- paste(ori_dir, tolower(signature),"/", signature, ".csv", sep="")
  
  data1 <- read.csv(sim_file)
  data2 <- read.csv(ori_file)
  
  colnames(data1) <- c("mut_type", "percentage")
  value_from <- rep("simulation", 96)
  data1$value_from <- value_from
  
  value_from <- rep("cosmic", 96)
  data2$value_from <- value_from
  
  cos_sim <- round(cosine(data1$percentage, data2$percentage),3)
  cos_sims <- append(cos_sims, cos_sim)
  
  mut_type_level = data1$mut_type
  
  data.combined <- bind_rows(data1, data2)
  
  y_max <- max(data.combined$percentage)
  
  
  p <- ggplot(data.combined, aes(x=factor(mut_type, levels = mut_type_level), y=percentage, fill=value_from)) +
    geom_bar(stat="identity", position="identity", alpha=0.3) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      plot.title = element_text(hjust = 0.5)) +
    ggtitle(signature) +
    xlab("Trinucleotide") + ylab("Percentage") +
    scale_fill_manual(labels = c("COSMIC", "Simulation"), values = c("blue", "red")) +
    guides(fill=guide_legend(title="Value from")) +
    annotate("text", x = 8, y = y_max, label = paste("Cosine similarity:", as.character(cos_sim)))
  
  out_file <- paste(out_dir, signature, "_validation.png", sep="")
  #ggsave(out_file, plot = p, width = 30, height = 8, units = "cm")
}

df.cos_sims <- data.frame(
  signature = signatures,
  cosine_similarity = cos_sims
)

p1 <- ggplot(df.cos_sims, aes(x=factor(signature, levels=signatures), y=cosine_similarity)) +
  geom_point() + ylim(0, 1) +
  ggtitle("Signature validation") +
  xlab("Signature") + ylab("Cosine similarity between COSMIC and simulation") +
  theme_minimal()  +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  )
  
out_file <- paste(out_dir, "signatures_validation.pdf", sep="")
ggsave(out_file, plot = p1, width = 180, height = 100, units = "mm")

out_file <- paste(out_dir, "signatures_validation.png", sep="")
ggsave(out_file, plot = p1, dpi = 600, width = 180, height = 100, units = "mm")


