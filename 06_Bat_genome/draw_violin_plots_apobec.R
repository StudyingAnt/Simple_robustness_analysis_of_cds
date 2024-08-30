# import package
library(tidyverse)
library(dplyr)
library(ggplot2)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/"

# read file
all_data <- list()
for (species in c("molMol2", "myoMyo6")) {
  for (mut_type in c("profile", "apobec", "apobec_diff")) {
    rum_file <- paste(base_path, species, "_rums_", mut_type, ".csv", sep="")
    all_data[[species]][[mut_type]] <- read.csv(rum_file)
  }
}

# data format
# Transcript_name Species Signature Robustness
transcripts <- c()
speciess <- c()
signatures <- c()
robustnesss <- c()
for (species in c("molMol2", "myoMyo6")) {
  for (mut_type in c("profile")) {
    n <- length(all_data[[species]][[mut_type]]$Transcript)
    
    transcripts <- c(transcripts, all_data[[species]][[mut_type]]$Transcript)
    speciess <- c(speciess, rep(species, n))
    signatures <- c(signatures, rep("SBS2", n))
    robustnesss <- c(robustnesss, all_data[[species]][[mut_type]]$SBS2)
  }
  
  for (mut_type in c("apobec")) {
    for (apobec in c("APOBEC3A", "APOBEC3C", "APOBEC1")) {
      n <- length(all_data[[species]][[mut_type]]$Transcript)
      
      transcripts <- c(transcripts, all_data[[species]][[mut_type]]$Transcript)
      speciess <- c(speciess, rep(species, n))
      signatures <- c(signatures, rep(paste(apobec, species, sep = "_"), n))
      robustnesss <- c(robustnesss, all_data[[species]][[mut_type]][,apobec])
    }
  }
  
  for (mut_type in c("apobec_diff")) {
    for (apobec in c("APOBEC3A", "APOBEC3C", "APOBEC1")) {
      n <- length(all_data[[species]][[mut_type]]$Transcript)
      
      transcripts <- c(transcripts, all_data[[species]][[mut_type]]$Transcript)
      speciess <- c(speciess, rep(species, n))
      if (species == "molMol2") {
        signatures <- c(signatures, rep(paste(apobec, "myoMyo6", sep = "_"), n))
      } else {
        signatures <- c(signatures, rep(paste(apobec, "molMol2", sep = "_"), n))
      }
      robustnesss <- c(robustnesss, all_data[[species]][[mut_type]][,apobec])
    }
  }
}

# human control SBS2
rum_file <- paste(base_path, "all_gencode_rums_profile.csv", sep="")
tmp.data <- read.csv(rum_file)

n <- length(tmp.data$Transcript)

transcripts <- c(transcripts, tmp.data$Transcript)
speciess <- c(speciess, rep("human", n))
signatures <- c(signatures, rep("SBS2", n))
robustnesss <- c(robustnesss, tmp.data$SBS2)

# human control bat apobecs
for (bat in c("molMol2", "myoMyo6")) {
  rum_file <- paste(base_path, "human_control_", bat, "_rums_apobec.csv", sep="")
  tmp.data <- read.csv(rum_file)
  
  n <- length(tmp.data$Transcript)
  
  for (apobec in c("APOBEC3A", "APOBEC3C", "APOBEC1")) {
    n <- length(tmp.data$Transcript)
    
    transcripts <- c(transcripts, tmp.data$Transcript)
    speciess <- c(speciess, rep("human", n))
    signatures <- c(signatures, rep(paste(apobec, bat, sep = "_"), n))
    robustnesss <- c(robustnesss, tmp.data[,apobec])
  }
}


data.plot <- data.frame(
  Transcript_name = transcripts,
  Species = speciess,
  Signature = signatures,
  Robustness = robustnesss
)

# Define the factor levels as normal text
data.plot$Signature <- factor(data.plot$Signature, 
                              levels = c("APOBEC1_molMol2", "APOBEC3A_molMol2", "APOBEC3C_molMol2",
                                         "APOBEC1_myoMyo6", "APOBEC3A_myoMyo6", "APOBEC3C_myoMyo6",
                                         "SBS2"))

data.plot$Species <- factor(data.plot$Species,
                            levels = c("human","molMol2","myoMyo6"),
                            labels = c("H. sapiens", "M. molossus", "M. myotis"))

fig_name <- c("APOBEC1_molMol2", "APOBEC3A_molMol2", "APOBEC3C_molMol2",
              "APOBEC1_myoMyo6", "APOBEC3A_myoMyo6", "APOBEC3C_myoMyo6",
              "SBS2")
sigs <- c("APOBEC1","APOBEC3A","APOBEC3C","APOBEC1","APOBEC3A","APOBEC3C","SBS2")
spes <- c("M. molossus", "M. molossus", "M. molossus", "M. myotis", "M. myotis", "M. myotis", "H. sapiens")
for (i in 1:7) {
  big.plot <- ggplot(data.plot %>% filter(Signature == fig_name[i]), aes(x=Species, y=Robustness, fill = Species)) + 
    geom_violin(width = 0.8, linewidth = 0.2) +
    geom_boxplot(width = 0.3, size = 0.2, outlier.size = 0.3, outlier.stroke = 0) +
    ylim(0,35) +
    #ggtitle(bquote(.(sigs[i]) ~ "("*italic(.(spes[i]))*")")) +
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
  
  small.plot <- ggplot(data.plot %>% filter(Signature == fig_name[i]), aes(x=Species, y=Robustness, fill = Species)) + 
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
  
  
  ggsave(paste(base_path, fig_name[i], "_allRange.pdf", sep = ""), plot = big.plot, width = 30, height = 60, units = "mm")
  ggsave(paste(base_path, fig_name[i], "_smallRange.pdf", sep = ""), plot = small.plot, width = 25, height = 47, units = "mm")
  
}

#######
origin.data.plot <- data.plot %>% filter((Species == "H. sapiens" & Signature == "SBS2") |
                                           (Species == "M. molossus" & Signature == "APOBEC3A_molMol2") |
                                           (Species == "M. myotis" & Signature == "APOBEC3A_myoMyo6") ) 


small_plot <- ggplot(origin.data.plot, aes(x=Species, y=Robustness, fill = Species)) + 
  geom_violin(width = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.3, size = 0.3, outlier.size = 0.5, outlier.stroke = 0) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_x_discrete(labels = c(expression(italic("H. sapiens")),
                              expression(italic("M. molossus")),
                              expression(italic("M. myotis"))
  )) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),  # Remove x-axis title
    axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),   # Remove x-axis text
    axis.text.y = element_text(size = 5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    plot.background = element_rect(colour = "red", fill = NA, size = 0.7)
  )

big_plot <- ggplot(origin.data.plot, aes(x=Species, y=Robustness, fill = Species)) + 
  geom_violin(width = 0.8, linewidth = 0.2) +
  geom_boxplot(width = 0.3, size = 0.2, outlier.size = 0.3, outlier.stroke = 0) +
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

ggsave(paste(base_path,"native_A3A_allRange.pdf", sep = ""), plot = big_plot, width = 30, height = 60, units = "mm")
ggsave(paste(base_path,"native_A3A_smallRange.pdf", sep = ""), plot = small_plot, width = 60, height = 90, units = "mm")

# statistical test

human.sbs2 <- origin.data.plot %>% filter(Species == "human")
molMol2.a3a <- origin.data.plot %>% filter(Species == "molMol2", Signature == "APOBEC3A_molMol2") 
myoMyo6.a3a <- origin.data.plot %>% filter(Species == "myoMyo6", Signature == "APOBEC3A_myoMyo6") 

wilcox.test(human.sbs2$Robustness, molMol2.a3a$Robustness, alternative = "greater")
wilcox.test(human.sbs2$Robustness, myoMyo6.a3a$Robustness, alternative = "less")
wilcox.test(molMol2.a3a$Robustness, myoMyo6.a3a$Robustness, alternative = "less")


