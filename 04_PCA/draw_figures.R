library(tidyverse)
library(ggplot2)
library(dplyr)

# file directories
scripts_dir <- getwd()
scripts_dir
dir_names <- strsplit(scripts_dir, split="/")[[1]]
exp_dir <- paste(head(dir_names, -1), collapse = "/") 
in_dir <- paste(exp_dir, "in", sep = "/")
out_dir <- paste(exp_dir, "out", sep = "/")
fig_dir <- paste(exp_dir, "figures", sep = "/")

input_file <- paste(out_dir, "rus_profile_umap.csv", sep="/")

data <- read.csv(input_file)
p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point()
p

xy_file <- paste(out_dir, "xy.recol.csv", sep = "/")
pca_xy <- read.csv(xy_file)
z_file <- paste(out_dir, "z.recol.csv", sep = "/")
pca_z <- read.csv(z_file, header = FALSE)

x <- pca_xy$PC
y <- pca_xy$Sorted_SBS
data <- expand.grid(X=x,Y=y)
data$Z <- pca_z[[1]]

top5 <- c()
for (i in 1:5){
  top5[i] <- paste("PC", as.character(i), sep = "")
}

top10 <- c()
for (i in 1:10){
  top10[i] <- paste("PC", as.character(i), sep = "")
}

top12 <- c()
for (i in 1:12){
  top12[i] <- paste("PC", as.character(i), sep = "")
}

top17 <- c()
for (i in 1:17){
  top17[i] <- paste("PC", as.character(i), sep = "")
}

top21 <- c()
for (i in 1:21){
  top21[i] <- paste("PC", as.character(i), sep = "")
}

p <- ggplot(data, aes(factor(X),factor(Y),fill=Z)) +
  geom_tile() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                        high="red") +
  #xlim(top17)+
  xlab("Principal Component") +
  ylab("Signature Sorted by GC Bias") +
  labs(fill="Coefficient")  +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    axis.text.x = element_text(size = 5, family = "Arial", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 5, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 5, family = "Arial"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    )

p

out_file <- paste(fig_dir, "pc_all.recol.v2.png", sep="/")
ggsave(out_file, plot = p, width = 108, height = 95, units = "mm")

################################################################################
pca_file <- "/home/augustine/Desktop/Labs/0_Dry/projects/Cancer_Evolution_and_Mutational_Signature/exp/20221214_all_transcripts/out/data1.csv"

pca_data <- read.csv(pca_file)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=gc_wobble)) +
  geom_point(size = 0.0001, shape=20) +
  scale_color_gradientn(colors = jet.colors(7), limits = c(0, 1)) +
  xlab("PC1") +
  ylab("PC2") +
  labs(color="GC Wobble")  +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 5, family = "Arial"),
    axis.text.y = element_text(size = 5, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 5, family = "Arial"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

p <- ggplot(pca_data, aes(x=PC1, y=PC2)) +
  geom_point(size = 0.0001, shape=20) +
  xlab("PC1") +
  ylab("PC2") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 5, family = "Arial"),
    axis.text.y = element_text(size = 5, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

p
fig_dir <- "/home/augustine/Desktop/Labs/0_Dry/projects/Cancer_Evolution_and_Mutational_Signature/exp/20221214_all_transcripts/figures"
out_file <- paste(fig_dir, "pca_only.png", sep="/")
ggsave(out_file, plot = p, dpi=900, width = 100, height = 85, units = "mm")
