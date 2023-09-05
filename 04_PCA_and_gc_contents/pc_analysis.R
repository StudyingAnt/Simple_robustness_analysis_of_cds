library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(egg)
library(readxl)
library(xlsx)

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

aetiology_file <- paste(base_path, "Table_S1.xlsx", sep="")
aetiologies <- read_excel(aetiology_file)
aetiologies <- aetiologies %>% filter(!row_number() %in% c(1,2,8,12,18,20,26,34))

colnames(aetiologies)[1] = "Signature"
colnames(aetiologies)[2] = "Aetiology"

wb <- createWorkbook(type="xlsx")

for (i in 1:length(signatures_noseqerr)){
  pc_tosee <- paste("PC", as.character(i), sep = "")
  
  table_data <- pc_matrix %>% filter(Principal_component == pc_tosee) %>% 
    filter(Signature %in% signatures_aetiology) %>% 
    mutate(Signature = factor(Signature, signatures_aetiology)) %>%
    arrange(Signature)
  table_data$Aetiology <- aetiologies$Aetiology
  table_data <- table_data %>% arrange(desc(Coefficient))
  
  sheet <- createSheet(wb, pc_tosee)
  addDataFrame(table_data, sheet = sheet, row.names = FALSE)

  plot_data <- pc_matrix %>% filter(Principal_component == pc_tosee) %>% arrange(Coefficient) %>% filter(Signature %in% signatures_aetiology)
  sig_order <- plot_data$Signature
  
  sig_plot <- ggplot(plot_data, aes(x=factor(Signature, levels = sig_order), y=Coefficient, fill=Coefficient)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="#faedcd",
                         high="red") +
    xlab("Signature") +
    ylab("Coefficient") +
    labs(fill="Coefficient") +
    theme_light() +
    theme(
      axis.text.x = element_text(size = 6, family = "Arial", angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 8, family = "Arial"),
      legend.title = element_text(size = 7, family = "Arial"),
      legend.text = element_text(size = 6, family = "Arial"),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.box.margin = margin(-10,-10,-10,-10),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) 
  
  
  out_file <- paste(base_path, pc_tosee, "_sorted_signatures.png", sep="")
  ggsave(out_file, plot = sig_plot, dpi = 1200, width = 50, height = 115, units = "mm")
}

wb_file <- paste(base_path, "Table_S2.xlsx", sep = "")
saveWorkbook(wb, wb_file)
