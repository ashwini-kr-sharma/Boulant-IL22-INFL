#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/"
resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))
rlog  = readRDS(paste0(path, "analysis/Exploratory/rlogTransformation.RDS"))

#-------------------------------------------------------------------------------
# IL22 and IFN-L signalling related proteins logFC
#-------------------------------------------------------------------------------

infr = resFC[rownames(resFC) %in% c("JAK1", "TYK2", 
                                    "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
                                    "IFNLR1", "IL10RB",
                                    "IL22RA1", "IL22RA2"),]

# # Keeping only those genes with abs log2FC > 1 in at least one condition
# sel = apply(infr, 1, function(x){
#   sum(abs(x) > 0.5)
# })
# infr = infr[sel > 0 ,]
# rm(sel)

myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
br = c(seq(min(infr), 0, length.out=ceiling(100/2) + 1), 
       seq(max(infr)/100, max(infr), length.out=floor(100/2)))

infr = infr[c("JAK1", "TYK2", 
              "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
              "IFNLR1", "IL10RB",
              "IL22RA1", "IL22RA2")
            ,
            
            c("IL22_3h", "IL22_6h", "IL22_12h", "IL22_24h", 
               "IFNL_3h", "IFNL_6h", "IFNL_12h", "IFNL_24h",
               "IL22_IFNL_3h", "IL22_IFNL_6h", "IL22_IFNL_12h", "IL22_IFNL_24h")]

col_anno <- data.frame(Treatment = c(rep("IL22",4), rep("IFNL", 4), rep("IL22_IFNL", 4)),
                       Time = rep(c("3h", "6h", "12h", "24h"), 3), row.names = colnames(infr))

ann_colors = list(
  Treatment = c(IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
  Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")
)

pheatmap(infr,
         breaks = br,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 6,
         border_color = "white", 
         col = myColor,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         filename = paste0(path, "analysis/DiffExp/inf_response_receptors_logFC.pdf"), width = 2.8, height = 2.5)

rm(br, myColor)

#-------------------------------------------------------------------------------
# IL22 and IFN-L signalling related proteins expression levels
#-------------------------------------------------------------------------------

infrexp = melt(rlog)
infrexp = cbind(infrexp, do.call("rbind", strsplit(as.character(infrexp$Var2), "-", fixed=T))[,1:2])
infrexp = infrexp[,c(1,4,5,3)]
colnames(infrexp) = c("Gene", "Treatment", "Time", "Expression")

infrexp$Time = factor(as.character(infrexp$Time), 
                      levels = c("3h", "6h", "12h", "24h"))

infrexp$Treatment = factor(as.character(infrexp$Treatment), 
                           levels = c("Mock", "IFNL", "IL22", "IL22_IFNL"))

infrexp$Gene[! infrexp$Gene %in% c("JAK1", "TYK2", 
                                   "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
                                   "IFNLR1", "IL10RB",
                                   "IL22RA1", "IL22RA2")] = NA

infrexp$Gene = factor(as.character(infrexp$Gene), 
                      levels = c("JAK1", "TYK2", 
                                 "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
                                 "IFNLR1", "IL10RB",
                                 "IL22RA1", "IL22RA2"))

p1 = ggplot(infrexp, aes(x = Treatment, y = Expression)) + 
  theme_classic(base_size = 9) +
  labs(x = "", y = "rlog expression (~18,500 genes)") +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(aes(colour = Gene), width = 0.1) + 
  scale_color_brewer(palette="Set3", na.translate=FALSE) +
  facet_grid(~Time) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

ggsave(filename = paste0(path, "analysis/DiffExp/inf_response_receptors_absolute_expression.pdf"), plot = p1, width = 5, height = 3)
