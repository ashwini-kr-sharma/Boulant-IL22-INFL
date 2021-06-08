#library(UpSetR)
#devtools::install_github("yanlinlin82/ggvenn")
#devtools::install_github("wjawaid/enrichR")

library(DESeq2)
library(ggvenn)
library(tidyverse)
library(patchwork)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(rio)

library(enrichR)
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs = c("TRRUST_Transcription_Factors_2019", 
        "MSigDB_Hallmark_2020")

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/"
dat = readRDS(paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))

#-------------------------------------------------------------------------------
# Top up-regulated or down-regulated genes
#-------------------------------------------------------------------------------

topgenes = lapply(dat, function(x){
  up = rownames(x)[x$log2FoldChange > 1 & x$padj < 0.05]
  up = up[!is.na(up)]
  
  dw = rownames(x)[x$log2FoldChange < -1 & x$padj < 0.05]
  dw = dw[!is.na(dw)]
  
  return(list(up = up, down = dw))
})
names(topgenes) = gsub("res_", "", names(topgenes))

up = lapply(topgenes, function(x) x$up)
dw = lapply(topgenes, function(x) x$down)

rm(topgenes)

#-------------------------------------------------------------------------------
# Venn diagrams UP
#-------------------------------------------------------------------------------

all.up = list(IL22 = unique(unlist(up[1:4])), 
           IFNL = unique(unlist(up[5:8])), 
           `IL22+IFNL` = unique(unlist(up[9:12])) 
           )
  
p_IL22.up = ggvenn(up[1:4], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.up$IL22)))

p_IFNL.up = ggvenn(up[5:8], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.up$IFNL)))

p_comb.up = ggvenn(up[9:12], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.up$`IL22+IFNL`)))

p_all.up = ggvenn(all.up, fill_color = c("#fdb462", "#80b1d3", "#b3de69"), 
                  stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F)

#-------------------------------------------------------------------------------
# Venn diagrams DW
#-------------------------------------------------------------------------------

all.dw = list(IL22 = unique(unlist(dw[1:4])), 
              IFNL = unique(unlist(dw[5:8])), 
             `IL22+IFNL` = unique(unlist(dw[9:12])) 
)

p_IL22.dw = ggvenn(dw[1:4], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                   stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.dw$IL22)))

p_IFNL.dw = ggvenn(dw[5:8], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                   stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.dw$IFNL)))

p_comb.dw = ggvenn(dw[9:12], fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                   stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
  labs(subtitle = paste0("n=",length(all.dw$`IL22+IFNL`)))

p_all.dw = ggvenn(all.dw, fill_color = c("#fdb462", "#80b1d3", "#b3de69"), 
                  stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F)

#-------------------------------------------------------------------------------
# Overlapping genes
#-------------------------------------------------------------------------------

# Enrichment function ----------------------------------------------------------

enrichrAnalysis = function(genes, type){
  
  if(length(genes) >= 20){
    
    if(type == "IL22"){
      colID = "#fdb462"
    }else if(type == "IFNL"){
      colID = "#80b1d3"
    }else if(type == "IL22+IFNL"){
      colID = "#b3de69"
    }
    
    res <- melt(lapply(enrichr(genes, dbs), function(x){
      x = x[,c(1,4,7)]
      x = x[x$Adjusted.P.value < 0.01,]
    }), id = c("Term", "Adjusted.P.value", "Odds.Ratio"))
    
    if(nrow(res) != 0){
      res$Adjusted.P.value = -log10(res$Adjusted.P.value)
      res$Odds.Ratio = log10(res$Odds.Ratio)
      res = res[grep("mouse", res$Term, invert = T),]
      res$Adjusted.P.value[res$Adjusted.P.value > 10] = 10
      
      p <- ggplot(data = res, aes(x = Adjusted.P.value, y = reorder(Term, Adjusted.P.value))) + 
        theme_bw(base_size = 8) + labs(y = "", subtitle = paste0(type,"(n=", length(genes),")")) + geom_bar(stat="identity", fill = colID) + xlim(0,10) +
        facet_wrap(~ L1, scales ="free", ncol = 1) + geom_vline(xintercept = -log10(0.01)) +
        theme(panel.grid = element_blank(), axis.line = element_blank(), strip.background = element_blank())
    }else{
      p = NULL
    }
  }else{
    p = NULL
  }
  return(p)
}

#-------------------------------------------------------------------------------

#######################
# Over expressed genes
#######################

a = intersect(all.up$IFNL, all.up$IL22)
b = all.up$IFNL[! all.up$IFNL %in% a] # IFNL modulated
c = all.up$IL22[! all.up$IL22 %in% a] # IL22 modulated
d = intersect(a, all.up$`IL22+IFNL`)
e = all.up$`IL22+IFNL`[! all.up$`IL22+IFNL` %in% c(all.up$IFNL, all.up$IL22)] # IL22 + IFNL modulated

enr_IFNL_up = enrichrAnalysis(genes = b, type = "IFNL") # IFNL modulated
enr_IL22_up = enrichrAnalysis(genes = c, type = "IL22") # IL22 modulated
enr_IL22_IFNL_up = enrichrAnalysis(genes = e, type = "IL22+IFNL") # IL22 + IFNL modulated

rm(a,b,c,d,e)

#######################
# Down expressed genes
#######################

a = intersect(all.dw$IFNL, all.dw$IL22)
b = all.dw$IFNL[! all.dw$IFNL %in% a] # IFNL modulated
c = all.dw$IL22[! all.dw$IL22 %in% a] # IL22 modulated
d = intersect(a, all.dw$`IL22+IFNL`)
e = all.dw$`IL22+IFNL`[! all.dw$`IL22+IFNL` %in% c(all.dw$IFNL, all.dw$IL22)] # IL22 + IFNL modulated

enr_IFNL_dw = enrichrAnalysis(genes = b, type = "IFNL") # IFNL modulated -  Only 4 genes - "GPS2", "TNFRSF10D", "MAOB", "ACO1"  
enr_IL22_dw = enrichrAnalysis(genes = c, type = "IL22") # IL22 modulated
enr_IL22_IFNL_dw = enrichrAnalysis(genes = e, type = "IL22+IFNL") # IL22 + IFNL modulated

# There were no significant Enrichment observed

rm(a,b,c,d,e,enr_IL22_dw,enr_IL22_IFNL_dw,enr_IFNL_dw)

#-------------------------------------------------------------------------------
# Consolidated plots
#-------------------------------------------------------------------------------

# Over expressed gene
p_up = p_IL22.up | p_IFNL.up | p_comb.up | p_all.up
enr_up = enr_IL22_up | enr_IFNL_up | enr_IL22_IFNL_up

p = p_up / enr_up + plot_layout(heights = c(1, 1.5))
ggsave(filename = paste0(path, "analysis/DiffExp/OverExpressedGenes.pdf"), plot = p, width = 11, height = 7)

rm(p_IL22.up, p_IFNL.up, p_comb.up, p_all.up, enr_IFNL_up, enr_IL22_up, enr_IL22_IFNL_up)

# Under expressed gene
p_dw = p_IL22.dw | p_IFNL.dw | p_comb.dw | p_all.dw
ggsave(filename = paste0(path, "analysis/DiffExp/UnderExpressedGenes.pdf"), plot = p_dw, width = 11, height = 4)

rm(p_IL22.dw, p_IFNL.dw, p_comb.dw, p_all.dw)

#-------------------------------------------------------------------------------
# Fold change analysis between IL22 and IFN-L cases
#-------------------------------------------------------------------------------

resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

# Sergio list
# isg = read.table("/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/data/Sergio_ISG_list.txt")
# isg = unique(c(isg[,1], rownames(resFC)[grepl("^IFN", rownames(resFC))],
#                rownames(resFC)[grepl("^ISG", rownames(resFC))],
#                rownames(resFC)[grepl("^IFI", rownames(resFC))],
#                rownames(resFC)[grepl("^IRG", rownames(resFC))]))

# ISGs from https://doi.org/10.1038/s41590-019-0323-3 (Supplementary Table 1)
isg = import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0323-3/MediaObjects/41590_2019_323_MOESM3_ESM.xlsx", sheet = 3)
isg = isg[,1:4]
colnames(isg) = isg[1,]
isg = isg[-1,]
isg = isg$ISG
#isg = isg$ISG[which(isg$`Used as Bait` == "+")]

# For clearer visualization abls(log2FC) > 2 has been floored to +/- 2
resFC[resFC > 2] = 2
resFC[resFC < -2] = -2

# isg = resFC[rownames(resFC) %in% isg,]
# isg = isg[rowSums(abs(isg) > 0.2) > 0,]
# isg = rownames(isg)

colnames(resFC) = paste0(colnames(resFC),"_logFC")
resFC = data.frame(Gene = rownames(resFC), isISG = "Non-ISG", resFC[,1:8])
resFC$isISG[rownames(resFC) %in% isg] = "ISG"

p_3h = ggplot(data = resFC, aes(x = IL22_3h_logFC, y = IFNL_3h_logFC, color = isISG)) + 
  theme_bw(base_size = 8) + labs(color = "") +
  geom_point(data = subset(resFC, isISG == "Non-ISG"), size = 0.3) + 
  geom_point(data = subset(resFC, isISG == "ISG"), size = 0.5) + 
  ylim(c(-2,2)) + xlim(c(-2,2)) + scale_color_manual(values = c("firebrick", "grey80")) +
  geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = c(-1,1)) + 
  theme(panel.grid = element_blank(), legend.position = "top")

p_6h = ggplot(data = resFC, aes(x = IL22_6h_logFC, y = IFNL_6h_logFC, color = isISG)) + 
  theme_bw(base_size = 8) + labs(color = "") +
  geom_point(data = subset(resFC, isISG == "Non-ISG"), size = 0.3) + 
  geom_point(data = subset(resFC, isISG == "ISG"), size = 0.5)  + 
  ylim(c(-2,2)) + xlim(c(-2,2)) + scale_color_manual(values = c("firebrick", "grey80")) +
  geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = c(-1,1)) + 
  theme(panel.grid = element_blank(), legend.position = "top")

p_12h = ggplot(data = resFC, aes(x = IL22_12h_logFC, y = IFNL_12h_logFC, color = isISG)) + 
  theme_bw(base_size = 8) + labs(color = "") +
  geom_point(data = subset(resFC, isISG == "Non-ISG"), size = 0.3) + 
  geom_point(data = subset(resFC, isISG == "ISG"), size = 0.5) + 
  ylim(c(-2,2)) + xlim(c(-2,2)) + scale_color_manual(values = c("firebrick", "grey80")) +
  geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = c(-1,1)) + 
  theme(panel.grid = element_blank(), legend.position = "top")

p_24h = ggplot(data = resFC, aes(x = IL22_24h_logFC, y = IFNL_24h_logFC, color = isISG)) + 
  theme_bw(base_size = 8) + labs(color = "") +
  geom_point(data = subset(resFC, isISG == "Non-ISG"), size = 0.3) + 
  geom_point(data = subset(resFC, isISG == "ISG"), size = 0.5) + 
  ylim(c(-2,2)) + xlim(c(-2,2)) + scale_color_manual(values = c("firebrick", "grey80")) +
  geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = c(-1,1)) + 
  theme(panel.grid = element_blank(), legend.position = "top")

p = p_3h + p_6h + p_12h + p_24h + plot_layout(ncol = 2, nrow = 2, guides = "collect") & theme(legend.position='top')
rm(p_3h, p_6h, p_12h, p_24h)

ggsave(filename = paste0(path, "analysis/DiffExp/IL22_IFNL_logFC_compare.pdf"), plot = p, width = 5, height = 4)

#-------------------------------------------------------------------------------
# Plot heatmap of the top genes from each condition
#-------------------------------------------------------------------------------

cnt = readRDS(paste0(path, "data/counts/T84_IL22_INFL_filtered_counts.RDS"))
diffExp = readRDS(paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))

dds = DESeqDataSetFromMatrix(countData = cnt$counts,
                             colData   = cnt$sampanno,
                             rowData   = cnt$geneanno,
                             design    = ~ Condition)

rm(cnt)
rld = rlog(dds, blind = TRUE)
rld = assay(rld)
sample_order = c(paste0("Mock-3h-Rep-",1:3), paste0("Mock-6h-Rep-",1:3),  paste0("Mock-12h-Rep-",1:3), paste0("Mock-24h-Rep-",1:3), 
                 paste0("IFNL-3h-Rep-",1:3),  paste0("IFNL-6h-Rep-",1:3),  paste0("IFNL-12h-Rep-",1:3),  paste0("IFNL-24h-Rep-",1:3),
                 paste0("IL22-3h-Rep-",1:3), paste0("IL22-6h-Rep-",1:3),  paste0("IL22-12h-Rep-",1:3), paste0("IL22-24h-Rep-",1:3),
                 paste0("IL22_IFNL-3h-Rep-",1:3), paste0("IL22_IFNL-6h-Rep-",1:3), paste0("IL22_IFNL-12h-Rep-",1:3), paste0("IL22_IFNL-24h-Rep-",1:3))
rld = rld[,sample_order]

topgenes = lapply(diffExp, function(x){
  
  x = as.data.frame(x)
  x = na.omit(x)
  
  up = x[x$log2FoldChange > 1 & x$padj < 0.01,]
  up = up[order(up$log2FoldChange, decreasing = T),]
  up = head(up, 20)
  up = rownames(up)
  
  # dw = x[x$log2FoldChange < -1 & x$padj < 0.01,]
  # dw = dw[order(dw$log2FoldChange, decreasing = F),]
  # dw = head(dw, 25)
  # dw = rownames(dw)
  # 
  # selgenes = c(up, dw)
  
  selgenes = up
  
  return(selgenes)
})

tmpdat1 = rld[rownames(rld) %in% unique(unlist(topgenes[1:8])),]
dim(tmpdat1)
#[1] 76  48

orderID = sort(rowMeans(tmpdat1[,13:24])/rowMeans(tmpdat1[,25:36]), decreasing = T)
tmpdat1 = tmpdat1[names(orderID),]

datanno = data.frame(do.call("rbind", strsplit(colnames(tmpdat1), "-"))[,c(2,1)])
colnames(datanno) = c("Time", "Treatment")
rownames(datanno) = colnames(tmpdat1)

ann_colors = list(Treatment = c(Mock = "#fb8072", IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
                  Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525"))

pheatmap(t(tmpdat1), scale = "column", cluster_rows = F, cluster_cols = F, 
         show_rownames = F, show_colnames = T,
         fontsize = 6, border_color = "white",
         gaps_row = c(12, 24, 36),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
         annotation_row = datanno, annotation_colors = ann_colors,
         filename = paste0(path, "analysis/DiffExp/mostUPregulated_top20genes_from_IL22_IFNL_groups.pdf"),
         width = 7.5, height = 2.5)

rm(dds, tmpdat1, orderID, datanno, ann_colors, rld, sample_order, diffExp, topgenes)
