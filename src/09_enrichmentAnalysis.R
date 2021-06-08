#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(tidyverse)
library(fgsea)
library(msigdbr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/"
resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

#-------------------------------------------------------------------------------
# Enrichment analysis - Hallmarks and GO-BP
#-------------------------------------------------------------------------------

# Enrichment analysis function -------------------------------------------------

performEnrichment = function(gene_sets, type)
{
  gene_sets = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  
  if(type == "Hallmarks"){
    pval_cutoff = 0.01
    w = 4
    h = 3
  }else if(type == "GOBP"){
    pval_cutoff = 0.001
    w = 8.5
    h = 15
  }
  
  # Enrichment analysis
  fgseaRes <- lapply(resFC_list, function(x) {
    enr = fgsea(pathways = gene_sets, 
                stats = x,
                eps = 0,
                minSize = 15,
                maxSize = 500)
    enr = enr[enr$padj < pval_cutoff,]
    enr = enr[order(enr$padj),][,c(1,2,3,6,8)]
    return(enr)
  })
  
  # Enrichment heatmap
  gseaHM = melt(lapply(fgseaRes, function(x) x[,c("pathway", "NES")]))
  gseaHM = gseaHM[,c(1,4,3)]
  gseaHM = dcast(data = gseaHM, formula = pathway ~ L1)
  
  gseaHM = gseaHM[apply(gseaHM, 1, function(x) sum(is.na(x))) < 11,]
  
  gseaHM = data.frame(gseaHM, row.names = 1)
  gseaHM[is.na(gseaHM)] = 0
  
  # Arranging the columns and rows as we want
  gseaHM = gseaHM[,c("IL22_3h", "IL22_6h", "IL22_12h", "IL22_24h", 
                     "IFNL_3h", "IFNL_6h", "IFNL_12h", "IFNL_24h",
                     "IL22_IFNL_3h", "IL22_IFNL_6h", "IL22_IFNL_12h", "IL22_IFNL_24h")]
  gseaHM = gseaHM[order(rowSums(gseaHM), decreasing = T),]
  
  # Heatmap design 
  myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
  br = c(seq(min(gseaHM), 0, length.out=ceiling(100/2) + 1), 
         seq(max(gseaHM)/100, max(gseaHM), length.out=floor(100/2)))
  
  if(type == "Hallmarks"){
    rownames(gseaHM) = gsub("HALLMARK_", "", rownames(gseaHM))
  }
  
  if(type == "GOBP"){
    rownames(gseaHM) = gsub("GO_", "", rownames(gseaHM))
  }
  
  col_anno <- data.frame(Treatment = c(rep("IL22",4), rep("IFNL", 4), rep("IL22_IFNL", 4)),
                         Time = rep(c("3h", "6h", "12h", "24h"), 3), row.names = colnames(gseaHM))
  
  ann_colors = list(
    Treatment = c(IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
    Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")
  )
  
  pheatmap(gseaHM, breaks = br,
           cluster_cols = F,
           cluster_rows = F,
           fontsize = 6,
           labels_row = gsub("_", " ", rownames(gseaHM), fixed = T),
           border_color = "white", 
           col = myColor,
           annotation_col = col_anno,
           annotation_colors = ann_colors,
           filename = paste0(path, "analysis/DiffExp/enriched_pathways_", type, ".pdf"), width = w, height = h)
  
  return(NULL)
}

# ------------------------------------------------------------------------------

# Ranked list for GSEA
resFC_list = apply(resFC, 2, function(x) {
  x = sort(setNames(x, rownames(resFC)))
  return(list(x))
})
resFC_list = lapply(resFC_list, unlist)

# Gene sets
hl_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
go_gene_sets = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")

performEnrichment(gene_sets = hl_gene_sets, type ="Hallmarks")
performEnrichment(gene_sets = go_gene_sets, type ="GOBP")


