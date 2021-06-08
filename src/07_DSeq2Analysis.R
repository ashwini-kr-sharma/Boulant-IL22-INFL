#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(DESeq2)
library(WriteXLS)
library(ashr)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/"
dat  = readRDS(paste0(path, "data/counts/T84_IL22_INFL_filtered_counts.RDS"))

# Remove the outlier sample
dat$counts = dat$counts[, !colnames(dat$counts) == "IL22−24h−Rep−3"]
dat$sampanno = dat$sampanno[!dat$sampanno$SampleID == "IL22−24h−Rep−3",]

#-------------------------------------------------------------------------------
# Create directories and download two independent signatures for stemness
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/DiffExp"))){
  dir.create(paste0(path, "analysis/DiffExp"),recursive = T)
}

#-------------------------------------------------------------------------------
# DESeq2 object - All conditions logFC
#-------------------------------------------------------------------------------

ddsFC = DESeqDataSetFromMatrix(countData = dat$counts,
                               colData = dat$sampanno,
                               rowData = dat$geneanno,
                               design = ~ Condition)
rm(dat)

ddsFC = DESeq(ddsFC)

# Function to compute all possible contrasts with Mock
diffexp = function(cond1, cond2){
  dat = results(ddsFC, contrast = c("Condition", cond1, cond2), alpha = 0.01)
  dat = lfcShrink(ddsFC, res = dat, contrast = c("Condition", cond1, cond2), type="ashr")
}

# IL22
res_3hr_IL22  = diffexp(cond1 = "IL22-3h", cond2 = "Mock-3h")
res_6hr_IL22  = diffexp(cond1 = "IL22-6h", cond2 = "Mock-6h")
res_12hr_IL22 = diffexp(cond1 = "IL22-12h", cond2 = "Mock-12h")
res_24hr_IL22 = diffexp(cond1 = "IL22-24h", cond2 = "Mock-24h")

# IFN-L
res_3hr_IFNL  = diffexp(cond1 = "IFNL-3h", cond2 = "Mock-3h")
res_6hr_IFNL  = diffexp(cond1 = "IFNL-6h", cond2 = "Mock-6h")
res_12hr_IFNL = diffexp(cond1 = "IFNL-12h", cond2 = "Mock-12h")
res_24hr_IFNL = diffexp(cond1 = "IFNL-24h", cond2 = "Mock-24h")

# IL22 + IFN-L
res_3hr_IL22_IFNL  = diffexp(cond1 = "IL22_IFNL-3h", cond2 = "Mock-3h")
res_6hr_IL22_IFNL  = diffexp(cond1 = "IL22_IFNL-6h", cond2 = "Mock-6h")
res_12hr_IL22_IFNL = diffexp(cond1 = "IL22_IFNL-12h", cond2 = "Mock-12h")
res_24hr_IL22_IFNL = diffexp(cond1 = "IL22_IFNL-24h", cond2 = "Mock-24h")

# Save differential gene expression results

dge = list(res_3hr_IL22 = res_3hr_IL22, res_6hr_IL22 = res_6hr_IL22, res_12hr_IL22 = res_12hr_IL22, res_24hr_IL22 = res_24hr_IL22,
           res_3hr_IFNL = res_3hr_IFNL, res_6hr_IFNL = res_6hr_IFNL, res_12hr_IFNL = res_12hr_IFNL, res_24hr_IFNL = res_24hr_IFNL,
           res_3hr_IL22_IFNL = res_3hr_IL22_IFNL, res_6hr_IL22_IFNL = res_6hr_IL22_IFNL, res_12hr_IL22_IFNL = res_12hr_IL22_IFNL, res_24hr_IL22_IFNL = res_24hr_IL22_IFNL)

saveRDS(dge, paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))

dge = lapply(dge, function(x){
  x = data.frame(x)
  x[,c(1:3)] = round(x[,c(1:3)], 2)
  x$pvalue = signif(x$pvalue, 2)
  x$padj = signif(x$padj, 2)
  return (x)
})

WriteXLS(lapply(dge, function(x){data.frame(Gene = rownames(x), x)}), 
         ExcelFileName = paste0(path, "analysis/DiffExp/DGEtables.xls"),
         AdjWidth = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1, 
         SheetNames = gsub("res_", "", names(dge)))

rm(dge)

# Consolidated FC 
resFC = data.frame(IL22_3h = res_3hr_IL22$log2FoldChange, 
                   IL22_6h = res_6hr_IL22$log2FoldChange, 
                   IL22_12h = res_12hr_IL22$log2FoldChange, 
                   IL22_24h = res_24hr_IL22$log2FoldChange,
                   
                   IFNL_3h = res_3hr_IFNL$log2FoldChange,
                   IFNL_6h = res_6hr_IFNL$log2FoldChange,
                   IFNL_12h = res_12hr_IFNL$log2FoldChange,
                   IFNL_24h = res_24hr_IFNL$log2FoldChange,
                   
                   IL22_IFNL_3h = res_3hr_IL22_IFNL$log2FoldChange,
                   IL22_IFNL_6h = res_6hr_IL22_IFNL$log2FoldChange,
                   IL22_IFNL_12h = res_12hr_IL22_IFNL$log2FoldChange,
                   IL22_IFNL_24h = res_24hr_IL22_IFNL$log2FoldChange,
                   
                   row.names = rownames(ddsFC))

saveRDS(resFC, paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

rm(res_3hr_IL22, res_6hr_IL22, res_12hr_IL22, res_24hr_IL22,
   res_3hr_IFNL, res_6hr_IFNL, res_12hr_IFNL, res_24hr_IFNL,
   res_3hr_IL22_IFNL, res_6hr_IL22_IFNL, res_12hr_IL22_IFNL, res_24hr_IL22_IFNL)

#-------------------------------------------------------------------------------
# DESeq2 object - all Time course
#-------------------------------------------------------------------------------

# ddsTC = DESeqDataSetFromMatrix(countData = dat$counts,
#                                colData   = dat$sampanno,
#                                rowData   = dat$geneanno,
#                                design    = ~ Type + Time + Type:Time)
# ddsTC = DESeq(ddsTC, test = "LRT", reduced = ~ Type + Time)
# 
# if(identical(rownames(resTC), rowData(ddsTC)$GeneID)){
#   resTC = results(ddsTC, alpha = 0.01)
#   resTC$Symbol = rowData(ddsTC)$gene_name
#   resTC = as.data.frame(resTC)
#   resTC = resTC[order(resTC$log2FoldChange, resTC$padj),]
# }
# 
# resTC[resTC$Symbol %in% hif, ]
# 
# sel = resTC[which(resTC$Symbol %in% hypGO &
#             resTC$log2FoldChange > 1 & 
#             resTC$padj < 0.01), ]
# 
# fiss <- plotCounts(ddsTC, gene = "ENSG00000112715", 
#                    intgroup = c("Time","Type"), returnData = TRUE)
# 
# ggplot(fiss,
#        aes(x = Time, y = count, color = Type, group = Type)) + 
#   geom_point() + stat_summary(fun=median, geom="line") +
#   scale_y_log10()
