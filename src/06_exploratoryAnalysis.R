#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(progeny)
library(dorothea)
library(patchwork)
library(ggpubr)
library(reshape2)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/T84_IL22_IFNL/"
dat  = readRDS(paste0(path, "data/counts/T84_IL22_INFL_filtered_counts.RDS"))

#-------------------------------------------------------------------------------
# Create directories and download two independent signatures for stemness
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/Exploratory"))){
  dir.create(paste0(path, "analysis/Exploratory"),recursive = T)
}

#-------------------------------------------------------------------------------
# DESeq2 object
#-------------------------------------------------------------------------------

dds = DESeqDataSetFromMatrix(countData = dat$counts,
                             colData   = dat$sampanno,
                             rowData   = dat$geneanno,
                             design    = ~ Condition)

rm(dat)

#-------------------------------------------------------------------------------
# Data transformation - rlog
#-------------------------------------------------------------------------------

rld = rlog(dds, blind = TRUE)

#-------------------------------------------------------------------------------
# Sample similarity heatmap based on Euclidean distance
#-------------------------------------------------------------------------------

# Compute distance
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colData(dds)$SampleID
colnames(sampleDistMatrix) = colData(dds)$SampleID

# Annotations for plotting
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(10)
col_anno = as.data.frame(colData(dds), row.names = 1 )[,1:2]
col_anno$Library_Size = colSums(dds@assays@data$counts)/10^6

identical(colnames(sampleDistMatrix), rownames(col_anno))
#[1] TRUE

ann_colors = list(
  Treatment = c(Mock = "#fb8072", IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
  Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")
)

# Heatmap
pheatmap(sampleDistMatrix, clustering_method = "ward.D2",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cutree_rows = 5, cutree_cols = 5, fontsize = 6,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         border_color = NA,
         col = colors,
         treeheight_row = 15, treeheight_col = 15,
         filename = paste0(path, "analysis/Exploratory/sample_similarity_heatmap_cluster.pdf"), width = 7, height = 6)

# Heatmap
sample_order = c(paste0("Mock-3h-Rep-",1:3), paste0("IL22-3h-Rep-",1:3), paste0("IFNL-3h-Rep-",1:3), paste0("IL22_IFNL-3h-Rep-",1:3),
                 paste0("Mock-6h-Rep-",1:3), paste0("IL22-6h-Rep-",1:3), paste0("IFNL-6h-Rep-",1:3), paste0("IL22_IFNL-6h-Rep-",1:3),
                 paste0("Mock-12h-Rep-",1:3), paste0("IL22-12h-Rep-",1:3), paste0("IFNL-12h-Rep-",1:3), paste0("IL22_IFNL-12h-Rep-",1:3),
                 paste0("Mock-24h-Rep-",1:3), paste0("IL22-24h-Rep-",1:3), paste0("IFNL-24h-Rep-",1:3), paste0("IL22_IFNL-24h-Rep-",1:3))

# Sanity check !!
identical(sort(sample_order), sort(colnames(sampleDistMatrix)))
#[1] TRUE

sampleDistMatrix = sampleDistMatrix[sample_order, sample_order]

# Sanity check !!
identical(sort(sample_order), sort(rownames(col_anno)))
#[1] TRUE

col_anno = col_anno[sample_order,]

# Sanity check !!
identical(colnames(sampleDistMatrix), rownames(col_anno))
#[1] TRUE

pheatmap(sampleDistMatrix, 
         cluster_rows = F, cluster_cols = F,
         fontsize = 6,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         border_color = NA, 
         col = colors,
         filename = paste0(path, "analysis/Exploratory/sample_similarity_heatmap_orderedBy_Name.pdf"), width = 7, height = 6)

rm(colors, ann_colors, sample_order, sampleDists)

#-------------------------------------------------------------------------------
# Multidimensional scaling (MDS) plot
#-------------------------------------------------------------------------------

mds <- cbind(col_anno, cmdscale(sampleDistMatrix))

ggplot(mds, aes(x = `1`, y = `2`, color = Time, shape = Treatment)) + 
theme_bw(base_size = 8) + labs(x = "MDS dimension 1", y = "MDS dimension 2") +
geom_point(size = 2) +  ggtitle("MDS - rlog transformed data") +
scale_color_manual(values = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")) +
theme(panel.grid = element_blank(),
      axis.text = element_text(colour="black"), 
      axis.line = element_blank(),
      legend.position = "right") +
ggsave(filename = paste0(path, "analysis/Exploratory/sample_similarity_MDS.pdf"), width = 4, height = 3)

rm(sampleDistMatrix, mds, col_anno)

#-------------------------------------------------------------------------------
# Removing the outlier sample
#-------------------------------------------------------------------------------

expr = assay(rld)
expr = expr[, !colnames(expr) %in% "IL22−24h−Rep−3"]
saveRDS(object = expr, file = paste0(path, "analysis/Exploratory/rlogTransformation.RDS"))

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Compute pathway activity scores
#-------------------------------------------------------------------------------

pr_res = progeny(expr)

write.csv(data.frame(Samples = rownames(pr_res), pr_res, stringsAsFactors = F),
          paste0(path, "analysis/Exploratory/progeny_all_results.csv"), 
          quote = FALSE, row.names = F)

pr_res_summ = split(data.frame(pr_res), 
               sapply(strsplit(rownames(pr_res), "-Rep-"), function(x) x[1]))

pr_res_summ = t(sapply(pr_res_summ, function(x){
  apply(x, 2, median)
}))

pr_res_summ = pr_res_summ[c("Mock-3h", "Mock-6h", "Mock-12h", "Mock-24h",
                  "IL22-3h", "IL22-6h", "IL22-12h", "IL22-24h",
                  "IFNL-3h", "IFNL-6h", "IFNL-12h", "IFNL-24h",
                  "IL22_IFNL-3h", "IL22_IFNL-6h", "IL22_IFNL-12h", "IL22_IFNL-24h"),]

paletteLength = 100
cols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(paletteLength)
brk = c(seq(min(pr_res_summ), 0, length.out=ceiling(paletteLength/2) + 1),
        seq(max(pr_res_summ)/paletteLength,  max(pr_res_summ), length.out=floor(paletteLength/2)))

col_anno <- data.frame(do.call("rbind", strsplit(rownames(pr_res_summ), "-")), row.names = rownames(pr_res_summ))
colnames(col_anno) = c("Treatment", "Time")

# Sanity check
identical(rownames(col_anno), rownames(pr_res_summ))

ann_colors = list(
  Treatment = c(Mock = "#fb8072", IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
  Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")
)

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Heatmap of activity scores
#-------------------------------------------------------------------------------

pheatmap(pr_res_summ, clustering_method = "ward.D2",
         cluster_rows = F,
         cutree_cols = 3, 
         fontsize = 6,
         annotation_row = col_anno,
         annotation_colors = ann_colors,
         border_color = "white", 
         col = cols,
         breaks = brk,
         treeheight_row = 15, treeheight_col = 15,
         filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_heatmap.pdf"), width = 3.5, height = 2.5)

rm(cols, brk, paletteLength, col_anno, ann_colors, pr_res_summ)

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Violin plots of activity scores
#-------------------------------------------------------------------------------

if(identical(rownames(pr_res), colData(rld)$SampleID))
{
  pr_res_dat = cbind(pr_res, colData(rld)[,2:3])
  pr_res_dat = melt(as.data.frame(pr_res_dat))
  
  p1 = ggplot(pr_res_dat, aes(x = Time, y = value, fill = Treatment)) + 
    geom_boxplot(lwd = 0.1) +
    labs(y = "Signalling score", x = "") +
    scale_fill_manual(values = c("#fb8072", "#b3de69", "#fdb462", "#80b1d3")) +
    #geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic(base_size = 8) + facet_wrap(~variable , scales = "free", ncol = 4) + theme(strip.background = element_rect(size = 0.1))
  
  ggsave(filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_violin_byTreatment.pdf"), plot = p1, width = 7.5, height = 6)
  
  p2 = ggplot(pr_res_dat, aes(x = Treatment, y = value, fill = Time)) + 
    geom_boxplot(lwd = 0.1) +
    labs(y = "Signalling score", x = "") +
    scale_fill_manual(values = c("#cccccc", "#969696", "#636363", "#252525")) +
    #geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic(base_size = 8) + facet_wrap(~variable , scales = "free", ncol = 4) + theme(strip.background = element_rect(size = 0.1))
  
  ggsave(filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_violin_byTime.pdf"), plot = p2, width = 7.5, height = 6)
}

rm(p1, p2, pr_res_dat)

#-------------------------------------------------------------------------------
# TF regulons from Dorothea
#-------------------------------------------------------------------------------

# Extracting literature curatded TF regulons from Dorothea
data(dorothea_hs, package = "dorothea")

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

length(unique(regulons$tf))
#[1] 292

# Computing TF activities using viver
tf_activities <- run_viper(expr, regulons, 
                           options =  list(method = "scale", minsize = 25, eset.filter = FALSE, 
                                           cores = 1, verbose = FALSE, nes=T))
nrow(tf_activities)
#[1] 181

write.csv(data.frame(Samples = rownames(tf_activities), tf_activities, stringsAsFactors = F),
          paste0(path, "analysis/Exploratory/tfactivity_all_results.csv"), 
          quote = FALSE, row.names = F)

# Selecting the top 25% most variable TFs
v = apply(tf_activities, 1, mad)
tf_activities = t(tf_activities[v > quantile(v, 0.75),])

# Summarization
tf_activities.summ = split(data.frame(tf_activities), 
                    sapply(strsplit(rownames(tf_activities), "-Rep-"), function(x) x[1]))

tf_activities.summ = t(sapply(tf_activities.summ, function(x){
  apply(x, 2, median)
}))

tf_activities.summ = tf_activities.summ[c("Mock-3h", "Mock-6h", "Mock-12h", "Mock-24h",
                                          "IL22-3h", "IL22-6h", "IL22-12h", "IL22-24h",
                                          "IFNL-3h", "IFNL-6h", "IFNL-12h", "IFNL-24h",
                                          "IL22_IFNL-3h", "IL22_IFNL-6h", "IL22_IFNL-12h", "IL22_IFNL-24h"),]
ncol(tf_activities.summ)
#[1] 45

# Annotations for plotting
paletteLength = 100
cols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(paletteLength)
brk = c(seq(min(tf_activities.summ), 0, length.out=ceiling(paletteLength/2) + 1),
        seq(max(tf_activities.summ)/paletteLength,  max(tf_activities.summ), length.out=floor(paletteLength/2)))

row_anno <- data.frame(do.call("rbind", strsplit(rownames(tf_activities.summ), "-")), row.names = rownames(tf_activities.summ))
colnames(row_anno) = c("Treatment", "Time")

ann_colors = list(
  Treatment =  c(Mock = "#fb8072", IL22 = "#b3de69", IFNL = "#fdb462", IL22_IFNL = "#80b1d3"),
  Time = c(`3h` = "#cccccc", `6h` = "#969696", `12h` = "#636363", `24h` = "#252525")
)

pheatmap(tf_activities.summ,
         clustering_method = "ward.D2",
         cluster_rows = F,
         gaps_row = c(4,8,12),
         cutree_cols = 3,
         treeheight_col = 15,
         fontsize = 6,
         border_color = "white", 
         col = cols,
         breaks = brk,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         filename = paste0(path, "analysis/Exploratory/TFactivity_Dorothea_Viper.pdf"), width = 5.8, height = 2.5)

rm(cols, brk, v, paletteLength, row_anno, ann_colors, dorothea_hs, regulons, tf_activities.summ)

#-------------------------------------------------------------------------------
# TF - TF interactions
#-------------------------------------------------------------------------------

# #int = import_PathwayExtra_Interactions(filter_databases = c("BioGRID","STRING"), select_organism = 9606)
# int = import_PathwayExtra_Interactions(select_organism = 9606)
# int = int[int$source_genesymbol %in% rownames(tf_activities) & int$target_genesymbol %in% rownames(tf_activities),]
# int = int[,c(1:10, 13)]
