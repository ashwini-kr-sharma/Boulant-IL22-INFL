options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

bamdir  = as.character(args[1])
gtffile = as.character(args[2])
outdir  = as.character(args[3])

print(bamdir)
print(gtffile)
print(outdir)

library(Rsubread)

counts = featureCounts(files = list.files(bamdir, full.names = T, pattern=".*bam$"),
                       annot.ext = gtffile,
                       isGTFAnnotationFile = T,
                       GTF.attrType.extra = c("gene_name", "gene_biotype"),
                       
                       isPairedEnd = T,
                       checkFragLength = T,
                       requireBothEndsMapped = T,
                       countChimericFragments = F,
                       
                       nthreads = 20)

saveRDS(counts, file = paste0(outdir, "T84_IL22_INFL_counts.RDS"))
