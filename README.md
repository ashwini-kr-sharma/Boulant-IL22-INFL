## T84 cell treated with IL22 and IFNL

- This repository hosts the source code used to analyze the data from the paper - [IL-22 promotes cell proliferation to combat virus infection in human intestinal epithelial cells]()
- The companion website for this paper is hosted [here](https://ashwini-kr-sharma.github.io/T84-IL22-IFNL/)

### Scripts
All scripts are located in `/src/` directory

0. `/src/00_submitJobs.sh` - This is the `bash` script to run analysis `1,2,3,4` in the cluster environment
1. `/src/01_fastQC.sh` - Script to perform `fastQC` analysis on the raw `fastq` files
2. `/src/02_alignmentRsubread.R` - Script to perform `read alignment` analysis on the raw `fastq` files using `Rsubread`
3. `/src/03_countsRsubread.R` - Script to `count reads` overlapping genes
4. `/src/04_multiQC.sh` - Script to perform `QC analysis` of raw reads and alignment statistics
5. `/src/05_postProcessing.R` - Script to preprocess and filter the raw count data
6. `/src/06_exploratoryAnalysis.R` - Script to perform various `Exploratory analysis
7. `/src/07_DSeq2Analysis.R` - Script to perform `Differential gene expression analysis (DGE)` using DSeq2
8. `/src/07_DSeq2_Rmarkdown/` - Script to `visualize` the results of Dseq2 as interactive data tables
9. `/src/08_differentialGeneExpressionAnalysis.R` - Script to `analyze DGE results` focusing on specific genes of interest
10. `/src/09_enrichmentAnalysis.R` - Script to perform `Enrichment analyses`
11. `/src/10_plotGeneModulesExpression.R` - Script to plot the expression of `JAK-STAT`, `IL22` and `IFN-L` signaling genes
