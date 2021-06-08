## T84 cell treated with IL22 and IFNL

- This repository hosts the souce code used to analyze the data from the paper - [XXX]()
- The companion website for this paper is hosted [here](https://ashwini-kr-sharma.github.io/Boulant-IL22-INFL/)

### Scripts
All scripts are located in `/src/` directory

- `/src/00_submitJobs.sh` - This is the `bash` script to run the `fastQC`, `read alignment`, `gene counting` and `MultiQC` analysis in the cluster environment
- `/src/01_fastQC.sh` - Script to perform `fastQC` analysis on the raw `fastq` files
- `/src/02_alignmentRsubread.R` - Script to perform `read alignment` analysis on the raw `fastq` files using `Rsubread`
- `/src/03_countsRsubread.R` - Script to `count reads` overlapping genes
- `/src/04_multiQC.sh` - Script to perform `QC analysis` of raw reads and alignment statistics
- `/src/05_postProcessing.R` - Script to preprocess and filter the raw count data
- `/src/06_exploratoryAnalysis.R` - Script to perform various `Exploratory analysis`
- `/src/07_DSeq2Analysis.R` - Script to perform `Differential gene expression analysis (DGE)` using DSeq2
- `/src/07_DSeq2_Rmarkdown/` - Script to `visualize` the results of Dseq2 as interactive data tables
- `/src/08_differentialGeneExpressionAnalysis.R` - Script to `analyze DGE results` focusing on specific genes of interest
- `/src/09_enrichmentAnalysis.R` - Script to perform `Enrichment analyses`
- `/src/10_plotGeneModulesExpression.R` - Script to plot the expression of `JAK-STAT`, `IL22` and `IFN-L` signalling genes
