## T84 cells treated with IFNL and IL22

### Experimental design

T84 cells were treated with the cytokines, type III interferon lambda (IFN-L), interleukin 22 (IL22) and in combination at time points - 3h, 6h, 12h and 24h. Each time point and treatment were replicated three times. RNA sequencing was performed in a total of 48 samples consiting of untreated and treated samples.

![experimental design](/IL22_IFNL_Expdesign.png)

### Methodology

Raw RNA sequencing reads (fastq) were aligned to the ensembl human GRCh38 genome reference using [Rsubread(2.2.6)](https://doi.org/10.1093/nar/gkz114) with default settings. Read summariztion was done using [featureCounts](https://doi.org/10.1093/bioinformatics/btt656). Various quality metrics of the the raw reads and alignment statistics were analysed using [MultiQC](https://doi.org/10.1093/bioinformatics/btw354). Differential gene expression analysis was performed uing DSeq2, rlog transformed data was used for multi dimensional scaling and clustering analyses. Signalling programs were quantified using [PROGENy(1.10.0)](https://doi.org/10.1038/s41467-017-02391-6). Transcription factor activities were computed using [DoRothEA(1.0.1)](https://doi.org/10.1101/gr.240663.118) and [VIPER](https://doi.org/10.1038/ng.3593). Enrichment analysis on the most differentially expressed genes (-1< logFC >+1 and adjusted p value < 0.05) was performed using [enrichR(3.0)](https://doi.org/10.1093/nar/gkw377).

### Data download
Raw and processsed data from this study can be downloaded below. The source code for this study is available [here](https://github.com/ashwini-kr-sharma/Boulant-IL22-IFNL)

- [Raw fastq files](https://www.ncbi.nlm.nih.gov/gds)
- Raw counts - [RDS](/data/T84_IL22_INFL_filtered_counts.RDS), [csv](/data/T84_IL22_INFL_filtered_counts.csv)
- rlog transformed data - [RDS](/data/rlogTransformation.RDS)
- [Quality control - MultiQC](/results/MultiQC/multiqc_report.html)
- Differentially expresed genes (All) - [RDS](/data/diffExpGenes.RDS), [xls](/data/DGEtables.xls)
- Differentially expressed genes - IL22
  - [IL22 vs Mock (3hr)](/src/07_DSeq2_Rmarkdown/IL22_3hr_vs_Mock_3hr.html)
  - [IL22 vs Mock (6hr)](/src/07_DSeq2_Rmarkdown/IL22_6hr_vs_Mock_6hr.html)
  - [IL22 vs Mock (12hr)](/src/07_DSeq2_Rmarkdown/IL22_12hr_vs_Mock_12hr.html)
  - [IL22 vs Mock (24hr)](/src/07_DSeq2_Rmarkdown/IL22_24hr_vs_Mock_24hr.html)
- Differentially expressed genes - INFL
  - [INFL vs Mock (3hr)](/src/07_DSeq2_Rmarkdown/IFNL_3hr_vs_Mock_3hr.html)
  - [INFL vs Mock (6hr)](/src/07_DSeq2_Rmarkdown/IFNL_6hr_vs_Mock_6hr.html)
  - [INFL vs Mock (12hr)](/src/07_DSeq2_Rmarkdown/IFNL_12hr_vs_Mock_12hr.html)
  - [INFL vs Mock (24hr)](/src/07_DSeq2_Rmarkdown/IFNL_24hr_vs_Mock_24hr.html)
- Differentially expressed genes - IL22 + INFL
  - [IL22 + INFL vs Mock (3hr)](/src/07_DSeq2_Rmarkdown/IL22_IFNL_3hr_vs_Mock_3hr.html)
  - [IL22 + INFL vs Mock (6hr)](/src/07_DSeq2_Rmarkdown/IL22_IFNL_6hr_vs_Mock_6hr.html)
  - [IL22 + INFL vs Mock (12hr)](/src/07_DSeq2_Rmarkdown/IL22_IFNL_12hr_vs_Mock_12hr.html)
  - [IL22 + INFL vs Mock (24hr)](/src/07_DSeq2_Rmarkdown/IL22_IFNL_24hr_vs_Mock_24hr.html)
- [log2 fold changes](/src/07_DSeq2_Rmarkdown/log2_fold_change.html) [(RDS)](/data/diffExpLogFCmatrix.RDS)
- [Signalling programs](/data/progeny_all_results.csv)
- [Transcription factor activities](data/tfactivity_all_results.csv)
