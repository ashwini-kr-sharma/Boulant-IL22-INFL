## T84 cells treated with INFL and IL22

### Experimental design
![experimental design](/IL22_IFNL_Expdesign.png)

### Data download
Raw and processsed data from this study can be downloaded below. The source code for this study is available [here](https://github.com/ashwini-kr-sharma/Boulant-IL22-INFL)

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
- [log2 fold changes](/data/DGE/log2_fold_change.html) [(RDS)](/data/diffExpLogFCmatrix.RDS)
- [Signalling programs](/data/progeny_all_results.csv)
- [Transcription factor activities](data/tfactivity_all_results.csv)
