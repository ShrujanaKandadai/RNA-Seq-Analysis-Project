# RNA-Seq-Analysis-Project
Repository for RNA-Seq Analysis of a dataset obtained from Patients infected with Leishmania braziliensis. This dataset was analyzed to understand and recreate the findings created in the original paper referencing this dataset. 

**Analysis status:** [active/completed]  
**Last updated:** [2025-07-11]

🌐 [Live site](https://shrujanakandadai.github.io/RNA-Seq-Analysis-Project/)

![report preview](RNA-Seq-Analysis-Project_files/figure-html/step-1.png)

## Project Overview

In this analysis, I have worked on recreating the analysis done in the the study by [Amorim et al., 2019](https://doi.org/10.1126/scitranslmed.aax4204) Here, I have investigated transcriptomic differences between **HS (healthy skin) vs. CL lesions** using publicly available RNA-seq data from **GSE127831** from the study mentioned. The dataset comprises 28 samples,from human skin biopsies of healthy individuals(n=7) and from lesions of patients infected with Leishmania braziliensis (n=21).This analysis revealed that distinct transcriptional signatures differentiate HS and CL Lesions and that within the CL, gene expression profiles correlate with treatment outcomes (cure vs. failure). Novel additions in my analysis include identification of 5 novel target genes that are highly upregulated in CL compared to HS, and demonstration of cytotoxic genes (e.g., GZMB, PRF1) that are highly expressed in the cases of CL where treatment failure occurred.

---

## Dataset Summary
>GSE 127831 - The dataset comprises 28 samples,from human skin biopsies of healthy individuals(n=7) and from lesions of       patients infected with Leishmania braziliensis (n=21).
---
### Data Retrieval
>Instead of starting from raw FASTQ files, we retrieved **pre-aligned gene-level expression data** from the **ARCHS4** database.
---
### Key Analytical Steps
1. Filtering and Normalization : [As per the study - Normalized using the TMM method in EdgeR and genes with < 1 count per million (CPM) in 7 of the samples were filtered out.] 
2. Multivariate data analysis: [PCA Plot]
3. Differential Expression: [uNormalized, filtered data were variance-stabilized using the VOOM function in limma and differentially expressed genes were identified with linear modeling using limma (FDR ≤ 0.01; absolute logFC ≥ 1), after correcting for multiple testing using Benjamini-Hochberg]. 
4. GO Analysis: [GSEABase, Biobase, gprofiler2]
5. GSEA: [clusterProfiler and msigdbr]
---
## Key Findings
1. Distinct transcriptional signatures differentiate HS and CL lesions.
2. Within CL, gene expression profiles correlate with treatment outcomes (cure vs. failure).
3. Cytotoxic genes (e.g., GZMB, PRF1) are highly expressed in failure cases.
4. Pathway analyses reveal immune activation and suppression patterns relevant to disease state. 
---
## References
Amorim, C. F., Novais, F. O., Nguyen, B. T., Misic, A. M., Carvalho, L. P., Carvalho, E. M., Beiting, D. P., & Scott, P. (2019). Variable gene expression and parasite load predict treatment outcome in cutaneous leishmaniasis. Science translational medicine, 11(519), eaax4204. https://doi.org/10.1126/scitranslmed.aax4204

In this study, transcriptome profiling was carried out on biopsies obtained from healthy individuals and patients infected with Leishmania braziliensis resulting in chronic lesions. Apart from highlighting highly differential gene expression between healthy vs CL samples, this study also set forth to determine whether genes whose expression is highly variable correlated with treatment outcome. Amongst the most variable genes were components of the cytolytic pathway, the expression of which appeared to be driven by parasite load in the skin thus revealing that treatment failure can be directly linked to the cytolytic pathway activated during infection.

I have worked on recreating the analysis done in this study, and have further added a novel component, by highlighting the 5 most differentially regulated genes as novel targets for treatment of CL. 


