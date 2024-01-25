# CD70-in-Multiple-Myeloma
RNAseq data analysis about CD70 expression in primary samples from coMMpass (MMRF) IA19.

Bioinformatic analysis performed at Wiita Lab at the University of California, San Francisco (UCSF), on RNAseq data showed previously higher CD70 expression in relapsed myeloma patients. Here we analyzed CD70 expression in RNAseq data from coMMpass dataset (IA19) in CD138+ cells from primary samples.

This analysis was performed based on CD70 expression based on well-characterized tumor genotypes and showed that several high-risk myeloma karyotypic subtypes, specifically those with t(4;14), gain of Chr1q, and loss of Chr1p, showed significantly upregulated CD70 expression compared to all others in CoMMpass dataset. 

This analysis was performed using DESeq2 (v1.34.0) package in R (v4.1.2) as well, using 'apeglm' for LFC shrinkage and “gprofiler” (v0.2.2) for for gene list functional enrichment analysis and namespace conversion. 
