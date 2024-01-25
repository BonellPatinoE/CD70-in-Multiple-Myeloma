# CD70-in-Multiple-Myeloma
RNAseq data analysis about CD70 expression in primary samples from coMMpass (MMRF)

Bioinformatic analysis performed at Wiita Lab at University of California, San Francisco (UCSF), on RNAseq data showed previously higher CD70 expression in relapsed myeloma patients. Here we analyzed CD70 expression in RNAseq data from coMMpass dataset (MMRF)in CD138+ cells from primary samples.

This analysis was performed based on cytogenetic translocations and genomic rearrangements that defined different risks in multiple myeloma to determine if there were any associations with CD70 expression and different cytogenetic profiles. 

This analysis was performed using DESeq2 (v1.34.0) package in R (v4.1.2) as well, using 'apeglm' for LFC shrinkage and “gprofiler” (v0.2.2) for for gene list functional enrichment analysis and namespace conversion. 
