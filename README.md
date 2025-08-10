# CD70-in-Multiple-Myeloma

RNA-seq analysis of CD70 expression in primary myeloma samples from the CoMMpass (MMRF) IA19 dataset, as reported in Blood Journal: Targeting high-risk multiple myeloma genotypes with optimized anti-CD70 CAR-T cells (https://doi.org/10.1182/blood.2024025536).

Bioinformatic analysis conducted at the Wiita Lab, University of California, San Francisco (UCSF), USA. We examined CD70 expression in RNA-seq data from the CoMMpass IA19 dataset, focusing on CD138⁺ cells from primary patient samples, demonstrating elevated CD70 expression in relapsed multiple myeloma. 

Analysis by well-characterized tumor genotypes revealed that several high-risk cytogenetic subtypes—specifically t(4;14), gain of 1q, and loss of 1p—exhibited significantly upregulated CD70 expression compared with all other cases in the dataset.

Differential expression analysis was performed using the DESeq2 package (v1.34.0) in R (v4.1.2), with log₂ fold-change shrinkage via apeglm. Functional enrichment and namespace conversion were conducted using the gprofiler2 package (v0.2.2).
