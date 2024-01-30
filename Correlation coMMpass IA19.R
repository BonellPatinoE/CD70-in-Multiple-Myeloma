library(tidyverse)
library(dslabs)
library(dplyr)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(reshape2)
library(gprofiler2)
library(ggpubr)
library(naniar)

## Loading the coMMpass data set ##

coMMpass<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
coMMpass<-data.frame(coMMpass)

colnames(coMMpass) 
coMMpass1<-coMMpass

rownames(coMMpass1) <- coMMpass1[,1]
coMMpass1[,1] <- NULL
head(coMMpass1)

## gconvert to get gene names ##

coMMpass2<-gconvert(row.names(coMMpass1), organism = "hsapiens", target="HGNC", mthreshold = 1, filter_na = FALSE)
coMMpass2

coMMpass3<-coMMpass1
rownames(coMMpass3)<-NULL
coMMpass3$Gene_name<-coMMpass2$name
nrow(coMMpass3)

coMMpass4<-coMMpass3 %>% distinct(Gene_name, .keep_all = TRUE)

coMMpass4<-coMMpass4%>%filter(!is.na(Gene_name))

row.names(coMMpass4)<-coMMpass4$Gene_name
coMMpass5<-coMMpass4  
coMMpass5
coMMpass5$Gene_name<-NULL

coMMpass6<-t(coMMpass5)
coMMpass6<-data.frame(coMMpass6)

class(coMMpass6)
colnames(coMMpass6)
coMMpass6$CD70

coMMpass7<-log2(coMMpass6+1)

coMMpass6%>%select(CD70, LEF1)


pdf("Transcription factors and CD70 Correlation plots.pdf")

## Correlation plots ##

ggscatter(coMMpass6, x = "CD70", y = "TFAP2A", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "TFAP2A log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD70 vs TFAP2A (Pearson R value)")

ggscatter(coMMpass6, x = "CD70", y = "TFAP2A", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "TFAP2A log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD70 vs TFAP2A (Spearman R value)")

ggscatter(coMMpass6, x = "CD70", y = "PAX5", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "PAX5 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs PAX5 (Pearson R value)")

ggscatter(coMMpass6, x = "CD70", y = "PAX5", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "PAX5 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs PAX5 (Spearman R value)")

dev.off()

ggscatter(coMMpass6, x = "CD70", y = "ATF3", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "ATF3 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs ATF3 (Spearman R value)")

ggscatter(coMMpass6, x = "CD70", y = "ATF3", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "ATF3 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs ATF3 (Pearson R value)")


ggscatter(coMMpass6, x = "CD70", y = "LEF1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "LEF1 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs LEF1 (Pearson R value)")


ggscatter(coMMpass6, x = "CD70", y = "LEF1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "LEF1 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs LEF1 (Spearman R value)")


ggscatter(coMMpass6, x = "CD70", y = "JUN", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "JUN log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs JUN (Spearman R value)")

ggscatter(coMMpass6, x = "CD70", y = "JUN", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "JUN log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs JUN (Pearson R value)")

ggscatter(coMMpass6, x = "CD70", y = "FOXA1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "FOXA1 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs FOXA1 (Pearson R value)")


ggscatter(coMMpass6, x = "CD70", y = "FOXA1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "FOXA1 log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs FOXA1 (Spearman R value)")

ggscatter(coMMpass6, x = "CD70", y = "MYB", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD70 log2", ylab = "MYB log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs MYB (Spearman R value)")

ggscatter(coMMpass6, x = "CD70", y = "MYB", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD70 log2", ylab = "MYB log2")+
  xscale("log2", .format = TRUE)+
  yscale("log2", .format = TRUE)+
  ggtitle("CD70 vs MYB (Pearson R value)")


## References @@

?ggscatter

citation("ggpubr")
