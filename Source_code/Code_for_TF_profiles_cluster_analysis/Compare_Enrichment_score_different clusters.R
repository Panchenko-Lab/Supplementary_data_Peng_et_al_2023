library(ggplot2)
library(tidyr)
library(dplyr)
outname <-"Nucleosome_binding_in_solution_validation"

Odds_ratio <- read.table("Enrichment_score_NFR_0_Dnase_diff_0.2_NDR_Dnase_conserved_0.8_120_180_filter_PTF_non_PTF.csv",header = TRUE, sep=",",stringsAsFactors = FALSE)

cluster1 <- c("ELK1","KLF4","SPI1","BTF3","HMGA1","TFAP2A","TFAP2B","TFAP2C",
              "HNF1A","GATA4","RBPJ","ASCL","E12A","FOXA1","OCT4","POU5F1","BRN2","HMGN5","SOX2","HMGN1")
cluster2 <- c("NFACTC1","ESRRG","SOX5","ZNF250","NHLH2","CRX","CEBPA")
cluster3 <- c("CREM","IRF3","ELF2","SPDEF","MYOG","T","MYC","BRN5","PKNOX","TBX20")

## Cluster 1 strong binders to both naked DNA and nucleosomal DNA 
## Cluster 2 weak binders to both naked DNA and nucleosomal DNA (only one TF)
## Cluster 3: strong binders to naked DNA but weak binders to nucleosomal DNA. 
TF_microarray_cluster1 <- Odds_ratio %>% filter(TF_name %in% cluster1) 
TF_microarray_cluster2 <- Odds_ratio %>% filter(TF_name %in% cluster2) 
TF_microarray_cluster3 <- Odds_ratio %>% filter(TF_name %in% cluster3) 

TF_microarray_cluster1$cluster <- "Cluster1 and 2"
TF_microarray_cluster2$cluster <- "Cluster1 and 2"
TF_microarray_cluster3$cluster <- "Cluster3"

TF_microarray_cluster <- rbind(TF_microarray_cluster1,TF_microarray_cluster2,TF_microarray_cluster3)


ggplot(TF_microarray_cluster, aes(x=cluster, y=odds_ratio, fill=cluster)) +
  geom_violin(outlier.shape = NA ) + 
  geom_jitter(shape=16, position=position_jitter(0.1), size=2.5, alpha=0.4) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  ylab("Enrichment score")  +
  scale_y_continuous(limits = c(min(as.numeric(TF_microarray_cluster$odds_ratio)), max(as.numeric(TF_microarray_cluster$odds_ratio)))) + 
  theme(panel.border = element_rect(size =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title=element_text(size=16, vjust=3)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 16,face="bold"),
    axis.text.y = element_text(size = 14,face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16,face="bold"),
    axis.ticks.x = element_blank())
ggsave(paste(outname,".png",sep=""),width=6,height=3,units="in",dpi=300)

#Mann Whitney U test
wilcox.test(odds_ratio~cluster,
            data = TF_microarray_cluster[TF_microarray_cluster$cluster %in% c("Cluster1 and 2","Cluster3"),],
            exact = FALSE)