library(ggplot2)
library(tidyr)
library(dplyr)
library(qvalue)
library(forcats)

outname<-"Violinplot_H1_Odds_ratio_NFR_Dnase_diff_0.2_NDR_Dnase_conserved_0.8_120_180_Cell_Differentation_Reprogramming"
TF_NFR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_H1_Dnase_NFR_diff_0.2.txt",header = FALSE,sep="\t")
TF_NDR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_H1_Dnase_NDR_120_180_conserved_0.8.txt",header = FALSE,sep="\t")
All_NFR_len <-read.table("./length_of_NFR_NDR/Dnase_H1_NFR_diff_0.2.txt",header = FALSE)
All_NDR_len <-read.table("./length_of_NFR_NDR/Dnase_H1_NDR_120_180_conserved_0.8.txt",header = FALSE )

#outname<-"Violinplot_Odds_ratio_NFR_0_NDR_Dnase_120_180_Cell_Differentation_Reprogramming"
#TF_NFR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_NFR_0.txt",header = FALSE,sep="\t")
#TF_NDR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_Dnase_NDR_120_180.txt",header = FALSE,sep="\t")
#All_NFR_len <-read.table("./length_of_NFR_NDR/NFR_0_146_148.txt",header = FALSE)
#All_NDR_len <-read.table("./length_of_NFR_NDR/Dnase_NDR_120_180.txt",header = FALSE )

all_enrichment_fold <- data.frame(TF_name = character(), cell_line = character(),TF_cell_line = character(), odds_ratio = double(), overlap_bp_NFR = double(),overlap_bp_NDR = double(),stringsAsFactors = FALSE)

for (i in 1:nrow(TF_NFR)) {
  cell_line <- TF_NFR$V1[i] %>% as.character()
  TF_name <-TF_NFR$V2[i] %>% as.character()
  NFR_overlap_bp <- TF_NFR$V3[i] %>% as.numeric()
  NDR_overlap_bp <- filter(TF_NDR,V1==cell_line & V2==TF_name) %>% select(V3) %>% as.numeric()
  
  NFR_len <- filter(All_NFR_len,V1==cell_line ) %>% select(V2) %>% as.numeric()
  NDR_len <- filter(All_NDR_len,V1==cell_line ) %>% select(V2) %>% as.numeric()
  
  conti_table <- data.frame("NFR" = c(NFR_overlap_bp, (NFR_len-NFR_overlap_bp)), 
                            "NDR" = c(NDR_overlap_bp, (NDR_len-NDR_overlap_bp)), row.names = c("binding", "Non_binding"))

  odds_ratio <- (NFR_overlap_bp/(NFR_len-NFR_overlap_bp))/(NDR_overlap_bp/(NDR_len-NDR_overlap_bp))

  all_enrichment_fold <- rbind(all_enrichment_fold,
                               c(TF_name, cell_line, paste(TF_name,cell_line,sep="_"),odds_ratio,NFR_overlap_bp, NDR_overlap_bp), stringsAsFactors =FALSE)
  colnames(all_enrichment_fold) <- c("TF_name", "cell_line", "TF_cell_line", "odds_ratio", "overlap_bp_NFR","overlap_bp_NDR")
  
}



PTF1 <- read.table("./PTF_Cell_differentiation.txt",header = TRUE, sep = ",",stringsAsFactors = FALSE)
PTF1 <- PTF1$PTF_name

PTF2 <- read.table("./PTF_Cell_reprogramming.txt",header = TRUE, sep = ",",stringsAsFactors = FALSE)
PTF2 <- PTF2$PTF_name

data <- all_enrichment_fold %>% na.omit() %>%
  arrange(desc(as.numeric(overlap_bp_NFR))) 


data <- data %>% arrange(desc(as.numeric(odds_ratio))) %>%
  mutate(TF_cell_line = fct_reorder(TF_cell_line, desc(as.numeric(odds_ratio)))) #%>% filter(overlap_bp >100000)

data$type <- "Other TF"
data$type[data$TF_name %in% PTF1] <- "Differentiation"
data$type[data$TF_name %in% PTF2] <- "Reprogramming"
data$type <- factor(data$type , levels=c("Differentiation", "Reprogramming", "Other TF"))

data <- data %>% filter(overlap_bp_NDR >0)
data$odds_ratio <- as.numeric(data$odds_ratio)

ggplot(data, aes(x=type, y=odds_ratio, fill=type)) +
  geom_violin(outlier.shape = NA ) + 
  geom_jitter(shape=16, position=position_jitter(0.1), size=0.5, alpha=0.3) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  ylab("Enrichment score")  +

  scale_y_continuous(trans='log10',limits = c(min(as.numeric(data$odds_ratio)), 20)) + 
  theme(panel.border = element_rect(size =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title=element_text(size=16, vjust=3)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 14,face="bold"),
    axis.text.y = element_text(size = 14,face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 13,face="bold", angle = 0),
    axis.ticks.x = element_blank())
ggsave(paste("./", outname,".png",sep=""),width=6,height=3,units="in",dpi=300)
data$type

#Mann Whitney U test
wilcox.test(odds_ratio~type,
            data = data[data$type %in% c("Reprogramming","Other TF"),],
           exact = FALSE)