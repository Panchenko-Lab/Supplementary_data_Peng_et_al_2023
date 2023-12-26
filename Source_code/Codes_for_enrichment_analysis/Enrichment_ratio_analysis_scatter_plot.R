library(ggplot2)
library(tidyr)
library(dplyr)
library(qvalue)
library(forcats)

outname<-"Scatterplot_H1_Odds_ratio_NFR_Dnase_diff_0.2_NDR_Dnase_conserved_0.8_120_180_Cell_Differentation_Reprogramming"
TF_NFR <- read.table("../TF_motif_counts_NFR_NDR/anchor_H1_embryonic/All_TF_Dnase_NFR_diff_0.2.txt",header = FALSE,sep="\t")
TF_NDR <- read.table("../TF_motif_counts_NFR_NDR/anchor_H1_embryonic/All_TF_Dnase_NDR_120_180_conserved_0.8.txt",header = FALSE,sep="\t")
All_NFR_len <-read.table("../length_of_NFR_NDR/Dnase_NFR_diff_0.2.txt",header = FALSE)
All_NDR_len <-read.table("../length_of_NFR_NDR/Dnase_NDR_120_180_conserved_0.8.txt",header = FALSE )

#outname<-"Scatterplot_Odds_ratio_NFR_0_NDR_Dnase_120_180_Cell_Differentation_Reprogramming"
#TF_NFR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_NFR_0.txt",header = FALSE,sep="\t")
#TF_NDR <- read.table("./TF_motif_counts_NFR_NDR/All_TF_Dnase_NDR_120_180.txt",header = FALSE,sep="\t")
#All_NFR_len <-read.table("./length_of_NFR_NDR/NFR_0_146_148.txt",header = FALSE)
#All_NDR_len <-read.table("./length_of_NFR_NDR/Dnase_NDR_120_180.txt",header = FALSE )

all_enrichment_fold <- data.frame(TF_name = character(), cell_line = character(),TF_cell_line = character(), odds_ratio = double(), pvalue = double(), overlap_bp_NFR = double(),overlap_bp_NDR = double(),stringsAsFactors = FALSE)

for (i in 1:nrow(TF_NFR)) {
  cell_line <- TF_NFR$V1[i] %>% as.character()
  TF_name <-TF_NFR$V2[i] %>% as.character()
  NFR_overlap_bp <- TF_NFR$V3[i] %>% as.numeric()
  NDR_overlap_bp <- filter(TF_NDR,V1==cell_line & V2==TF_name) %>% select(V3) %>% as.numeric()
  
  NFR_len <- filter(All_NFR_len,V1==cell_line ) %>% select(V2) %>% as.numeric()
  NDR_len <- filter(All_NDR_len,V1==cell_line ) %>% select(V2) %>% as.numeric()
  
  conti_table <- data.frame("NFR" = c(NFR_overlap_bp, (NFR_len-NFR_overlap_bp)), 
                            "NDR" = c(NDR_overlap_bp, (NDR_len-NDR_overlap_bp)), row.names = c("binding", "Non_binding"))
  
  fisher_test <-fisher.test(conti_table,alternative = "greater")
  p_value<-fisher_test$p.value

  odds_ratio <- (NFR_overlap_bp/(NFR_len-NFR_overlap_bp))/(NDR_overlap_bp/(NDR_len-NDR_overlap_bp))

  all_enrichment_fold <- rbind(all_enrichment_fold,
                               c(TF_name, cell_line, paste(TF_name,cell_line,sep="_"),
                                 odds_ratio,p_value,NFR_overlap_bp, NDR_overlap_bp), stringsAsFactors =FALSE)
  colnames(all_enrichment_fold) <- c("TF_name", "cell_line", "TF_cell_line", "odds_ratio", "p_value", "overlap_bp_NFR","overlap_bp_NDR")
  
}

## get q-values
pvalues <- as.numeric(all_enrichment_fold$p_value)
qobj <- qvalue(p = pvalues)
qvalues <- qobj$qvalues
all_enrichment_fold <- cbind(all_enrichment_fold, qvalues)


PTF_Cell_differentation <- c("FOXA1", "FOXA2", "FOXA3","GATA1", "GATA2", "GATA3", "GATA4","CEBPA","CEBPB","NEUROD1","SPI1")
PTF_Reprogramming  <- c("NFYA", "NFYB","NFYC","ESRRB","POU5F1", "KLF4")


data <- all_enrichment_fold %>% na.omit() %>%
  arrange(desc(as.numeric(overlap_bp_NFR))) 

data <- data %>% arrange(TF_cell_line) %>%
  mutate(TF_cell_line = fct_reorder(TF_cell_line, TF_cell_line)) #%>% filter(overlap_bp >100000)

data <- data %>% filter(overlap_bp_NDR >0) %>% filter(odds_ratio >0)

ggplot(data[data$odds_ratio!=0,], aes(x=TF_cell_line, y=as.numeric(odds_ratio))) +
  geom_point(aes(colour = qvalues), 
      size=ifelse(data$TF_name %in% c(PTF_Cell_differentation), 3, 
                       ifelse(data$TF_name %in% c(PTF_Reprogramming), 3, 2)),
      shape=ifelse(data$TF_name %in% c(PTF_Cell_differentation), 15, 
                         ifelse(data$TF_name %in% c(PTF_Reprogramming), 17, 20)),
      alpha=ifelse(data$TF_name %in% c(PTF_Cell_differentation), 1, 
                          ifelse(data$TF_name %in% c(PTF_Reprogramming), 1, 0.3)),) + 
  
  scale_colour_gradient2(
       low = "darkred",
       mid = "yellow",
       high = "darkgreen",
       midpoint = 0.5,
       space = "Lab",
       na.value = "grey50",
       guide = "colourbar",
       aesthetics = "colour"
     ) + 
  theme_bw() +
  scale_fill_discrete(name="q-values") +
  theme(legend.position="right",legend.title = element_text(colour="black", size=18, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold") ) +
  xlab("Transcription Factor") +
  ylab("Enrichment score")  +
  scale_y_continuous(trans='log10',limits = c(min(as.numeric(data$odds_ratio))+0.01, 15)) + 
  #ylim(0,30) +
  #theme(axis.title.x=element_blank()) + 
  theme(panel.border = element_rect(size =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(axis.title.y=element_text(angle=90, vjust=0)) +
  theme(plot.title=element_text(size=16, vjust=3)) +
  theme(axis.ticks.x = element_text(hjust = -5))+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 20,face="bold"),
    axis.text.y = element_text(size = 15,face="bold"),
    axis.title.x = element_text(size = 20,face="bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  theme(axis.title.x=element_text(margin=margin(0,0,0,0)))
ggsave(paste("./", outname,".png",sep=""),width=8,height=4,units="in",dpi=300)
write.table(data,paste("./", outname,"2.csv",sep=""),quote = FALSE,row.names = FALSE,sep = ",")