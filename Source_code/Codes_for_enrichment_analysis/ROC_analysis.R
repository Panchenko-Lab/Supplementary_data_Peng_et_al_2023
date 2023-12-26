library(ggplot2)
library(plotROC)
library(PRROC)
library(mltools)
library(tidyr)
library(dplyr)

#file_name <-"NFR_0_Dnase_NDR_120_180_All_PTF_expressed"
file_name <-"NFR_Dnase_NFR_0_diff_0.2_NDR_conserved_0.8_120_180_All_PTF_Cell_differentiation"
#file_name <-"NFR_Enhancer_NFR_diff_0.2_NDR_conserved_0.8_All_PTF_expresed"

### Enrichment_scores of all TFs
Enrichment_ratio1 <- read.table("./Enrichment_scores/Enrichment_score_NFR_0_Dnase_diff_0.2_NDR_Dnase_conserved_0.8_120_180.csv",header = TRUE,sep=",")
Enrichment_ratio2 <- read.table("./Enrichment_scores/Enrichment_score_NFR_0_Dnase_diff_0.2_NDR_Dnase_conserved_0.8_120_180_expressed.csv",header = TRUE,sep=",")
Enrichment_ratio3 <- read.table("./Enrichment_scores/Enrichment_score_NFR_0_Dnase_NDR_120_180.csv",header = TRUE, sep=",")
Enrichment_ratio4 <- read.table("./Enrichment_scores/Enrichment_score_NFR_0_Dnase_NDR_120_180_expressed.csv",header = TRUE, sep=",")
Enrichment_ratio5 <- read.table("./Enrichment_scores/Enrichment_score_NFR_0_Enhancer_diff_0.2_NDR_Enhancer_conserved_0.8_120_180.csv",header = TRUE, sep=",")

PTF <- read.table("./PTF_Cell_differentiation.txt",header = TRUE, sep = ",",stringsAsFactors = FALSE)
PTF <- PTF$PTF_name

Enrichment_ratio1$flag <-0
Enrichment_ratio2$flag <-0
Enrichment_ratio3$flag <-0
Enrichment_ratio4$flag <-0
Enrichment_ratio5$flag <-0

Enrichment_ratio1[which(Enrichment_ratio1$TF_name %in% PTF),8] <-1  
Enrichment_ratio2[which(Enrichment_ratio2$TF_name %in% PTF),9] <-1  
Enrichment_ratio3[which(Enrichment_ratio3$TF_name %in% PTF),8] <-1  
Enrichment_ratio4[which(Enrichment_ratio4$TF_name %in% PTF),9] <-1  
Enrichment_ratio5[which(Enrichment_ratio5$TF_name %in% PTF),8] <-1  

### 
Enrichment_ratio <- Enrichment_ratio1

roc <- roc.curve(scores.class0 = Enrichment_ratio$odds_ratio, weights.class0 = Enrichment_ratio$flag,curve = TRUE)
pr <- pr.curve(scores.class0 = Enrichment_ratio$odds_ratio, weights.class0 = Enrichment_ratio$flag,curve = TRUE)

png(filename = paste("./ROC/", file_name,"_ROC.png"),width = 4,
    height = 3,units = "in", res = 300)
plot(roc)
dev.off()

png(filename = paste("./ROC/", file_name,"_prROC.png"),width = 4,
    height = 3,units = "in", res = 300)
plot(pr)
dev.off()

baseline1 = nrow(Enrichment_ratio1[Enrichment_ratio1$flag==1,])/nrow(Enrichment_ratio1)
baseline2 = nrow(Enrichment_ratio2[Enrichment_ratio2$flag==1,])/nrow(Enrichment_ratio2)
baseline3 = nrow(Enrichment_ratio3[Enrichment_ratio3$flag==1,])/nrow(Enrichment_ratio3)
baseline4 = nrow(Enrichment_ratio4[Enrichment_ratio4$flag==1,])/nrow(Enrichment_ratio4)
baseline5 = nrow(Enrichment_ratio5[Enrichment_ratio5$flag==1,])/nrow(Enrichment_ratio5)

######## Calculation of MCC score
df_mcc <- data.frame(cut_off = numeric(), mcc = numeric(), stringsAsFactors = FALSE)
for (cut_off in seq(0,0.7,by = 0.01)){
  preds <- Enrichment_ratio$odds_ratio
  actuals <- Enrichment_ratio$flag
  preds[which(preds < cut_off)] = 0
  preds[which(preds >= cut_off)] = 1
  mcc <- mcc(preds, actuals)
  df_mcc = rbind(df_mcc, c(cut_off,mcc),stringsAsFactors = FALSE)
  colnames(df_mcc) <- c("cut_off", "mcc")
}

df_mcc %>%
  ggplot( aes(x=cut_off, y=mcc)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("MCC") +
  xlab("Enrichment Ratio (Threshold)")  +
  theme(panel.border = element_rect(size =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 14,face="bold"),
    axis.text.y = element_text(size = 13,face="bold"),
    axis.title.x = element_text(size = 14,face="bold"),
    axis.text.x = element_text(size = 13,face="bold"))
ggsave(paste("./ROC/mcc_",file_name,".png",sep=''), width = 4, height = 3, units = "in", dpi = 300)
write.table(df_mcc,paste("./ROC/mcc_", file_name,".csv",sep=""),quote = FALSE,row.names = FALSE,sep = ",")