library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggpubr)

setwd("~/Desktop/VOR_eLife/PTF_DATA/TF_binding_profiles_Analysis/Codes")

## get list of TFs for each cell lines
df <- read.table("./ChIP_seq_TF_motif_file.Info",header = TRUE)
H1_TF <- df[df$Cellline=="H1",1]
HepG2_TF <- df[df$Cellline=="HepG2",1]
HeLa_TF <- df[df$Cellline=="HeLa-S3",1]
K562_TF <- df[df$Cellline=="K562",1]
MCF7_TF <- df[df$Cellline=="MCF-7",1]

CC_type<-  "pearson" # pearson or spearman
file_name  <- "End_preference_CC_0.4_ratio"

Read_TF_Celline <- function(cell_line, TFs) {
  all_TF_motif_counts <- data.frame()
  for (name in unique(TFs)) {
    tryCatch({
      motif_region_counts <-read.csv(paste("../TF_binding_motif_profiles_nucleosome_dyad_flank_1000/",cell_line,"_", name,"_NFR_0_flank_1000_146_148_motif_region_counts.txt",sep=''),header = TRUE)
      if (sum(motif_region_counts$counts_region[944:(nrow(motif_region_counts)-943)]) >=5000 ) {
        if (nrow(all_TF_motif_counts)==0) {
          all_TF_motif_counts <- select(motif_region_counts, "binding_sites_region", "counts_region")
          colnames(all_TF_motif_counts) <- c("binding_sites", paste(name,str_replace(cell_line, "-", ""),sep="_"))
        }
        else{
          all_TF_motif_counts <- cbind(all_TF_motif_counts, motif_region_counts$counts_region)
          colnames(all_TF_motif_counts) <- c(colnames(all_TF_motif_counts)[1:(ncol(all_TF_motif_counts)-1)], paste(name,str_replace(cell_line, "-", ""),sep="_"))
        }
      }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(all_TF_motif_counts)
}


all_TF_motif_counts <- Read_TF_Celline("K562",K562_TF)
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("HepG2",HepG2_TF)[,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("HeLa-S3",HeLa_TF)[,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("MCF-7",MCF7_TF)[,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("H1",H1_TF)[,-1])

## ## take binding motif profiles of regions of +/- 60 base pair from dyad 
all_TF_motif_counts_NFR <- all_TF_motif_counts %>% filter(binding_sites >=-60 & binding_sites <= 60)

## Calculate Pearson correlation coefficient of motif profiles between two symmetrical nucleosomal halves
all_TF_motif_CC <- data.frame(TF_cell_line=character(), CC=numeric(), stringsAsFactors = FALSE)

R_len =  61 #R_len = 61 ## length of nucleosome regions (half)

for (i in 2:ncol(all_TF_motif_counts_NFR)){
    all_TF_motif_CC <- rbind(all_TF_motif_CC, 
                             c(colnames(all_TF_motif_counts_NFR)[i], 
                               cor(all_TF_motif_counts_NFR[1:R_len,i], 
                                   rev(all_TF_motif_counts_NFR[R_len:(R_len*2-1),i]),
                                   method = c(CC_type))),stringsAsFactors = FALSE)
    colnames(all_TF_motif_CC) <- c("TF_cell_line","CC")
  }

TF_name_CC <- all_TF_motif_CC %>% 
  mutate(CC = as.double(CC)) %>%
  filter(CC>=0.4) %>%  ##removed TFs with Pearson correlation coefficient values less than 0.4.
  select("TF_cell_line")

### ## take binding motif profiles of regions of +/- 90 base pair from dyad 
all_TF_motif_counts_filter <-  all_TF_motif_counts %>% filter(binding_sites >=-90 & binding_sites <= 90)
data_mids_combine <- all_TF_motif_counts_filter[TF_name_CC$TF_cell_line] ## ##removed TF profiles with Pearson correlation coefficient values less than 0.4.
maxs <- apply(data_mids_combine, 2, max)
mins <- apply(data_mids_combine, 2, min)

data_mids_combine <- as.data.frame(scale(data_mids_combine, center = mins, scale = maxs - mins))

data_mids_combine$binding_sites <- all_TF_motif_counts_filter$binding_sites


data_mids_combine_pivot <- data_mids_combine %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")


Dyad_counts <- colSums(data_mids_combine[which(abs(data_mids_combine$binding_sites) <=15 ),])
End_counts <- colSums(data_mids_combine[which(abs(data_mids_combine$binding_sites) >=55 & abs(data_mids_combine$binding_sites) <=73),])
All_counts <- colSums(data_mids_combine[which(abs(data_mids_combine$binding_sites) >=0 & abs(data_mids_combine$binding_sites) <=70),])

TF_cellline_names <- names(Dyad_counts)
TF_names <-vector()
for (i in TF_cellline_names) {
  TF_names <- append(TF_names, str_split(i,"_")[[1]][1])
}

End_Dyad_Ratio <- data.frame(TF_cellline <-names(Dyad_counts), TF <- TF_names, 
                             End_Dyad_ratio <-End_counts/Dyad_counts)
colnames(End_Dyad_Ratio) <- c("TF_cellline", "TF", "End_Dyad_Ratio")


### classifications of TFs
TF_classfy <- read.table("../TF_classification.csv",header = TRUE, sep = ",") 
TF_classfy <- TF_classfy %>% select(TF, EMI.penetration..lig147.,Dyad.EMI.intensity,End.preference, Dyad.preference)
TF_classfy_Dyad_End_Ratio <- merge(End_Dyad_Ratio,TF_classfy, by= "TF") #%>% na.omit()
colnames(TF_classfy_Dyad_End_Ratio) <- c("TF", "TF_cellline","End_Dyad_Ratio","EMI.penetration", "Dyad.EMI.intensity","End.preference","Dyad.preference")

## Remove one outlier
TF_classfy_Dyad_End_Ratio <- TF_classfy_Dyad_End_Ratio %>% filter(End.preference=="Yes"&EMI.penetration<15)
   

# Make plots
ggplot(TF_classfy_Dyad_End_Ratio, aes(x = End_Dyad_Ratio, y = EMI.penetration)) + 
  geom_point() + 
  stat_cor(method="pearson", label.x.npc = 'middle') +
  xlab("End/Dyad binding ratio") +
  ylab("EMI.penetration") +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_text(size =16, face="bold"),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=15, 
                                   face="bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./",file_name,".png",sep=''), width = 6, height = 4, units = "in", dpi = 300)

write.csv(TF_classfy_Dyad_End_Ratio, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_files/cluster/TF_classfy_Dyad_End_Ratio.csv", row.names = FALSE)
write.csv(TF_classfy_Dyad_End_Ratio, file = "./cluster/TF_classfy_Dyad_End_Ratio.csv", row.names = FALSE)