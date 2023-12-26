library(tidyr)
library(dplyr)
library(ggplot2)
library(TTR)
library(reshape2)
library(RColorBrewer)
library(stringr)

df <- read.table("./ChIP_seq_TF_motif_file.Info",header = TRUE)
H1_TF <- df[df$Cellline=="H1",1]
HepG2_TF <- df[df$Cellline=="HepG2",1]
HeLa_TF <- df[df$Cellline=="HeLa-S3",1]
K562_TF <- df[df$Cellline=="K562",1]
MCF7_TF <- df[df$Cellline=="MCF-7",1]

out_file_name  <- "NFR_0_flank_1000_146_148_occupancy_CC"

Read_TF_Celline <- function(cell_line, TFs) {
  all_TF_motif_counts <- data.frame()
  for (name in unique(TFs)) {
    tryCatch({
      motif_region_counts <-read.csv(paste("../TF_binding_motif_profiles_nucleosome_dyad_flank_1000/",cell_line,"_", name,"_NFR_0_flank_1000_146_148_motif_region_counts.txt",sep=''),header = TRUE)
      motif_counts_smooths <- WMA(motif_region_counts$counts_region, n =31, wts = c(1:16,15:1)) # Smooth the profile with a 30 bp window
      motif_region_counts$counts_smooths_region[16:(length(motif_counts_smooths)-15)] <- motif_counts_smooths[31:length(motif_counts_smooths)]
      
      if (sum(motif_region_counts$counts_region) >=10000 ) {
        if (nrow(all_TF_motif_counts)==0) {
          all_TF_motif_counts <- select(motif_region_counts, "binding_sites_region", "counts_smooths_region")
          colnames(all_TF_motif_counts) <- c("binding_sites", paste(name,str_replace(cell_line, "-", ""),sep="_"))
        }
        else{
          all_TF_motif_counts <- cbind(all_TF_motif_counts, motif_region_counts$counts_smooths_region)
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

all_TF_motif_counts <- all_TF_motif_counts %>% filter(binding_sites >= -400 & binding_sites <= 400) ## Take the genomic regions centered on TF motif and flanked by 400 bp

## Calcuate the Correlation coerraltion between TF binding profile and ith nucleosome occupancy levels around nucleosome dyad
all_TF_motif_nucleosome_occupancy_CC <- data.frame(TF_cell_line=character(), CC= numeric(), pvalue=numeric(), stringsAsFactors = FALSE)
for (cell_line in c("MCF-7", "H1", "HepG2", "K562", "HeLa-S3")) {
  for (name in unique(df$TF_name)) {
    tryCatch({
      motif_region_counts <-read.csv(paste("../TF_binding_motif_profiles_nucleosome_dyad_flank_1000/",cell_line,"_", name,"_NFR_0_flank_1000_146_148_motif_region_counts.txt",sep=''),header = TRUE)
      motif_counts_smooths <- WMA(motif_region_counts$counts_region, n =51, wts = c(1:26,25:1))
      motif_region_counts$counts_smooths_region[26:(length(motif_counts_smooths)-25)] <- motif_counts_smooths[51:length(motif_counts_smooths)]
      
      motif_region_counts_Occupancy <-read.csv(paste("../Nucleosome_Occupancy_around_dyad_all_NFR/",cell_line,"_NFR_0_flank_1000_nucleosome_occupancy_all_NFR_motif_region_counts.txt",sep=''),header = TRUE)
      motif_counts_smooths_Occupancy <- WMA(motif_region_counts_Occupancy$counts_region, n =51, wts = c(1:26,25:1))
      motif_region_counts_Occupancy$counts_smooths_region[26:(length(motif_counts_smooths_Occupancy)-25)] <- motif_counts_smooths_Occupancy[51:length(motif_counts_smooths_Occupancy)]
      
      if (sum(motif_region_counts$counts_region) >=10000 ) {
      cc_oc <-cor.test(motif_region_counts$counts_smooths_region[604:(nrow(motif_region_counts)-603)], ## Take the genomic regions centered on TF motif and flanked by 400 bp
          motif_region_counts_Occupancy$counts_smooths_region[604:(nrow(motif_region_counts)-603)],
          method = c("pearson"))
      
      all_TF_motif_nucleosome_occupancy_CC <- rbind(all_TF_motif_nucleosome_occupancy_CC, 
                                                    c(paste(name,str_replace(cell_line, "-", ""),sep="_"), 
                                                      cc_oc$estimate, cc_oc$p.value), stringsAsFactors = FALSE)
      colnames(all_TF_motif_nucleosome_occupancy_CC) <- c("TF_cell_line","CC","pvalue")
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

all_TF_motif_nucleosome_occupancy_CC$CC <- as.numeric(all_TF_motif_nucleosome_occupancy_CC$CC)
all_TF_motif_nucleosome_occupancy_CC$pvalue <- as.numeric(all_TF_motif_nucleosome_occupancy_CC$pvalue)
all_TF_motif_nucleosome_occupancy_CC <- all_TF_motif_nucleosome_occupancy_CC %>% arrange(desc(CC))


#TFs with positive correlation coefficients between motif profile and nucleosome occupancy
TF_name_occupancy_CC <- all_TF_motif_nucleosome_occupancy_CC %>% 
  mutate(CC = as.double(CC),pvalue = as.numeric(pvalue)) %>%
  filter(CC>= 0.2 & pvalue <0.05) %>% 
  select("TF_cell_line")

##TFs with negative correlation coefficients between motif profile and nucleosome occupancy
TF_name_occupancy_CC_rev <- all_TF_motif_nucleosome_occupancy_CC %>% 
  mutate(CC = as.double(CC),pvalue = as.numeric(pvalue)) %>%
  filter(CC<=-0.4  & pvalue <0.05) %>% 
  select("TF_cell_line")

##TFs with weak correlation coefficients between motif profile and nucleosome occupancy
TF_name_occupancy_CC_weak <- all_TF_motif_nucleosome_occupancy_CC %>% 
  mutate(CC = as.double(CC),pvalue = as.numeric(pvalue)) %>%
  filter(CC >= -0.1 & CC <= 0.1 & pvalue < 0.05) %>% 
  select("TF_cell_line")

ggplot(all_TF_motif_nucleosome_occupancy_CC,  aes(x=CC)) +
  geom_histogram(fill="#69b3a2", color="black", alpha=0.8) +
  theme_bw() +
  #geom_vline(xintercept = -0.464,size=1, lty="dashed", color = "red") +
  ylab("Frequency") +
  xlab("Pearson correlation coefficient") + 
  theme(panel.border = element_rect(size =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 14,face="bold"),
    axis.text.y = element_text(size = 13,face="bold"),
    axis.title.x = element_text(size = 14,face="bold"),
    axis.text.x = element_text(size = 13,face="bold"))
ggsave(paste("./All_TF_nucleosome_occupancy_CC_density_all_NFR.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

all_TF_motif_counts_filter <- all_TF_motif_counts

## Normalize the TF binding profile and nucleosome occupancy levels to (0,1)
maxs <- apply(all_TF_motif_counts_filter, 2, max)
mins <- apply(all_TF_motif_counts_filter, 2, min)

all_TF_motif_counts_scaled <- as.data.frame(scale(all_TF_motif_counts_filter, center = mins, scale = maxs - mins))
all_TF_motif_counts_scaled$binding_sites <- all_TF_motif_counts_filter$binding_sites

#Binding motif profiles of TFs with positive correlation coefficients between motif profile and nucleosome occupancy
TF_cluster1 <- TF_name_occupancy_CC$TF_cell_line
##Binding motif profiles of TFs with negative correlation coefficients between motif profile and nucleosome occupancy
TF_cluster2  <- TF_name_occupancy_CC_rev$TF_cell_line
##Binding motif profiles of TFs with weak correlation coefficients between motif profile and nucleosome occupancy
TF_cluster3  <- TF_name_occupancy_CC_weak$TF_cell_line

data_mids_combine_cluster1 <- all_TF_motif_counts_scaled[TF_cluster1]
Mean_binding_cluster1 <- data.frame(binding_sites = all_TF_motif_counts_scaled$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster1))
data_mids_combine_cluster1$binding_sites <- all_TF_motif_counts_scaled$binding_sites

data_mids_combine_cluster2 <- all_TF_motif_counts_scaled[TF_cluster2]
Mean_binding_cluster2 <- data.frame(binding_sites = all_TF_motif_counts_scaled$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster2))
data_mids_combine_cluster2$binding_sites <- all_TF_motif_counts_scaled$binding_sites

data_mids_combine_cluster3 <- all_TF_motif_counts_scaled[TF_cluster3]
Mean_binding_cluster3 <- data.frame(binding_sites = all_TF_motif_counts_scaled$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster3))
data_mids_combine_cluster3$binding_sites <- all_TF_motif_counts_scaled$binding_sites


data_mids_combine_cluster1_pivot <- data_mids_combine_cluster1 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster2_pivot <- data_mids_combine_cluster2 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster3_pivot <- data_mids_combine_cluster3 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

#### Make Plots for each cluster
# Cluster1
ggplot(data_mids_combine_cluster1_pivot, aes(x = binding_sites, y = Counts)) + 
  geom_line(size = 0.25, color="grey") + 
  scale_x_continuous(breaks = c(-400,-300,-200,-100,0,100,200,300,400),
                     limits = c(-400,400)) +     #  ylim(0.00,0.006) + 
  geom_line(data=Mean_binding_cluster1, aes(x = binding_sites, y = Mean), color="#DC143C", size =1.5, linetype = "dotted") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  theme(axis.text = element_text(size =10, face="bold"),
        axis.title = element_text(size =10, face="bold"),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./",out_file_name,"_Cluster1_all_NFR.png",sep=''), width = 6, height = 2.5, units = "in", dpi = 300)


#Cluster2
ggplot(data_mids_combine_cluster2_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25, color="grey") + 
  scale_x_continuous(breaks = c(-400,-300,-200,-100,0,100,200,300,400),
                     limits = c(-400,400)) +     #  ylim(0.00,0.006) + 
  geom_line(data=Mean_binding_cluster2, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "dotted") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  theme(axis.text = element_text(size =10, face="bold"),
        axis.title = element_text(size =10, face="bold"),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 

ggsave(paste("./",out_file_name,"_Cluster2_all_NFR.png",sep=''), width = 6, height = 2.5, units = "in", dpi = 300)

#Cluster3
ggplot(data_mids_combine_cluster3_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25, color="grey") + 
  scale_x_continuous(breaks = c(-400,-300,-200,-100,0,100,200,300,400),
                     limits = c(-400,400)) +     #  ylim(0.00,0.006) + 
  geom_line(data=Mean_binding_cluster3, aes(x = binding_sites, y = Mean), color="blue", size =1.5, linetype = "dotted") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  theme(axis.text = element_text(size =12, face="bold"),
        axis.title = element_text(size =12, face="bold"),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),

        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./",out_file_name,"_Cluster3_all_NFR.png",sep=''), width = 6, height = 2.5, units = "in", dpi = 300)