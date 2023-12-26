library(cluster)
library(factoextra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(Rtsne)
library(stringr)

setwd("~/Desktop/VOR_eLife/PTF_DATA/TF_binding_profiles_Analysis/Codes")

df <- read.table("./ChIP_seq_TF_motif_file.Info",header = TRUE)

H1_TF <- df[df$Cellline=="H1",1]
HepG2_TF <- df[df$Cellline=="HepG2",1]
HeLa_TF <- df[df$Cellline=="HeLa-S3",1]
K562_TF <- df[df$Cellline=="K562",1]
MCF7_TF <- df[df$Cellline=="MCF-7",1]

file_name  <- "NFR_0_146_148_half"
flank_length<-0

Read_TF_Celline <- function(cell_line, TFs) {
  TF_motif_counts <- data.frame()
   for (TF_name in unique(TFs)) {
    tryCatch({
      motif_region_counts <-read.csv(paste("../TF_binding_motif_profiles_nucleosome_dyad_flank_1000/",cell_line,"_", TF_name,"_NFR_0_flank_1000_146_148_motif_region_counts.txt",sep=''),header = TRUE)
# Select TFs with the total genome-wide cumulative sum of ChiP-seq TF motifs on nucleosome regions (-60,60) > 500 base pairs
      if (sum(motif_region_counts$counts_region[944:(nrow(motif_region_counts)-943)]) >=500 ) {
      if (nrow(TF_motif_counts)==0) {
        TF_motif_counts <- select(motif_region_counts, "binding_sites_region", "counts_region")
        colnames(TF_motif_counts) <- c("binding_sites", paste(TF_name,str_replace(cell_line, "-", ""),sep="_"))
      }
      else{
        TF_motif_counts <- cbind(TF_motif_counts, motif_region_counts$counts_region)
        colnames(TF_motif_counts) <- c(colnames(TF_motif_counts)[1:(ncol(TF_motif_counts)-1)], paste(TF_name,str_replace(cell_line, "-", ""),sep="_"))
      }
      }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(list(TF_motif_counts))
}

all_TF_motif_counts <- data.frame()
## Read the raw TF binding profiles
all_TF_motif_counts <- Read_TF_Celline("K562",K562_TF)[[1]]
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("HepG2",HepG2_TF)[[1]][,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("HeLa-S3",HeLa_TF)[[1]][,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("MCF-7",MCF7_TF)[[1]][,-1])
all_TF_motif_counts <- cbind(all_TF_motif_counts, Read_TF_Celline("H1",H1_TF)[[1]][,-1])

all_TF_motif_counts_NFR <- all_TF_motif_counts %>% filter(binding_sites >=-60 & binding_sites <= 60) # take binding motif profiles of regions of +/- 60 base pair from dyad 

## Calculate Pearson correlation coefficient of motif profiles between two symmetrical nucleosomal halves
all_TF_motif_CC <- data.frame()
R_len = 61 ## length of nucleosome regions (half)

for (i in 2:ncol(all_TF_motif_counts_NFR)){
  if (nrow(all_TF_motif_CC)==0) {
    all_TF_motif_CC <- data.frame(TF_cell_line = colnames(all_TF_motif_counts_NFR)[i], 
                                  CC = cor(all_TF_motif_counts_NFR[1:R_len,i], 
                                           rev(all_TF_motif_counts_NFR[R_len:(R_len*2-1),i]), 
                                           method = c("pearson")),stringsAsFactors = FALSE) 
  }
  else{
    all_TF_motif_CC <- rbind(all_TF_motif_CC, 
                             c(colnames(all_TF_motif_counts_NFR)[i], 
                               cor(all_TF_motif_counts_NFR[1:R_len,i], 
                                   rev(all_TF_motif_counts_NFR[R_len:(R_len*2-1),i]),
                                   method = c("pearson"))), stringsAsFactors = FALSE)
    colnames(all_TF_motif_CC) <- c("TF_cell_line","CC")
    all_TF_motif_CC$TF_cell_line <-as.character(all_TF_motif_CC$TF_cell_line)
    all_TF_motif_CC$CC <-as.numeric(all_TF_motif_CC$CC)
  }
}

TF_name_CC <- all_TF_motif_CC %>% 
  mutate(CC = as.double(CC)) %>%
  filter(CC>=0.4) %>%  #removed TFs with Pearson correlation coefficient values less than 0.4.
  select("TF_cell_line")

## make plots for CC distributions
ggplot(all_TF_motif_CC,  aes(x=CC)) +
  geom_density(fill="#69b3a2", color="black", alpha=0.8) +
  theme_bw() +
  ylab("Density") +
  xlab("Pearson correlation coefficient") + 
  theme(panel.border = element_rect(linewidth =1.5, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.y = element_text(size = 14,face="bold"),
    axis.text.y = element_text(size = 13,face="bold"),
    axis.title.x = element_text(size = 14,face="bold"),
    axis.text.x = element_text(size = 13,face="bold"))
ggsave(paste("./cluster/CC_SHL_sym_density.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

write.csv(all_TF_motif_CC, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/TF_binding_profile_Correlation_two_halves.csv", row.names = FALSE)
write.csv(all_TF_motif_CC, file = "./cluster/TF_binding_profile_Correlation_two_halves.csv", row.names = FALSE)


## take binding motif profiles of regions of +/- 60 base pair from dyad 
all_TF_motif_counts_filter <-  all_TF_motif_counts %>%  filter(binding_sites >=-60 & binding_sites <= 60) %>% select("binding_sites", TF_name_CC$TF_cell_line) 

## Normalize the binding profile to range of (0,1)
maxs <- apply(all_TF_motif_counts_filter, 2, max)
mins <- apply(all_TF_motif_counts_filter, 2, min)

all_TF_motif_counts_scaled <- as.data.frame(scale(all_TF_motif_counts_filter, 
                                                  center = mins, scale = maxs - mins))


## Binding motif profiles between two symmetrical nucleosomal halves are combined for each TF.
all_TF_motif_counts_scaled_invert <- data.frame(t(all_TF_motif_counts_scaled))
colnames(all_TF_motif_counts_scaled_invert) <- all_TF_motif_counts_filter$binding_sites
half_length <- (ncol(all_TF_motif_counts_scaled_invert) +1)/2
all_TF_motif_counts_scaled_invert_half <- all_TF_motif_counts_scaled_invert[,half_length]

for (i in seq(1:(half_length-1)))
{
  all_TF_motif_counts_scaled_invert_half <- cbind(all_TF_motif_counts_scaled_invert_half,
                                                  (all_TF_motif_counts_scaled_invert[,half_length-i] + 
                                                     all_TF_motif_counts_scaled_invert[,half_length+i])/2 ) 
}

colnames(all_TF_motif_counts_scaled_invert_half) <- 0:(half_length-1)
rownames(all_TF_motif_counts_scaled_invert_half) <- rownames(all_TF_motif_counts_scaled_invert)
all_TF_motif_counts_scaled_half <- data.frame(t(all_TF_motif_counts_scaled_invert_half))
all_TF_motif_counts_scaled_half$binding_sites <- 0:(half_length-1)


set.seed(0)
# Apply t-distributed stochastic neighbor embedding (t-SNE) and projected all profiles onto two dimensions
tSNE_fit <- all_TF_motif_counts_scaled_invert_half[-1,] %>% Rtsne(perplexity=12,dims = 2,theta=0,initial_dims=61)

tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
           tSNE2="V2") %>%
mutate(ID=rownames(all_TF_motif_counts_scaled_invert[-1,]))
  rownames(tSNE_df) <- tSNE_df$ID
  
## Perform K-medoids clusering
pamx <- pam(tSNE_df[,1:2], 6)

## make plots to check cluster quality
png(file=paste("./cluster/",file_name,"_silhouette_distribution.png",sep=""),width=4,height=3,units="in",res=300)
fviz_silhouette(pamx, palette = "jco", ggtheme = theme_classic())
dev.off()

png(file=paste("./cluster/",file_name,"_silhouette.png",sep=""),width=4,height=3,units="in",res=300)
fviz_nbclust(tSNE_df[,1:2], pam, method ="silhouette",
             barfill = "white",
             barcolor = "white",
             linecolor = "black")
dev.off()

png(file=paste("./cluster/",file_name,"_cluster_distribution.png",sep=""),width=4,height=3,units="in",res=300)
fviz_cluster(pamx, geom="point",
             #palette =c("#007892","#D9455F","#00FFFF"),
             palette =c("#007892","#D9455F","#00FFFF", "#00FF00", "#FFFF00", "#008080","#800000"),
             ellipse.type ="euclid",
             repel =TRUE,
             ggtheme =theme_minimal())
dev.off()

## save each cluster and we renumbered each cluster to make it consistent with the plots
cluster1 <- names(pamx$clustering[which(pamx$clustering==2)])
cluster2 <- names(pamx$clustering[which(pamx$clustering==3)])
cluster3 <- names(pamx$clustering[which(pamx$clustering==5)])
cluster4 <- names(pamx$clustering[which(pamx$clustering==6)])
cluster5 <- names(pamx$clustering[which(pamx$clustering==4)])
cluster6 <- names(pamx$clustering[which(pamx$clustering==1)])

# members with silhouette width <= 0.25 were considered as outliers and 33 outliers removed.
silhouette_width <-data.frame(TF_name = rownames(silhouette(pamx)), sil_width=silhouette(pamx)[,3])
cluster_filter <- silhouette_width %>% filter(sil_width>=0.25) %>% select(TF_name) 
cluster1 <- intersect(cluster1,cluster_filter$TF_name) 
cluster2 <- intersect(cluster2,cluster_filter$TF_name)  
cluster3 <- intersect(cluster3,cluster_filter$TF_name)  
cluster4 <- intersect(cluster4,cluster_filter$TF_name)  
cluster5 <- intersect(cluster5,cluster_filter$TF_name)  
cluster6 <- intersect(cluster6,cluster_filter$TF_name)  

write.table(sort(cluster1),paste("./cluster/Cluster1_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sort(cluster2),paste("./cluster/Cluster2_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sort(cluster3),paste("./cluster/Cluster3_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sort(cluster4),paste("./cluster/Cluster4_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sort(cluster5),paste("./cluster/Cluster5_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sort(cluster6),paste("./cluster/Cluster6_",file_name,".txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


data_mids_combine_cluster1 <- all_TF_motif_counts_scaled_half[cluster1]
Mean_binding_cluster1 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster1))
data_mids_combine_cluster1$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster2 <- all_TF_motif_counts_scaled_half[cluster2]
Mean_binding_cluster2 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster2))
data_mids_combine_cluster2$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster3 <- all_TF_motif_counts_scaled_half[cluster3]
Mean_binding_cluster3 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster3))
data_mids_combine_cluster3$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster4 <- all_TF_motif_counts_scaled_half[cluster4]
Mean_binding_cluster4 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster4))
data_mids_combine_cluster4$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster5 <- all_TF_motif_counts_scaled_half[cluster5]
Mean_binding_cluster5 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster5))
data_mids_combine_cluster5$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster6 <- all_TF_motif_counts_scaled_half[cluster6]
Mean_binding_cluster6 <- data.frame(binding_sites = all_TF_motif_counts_scaled_half$binding_sites, 
                                    Mean=rowMeans(data_mids_combine_cluster6))
data_mids_combine_cluster6$binding_sites <- all_TF_motif_counts_scaled_half$binding_sites

data_mids_combine_cluster1_pivot <- data_mids_combine_cluster1 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster2_pivot <- data_mids_combine_cluster2 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster3_pivot <- data_mids_combine_cluster3 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster4_pivot <- data_mids_combine_cluster4 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster5_pivot <- data_mids_combine_cluster5 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")

data_mids_combine_cluster6_pivot <- data_mids_combine_cluster6 %>%
  pivot_longer(!binding_sites, names_to = "TFs", values_to = "Counts")


## make plots for each cluster
# Cluster1
ggplot(data_mids_combine_cluster1_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25) + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                                          limits = c(0,65)) + 
  geom_line(data=Mean_binding_cluster1, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  ggtitle("Cluster1") + 
  scale_color_manual(values = c(rep("grey",length(cluster1)))) +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_blank(),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=5, 
                                   face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./cluster/",file_name,"_Cluster1.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)


#Cluster2
ggplot(data_mids_combine_cluster2_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25) + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                     limits = c(0,65)) +    
  geom_line(data=Mean_binding_cluster2, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  ggtitle("Cluster2") + 
  scale_color_manual(values = c(rep("grey",length(cluster2)))) +
  #scale_color_discrete(name="")  +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_blank(),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=5, 
                                   face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./cluster/",file_name,"_Cluster2.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

#Cluster3
ggplot(data_mids_combine_cluster3_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25) + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                     limits = c(0,65)) +   
  geom_line(data=Mean_binding_cluster3, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  ggtitle("Cluster3") + 
  scale_color_manual(values = c(rep("grey",length(cluster3)))) +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_blank(),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=5, 
                                   face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./cluster/",file_name,"_Cluster3.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

#Cluster4
ggplot(data_mids_combine_cluster4_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25) + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                     limits = c(0,65)) +  
  geom_line(data=Mean_binding_cluster4, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  ggtitle("Cluster4") + 
  scale_color_manual(values = c(rep("grey",length(cluster4)))) +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_blank(),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=5, 
                                   face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./cluster/",file_name,"_Cluster4.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

#Cluster5
ggplot(data_mids_combine_cluster5_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
  geom_line(size = 0.25) + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                     limits = c(0,65)) +    
  geom_line(data=Mean_binding_cluster5, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
  xlab("Distance from Dyad (bp)") +
  ylab("Number of motif (normalized)") +
  theme_bw() +
  ggtitle("Cluster5") + 
  scale_color_manual(values = c(rep("grey",length(cluster5)))) +
  theme(axis.text = element_text(size =15, face="bold"),
        axis.title = element_blank(),
        plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
        legend.text = element_text(colour="black", size=5, 
                                   face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
ggsave(paste("./cluster/",file_name,"_Cluster5.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)

#Cluster6
  ggplot(data_mids_combine_cluster6_pivot, aes(x = binding_sites, y = Counts, color = TFs)) + 
    geom_line(size = 0.25) + 
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
                       limits = c(0,65)) +  
    geom_line(data=Mean_binding_cluster6, aes(x = binding_sites, y = Mean), color="black", size =1.5, linetype = "solid") + 
    xlab("Distance from Dyad (bp)") +
    ylab("Number of motif (normalized)") +
    theme_bw() +
    ggtitle("Cluster6") + 
    scale_color_manual(values = c(rep("grey",length(cluster6)))) +
    theme(axis.text = element_text(size =15, face="bold"),
          axis.title = element_blank(),
          plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5),
          legend.text = element_text(colour="black", size=5, 
                                     face="bold"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
    theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm")) 
  ggsave(paste("./cluster/",file_name,"_Cluster6.png",sep=''), width = 4, height = 3, units = "in", dpi = 300)
  
write.csv(data_mids_combine_cluster1, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster1.csv", row.names = FALSE)
write.csv(data_mids_combine_cluster2, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster2.csv", row.names = FALSE)
write.csv(data_mids_combine_cluster3, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster3.csv", row.names = FALSE)
write.csv(data_mids_combine_cluster4, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster4.csv", row.names = FALSE)
write.csv(data_mids_combine_cluster5, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster5.csv", row.names = FALSE)
write.csv(data_mids_combine_cluster6, file = "/Users/yunhuipeng/Desktop/VOR_eLife/Source_data/data_mids_combine_cluster6.csv", row.names = FALSE)
