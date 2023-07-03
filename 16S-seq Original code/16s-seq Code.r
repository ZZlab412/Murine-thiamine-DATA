

###PCOA
library(vegan)  
library(ggplot2)
library(plyr)
otu1=read.csv("Supplementary Fig.5e OTU_PCOA.csv") 
group=read.csv("group.csv")
row.names(otu1) <- otu1$x   
otu1=otu1[,-1]
otu <- data.frame(t(otu1))  
distance <- vegdist(otu, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = T)
plot_data <- data.frame({pcoa$point})[1:2]
plot_data$Sample_name <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig = pcoa$eig
plot_data$group <- c(rep('B_mND', 10), rep('A_mHFD', 10)) 
ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, color=group)) +
  geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  stat_ellipse(aes(fill=group),type="norm",geom="polygon",alpha=0.2,color=NA) 



###heatmap
rm(list = ls())
library(pheatmap)
x=read.csv("Supplementary Fig.5h Heatmap.csv")
row.names(x)=x[,1]
x=x[,-1]
x1=log2(x + 1)
group=c(rep('mND',10),
        rep('mHFD',10))
annotation_c <- data.frame(group)
rownames(annotation_c) <- colnames(x)

pheatmap(x1, color=colorRampPalette(c("#008B8B", "snow","#BF3EFF"))(100), 
           scale='row', border_color='gray', width=10, 
           height=18, fontsize=10, fontsize_row=10, show_rownames = T, cluster_rows=F, 
           cluster_cols=F, gaps_col = 10)


