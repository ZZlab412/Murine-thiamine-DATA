
###OPLS-DA
###install.packages ropls 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")
library(ropls)
x=read.csv("Extended Data Fig.3a mND_vs_mHFD_OPLSDA.csv")
row.names(x)=x[,1]
x=x[,-1]
x = t(x)
group=c(rep('mND',9),
        rep('mHFD',9))
sacurine.oplsda <- opls(x, group, predI = 1, orthoI = NA)

###PCA
library(ggplot2)
library(tidyr)
library(dplyr)
data=read.csv("Fig.3a mND_vs_mHFD_PCA.csv")
row.names(data)=data[,1]
data=data[,-1]
data <- t(data)  
group=c(rep('mND',9),
        rep('mHFD',9))

pca_data <- prcomp(data)
screeplot(pca_data, type = "lines")
summary(pca_data)
rownames(pca_data$x)
x_per <- round(summary(pca_data)$importance[2, 1]*100, 1)
y_per <- round(summary(pca_data)$importance[2, 2]*100, 1)
df_sample <- data.frame(name=rownames(pca_data$x), pca_data$x) %>%
  left_join(group, by = "name")
df_sample
theme_set(theme_bw())
ggplot(df_sample,aes(x=PC1,y=PC2,color = Group))+ geom_point(aes(color=Group),size=3)+
  stat_ellipse(level = 0.95, show.legend = F) + 
  annotate("text", label ="ND", x =0.003, y =0.001,size=5, colour ="#33CCCC") +        
   annotate("text", label ="HFD", x =0.025, y = 0,size=5, colour ="#f8766d") +
  xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
  ylab(paste("PC2","(", y_per,"%)",sep=" "))

ggplot(df_sample,aes(x=PC1,y=PC2,color = Group))+ geom_point(aes(color=Group),size=3)+
  stat_ellipse(aes(fill=Group),geom="polygon",alpha=0.2,color=NA, level = 0.95)  +
  xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
  ylab(paste("PC2","(", y_per,"%)",sep=" "))

###Volcano Plot
rm(list = ls())  

library(ggpubr)
library(ggthemes)   
data=read.csv("Extended Data Fig.3b mND_vs_mHFD_VolcanoPlot.csv")
head(data)
data$logP <- -1*log10(data$Pvalue)
ggscatter(data, x = "log2FC" , y = "logP") + theme_base()
ggscatter(data, x = "log2FC" , y = "logP", color = "regulated", 
          palette = c("sienna1","#DCDCDC","skyblue"), size = 1.5) +
  theme_base()+
  geom_hline(yintercept = -1*log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed")

###Differential metabolite heatmap
library("pheatmap")
x=read.csv("Fig.3c mND_vs_mHFD_diff_heatmap_Vitamin.csv")  ### or x=read.csv("mND_vs_mHFD_diff_heatmap_Vitamin.csv") 
row.names(x)=x[,1]
x=x[,-1]
group=c(rep('mND',9),
            rep('mHFD',9))
annotation_c <- data.frame(group)
rownames(annotation_c) <- colnames(x)
pheatmap(x, color=colorRampPalette(c("royalblue", "snow","orange"))(100), 
         scale='row', annotation_col =annotation_c, border_color=NA, width=10, 
         height=18, fontsize=10, fontsize_row=10, show_rownames = FALSE)

###KEGG
pathway=read.csv("Fig.3b KEGG.csv") 
ggplot(pathway,aes(X,Description))+
  geom_point(aes(size=-log10(qvalue),color=pvalue))+
  scale_color_gradient2(low = "#ff7300",
                        mid = "#d94dff", high = "#4169e1",midpoint =0.45)+
  theme_bw()  

####Metabolites and phenotypes----spearman Correlation Analysis
#Overall horizontal correlation
rm(list = ls())
data<-read.csv('Extended Data Fig.3c Metabolites_phenotypes_Correlation.csv', row.names = 1)
data <- t(data)
library(corrplot)
library(psych)
spearman1 <- corr.test(data, method = 'spearman')
spearman1[[11]]
spearman1$p
spearman1$r
col4 <- colorRampPalette(c( "#007FFF","#FFE4E1","#FF4500"))
corrplot(spearman1$r, method = 'square', type = 'lower', 
         p.mat = spearman1$p, insig = 'label_sig', sig.level = c(0.001,0.01,0.05),  #P值小于0.001,0.01,0.05显示*，*越多越显著
         pch.cex = 1, pch.col = "black", #改变星星大小及颜色
         diag = FALSE,  col = col4(100), tl.col = '#000000')

###One to one correlation
library(ggplot2)
library(corrplot)
library(ggpubr)
library(psych)
data <- read.delim('Fig.3c mND_vs_mHFD_diff_heatmap_Vitamin.csv', row.names = 1, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
data <- t(data)
spearman1 <- corr.test(data, method = 'spearman',adjust = "none")
spearman1[[11]] 
p <- spearman1$p
r <- spearman1$r
write.table(p,"p值.csv",row.names=T,col.names=TRUE,sep=",")
write.table(r,"r值.csv",row.names=T,col.names=TRUE,sep=",")
colnames(data)
data1 = data.frame(data)
p1 <- ggplot(data = data1, mapping = aes(x =Thiamine,y = 0dpp_cyst)) +     #Replace X and Y
  geom_point(colour = "#F08080", size = 4) +  geom_smooth(method = lm,colour= "#4169E1",fill= "#4169E1")

p1
p2 <- p1+ stat_cor(method = "spearman",label.x.npc = 'left',label.y.npc = 'top')
p2