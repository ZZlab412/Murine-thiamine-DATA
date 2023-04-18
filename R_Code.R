
########  Single-Library Analysis with Cell Ranger in server
########  CellRanger version=3.1.0

nohup /home/.../10xgenomics/cellranger-3.1.0/cellranger count --id=SampleName \
--localcores 20 \
--transcriptome=/home/.../10xgenomics/refdata-cellranger-mm10-3.0.0 \
--fastqs=/home/.../Singlecell_Rawdata   \
--sample=SampleName \
--chemistry=threeprime \
--force-cells=6000 &

#############################################################################################
         #########///////=========Seurat version=4.0.4========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)

###PD0_mND  ###Repeat the steps to read in P0_mHFD/P3_mND/P3_mHFD
PD0_Ctrl<-Read10X(data.dir =".../filtered_feature_bc_matrix")
colnames(x = PD0_Ctrl) <- paste('PD0_Ctrl', colnames(x = PD0_Ctrl), sep = '_')
PD0_Ctrl <- CreateSeuratObject(counts = PD0_Ctrl, project = "PD0_Ctrl", min.cells = 3, min.features = 200)
PD0_Ctrl[["percent.mt"]] <- PercentageFeatureSet(PD0_Ctrl, pattern = "^mt-")
VlnPlot(PD0_Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)  

###Value setting of filteration
PD0_Ctrl <- subset(PD0_Ctrl, subset =  nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 20) 
PD0_Ctrl <- NormalizeData(PD0_Ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
PD0_Ctrl <- FindVariableFeatures(PD0_Ctrl, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
all.genes <- rownames(PD0_Ctrl)
PD0_Ctrl <- ScaleData(PD0_Ctrl, features = all.genes)
PD0_Ctrl@meta.data$tech <- "PD0_Ctrl"
PD0_Ctrl@meta.data$con <- "PD0_Ctrl"
PD0_Ctrl <- RunPCA(PD0_Ctrl)
VizDimLoadings(PD0_Ctrl, dims = 1:5, reduction = "pca")
DimPlot(PD0_Ctrl, reduction = "pca")
DimHeatmap(PD0_Ctrl, dims = 1:20, cells = 500, balanced = TRUE)
PD0_Ctrl <- JackStraw(PD0_Ctrl, num.replicate = 100)
PD0_Ctrl <- ScoreJackStraw(PD0_Ctrl, dims = 1:20)
JackStrawPlot(PD0_Ctrl, dims = 1:20)
ElbowPlot(PD0_Ctrl)
PD0_Ctrl <- RunUMAP(PD0_Ctrl, dims = 1:20)

###======Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = list(PD0_Ctrl,PD0_HFD,PD3_Ctrl,PD3_HFD), dims = 1:20)  
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
head(Idents(immune.combined), 5)
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# Determine PCs, Examine and visualize PCA results a few different ways
print(immune.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
DimPlot(immune.combined, reduction = "pca")
# U-MAP and Clustering
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

DimPlot(immune.combined, reduction = "umap", group.by = "tech",pt.size = 0.5)
DimPlot(immune.combined, reduction = "umap", label = TRUE,pt.size = 0.5)


###find markers
DefaultAssay(PF.integration) <- "RNA"
PF.integration.All.Markers <- FindAllMarkers(PF.integration, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25

DefaultAssay(immune.combined) <- "integrated"
###Germ cells 
FeaturePlot(immune.combined, features = c("Ddx4", "Dazl"), min.cutoff = "q9")
###Granulosa cells
FeaturePlot(immune.combined, features = c("Amhr2", "Kitl"), min.cutoff = "q9")
###Epithelial cells
FeaturePlot(immune.combined, features = c("Mfap4", "Col1a1"), min.cutoff = "q9")
###Erythrocytes 
FeaturePlot(immune.combined, features = c("Alas2", "Rhd"), min.cutoff = "q9")
###Endothelial cells 
FeaturePlot(immune.combined, features = c("Egfl7", "Aplnr"), min.cutoff = "q9")
###Immune cells 
FeaturePlot(immune.combined, features = c("Cd53", "Cd52"), min.cutoff = "q9")

####Cell types are defined according to gene marker
Re_name <- RenameIdents(immune.combined, `0` = "...cells", ...., `12` = "...cells")
DimPlot(Re_name, reduction = "umap", label = TRUE,pt.size = 0.5)
saveRDS(Re_name, file = ".../Re_name.rds")


###subset single cell type
cell_type <- subset(Re_name, idents = c("cell_type"))
saveRDS(cell_type, file = "D:\\...\\cell_type.rds")

###Run non-linear dimensional reduction (UMAP/tSNE)
subGerm <- readRDS(file = "D:\\...\\Germ_cells.rds")
subGerm <- FindNeighbors(subGerm, reduction = "pca", dims = 1:20)
subGerm <- FindClusters(subGerm, resolution = 0.5)
subGerm <- RunTSNE(object = subGerm, dims.use = 1:20, do.fast = TRUE)
subGerm <- RunUMAP(subGerm, reduction = "pca", dims = 1:20)
DimPlot(subGerm, reduction = "umap", group.by = "tech",pt.size = 1)
DimPlot(subGerm, reduction = "umap", label = TRUE,pt.size = 1)


#############################################################################################
         #########///////=========Monocle version=2.16========\\\\\\\#############
#############################################################################################

library(monocle)
DefaultAssay(subGerm) <- "RNA"
subGerm <- FindVariableFeatures(subGerm, selection.method = "dispersion")
head(subGerm@meta.data)[1:5, 1:5]

data <- as(as.matrix(subGerm@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = subGerm@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monoGERM <- newCellDataSet(data,
                           phenoData = pd,
                           featureData = fd,
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
HSMM<-monoGERM
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    +                                     num_cells_expressed >= 10))
print(head(pData(HSMM)))

#Clustering cells without marker genes 
disp_table <- dispersionTable(HSMM)
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

#Trajectory step 3: order cells along the trajectory  
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "tech")
plot_cell_trajectory(HSMM, color_by = "State")

###Number of cells per cell state in each group
table(HSMM$tech,HSMM$State)

##########BEAM Function
HSMM <- orderCells(HSMM,root_state = 3)
BEAM_res <- BEAM(HSMM, branch_point = 1)
head(BEAM_res)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
head(BEAM_res)
write.csv(BEAM_res,"BEAM_res.csv")
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
table(BEAM_res$qval < 1e-4)
 ##   FALSE  TRUE 
      17977  3759
BEAM_subset <- subset(BEAM_res,qval < 1e-4)
library("colorRamps")
library("RColorBrewer")
display.brewer.pal(11,"RdYlBu")

my_branched_heatmap <- plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval < 1e-4)),],
                                                   branch_point = 1,
                                                   num_clusters = 4,
                                                   cores = 1,hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(65),
                                                   use_gene_short_name = T,
                                                   show_rownames = F,
                                                   return_heatmap = TRUE)

head(my_branched_heatmap$annotation_row)
dim(my_branched_heatmap$annotation_row)
[1] 4160    1
table(my_branched_heatmap$annotation_row$Cluster)
  1    2    3    4 
1445  729 1354  632 


head(my_branched_heatmap)
branch1 <- my_branched_heatmap$heatmap_matrix
head(branch1)
head(my_branched_heatmap$annotation_row$Cluster)
write.csv(branch1,file = ".../branch1.csv")

my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)
write.csv(my_row,file = ".../myrow.csv")

###Gene Set 1/2/3/4-GO analysis
library(org.Mm.eg.db)
library(clusterProfiler)
x=read.csv("GeneSet1.csv")
eg = bitr(x$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ego <- enrichGO(gene          = eg$ENTREZID,
                #universe      = eg_all$ENTREZID,
                OrgDb         = org.Mm.eg.db,	
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

dotplot(ego,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 5,
        size = NULL,
        split = NULL,
        font.size = 12)
write.table(ego,"GeneSet1-GO.csv", sep = ",", row.names = TRUE)


#############################################################################################
         ######///////=========Desingle software differential analysis========\\\\\\\########
#############################################################################################

### subset Granulosa cells
Granulosa <-subset(Re_name,idents=c("Granulosa cells"))
Granulosa <- readRDS("../Granulose.rds")
DimPlot(Granulosa, reduction = "umap", label = TRUE,pt.size = 0.5)

###DEsignal differential gene analysis
library(DEsingle)
library(Seurat)
subGranulosa <- readRDS(".../Granulosa.rds")
table(subGranulosa$tech,Idents(subGranulosa))
x=read.csv("1.csv")
counts <- as.matrix(subGranulosa@assays$RNA@counts)
colnames(counts)<- x$V1
head(counts)[1:5,1:5]
counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) 
group <- factor(gsub("(PD3_mND|PD3_mHFD).*", "\\1", colnames(counts)), levels = c("PD3_mND", "PD3_mHFD"))
deg<-DEsingle(counts, group, parallel = FALSE)
deglist<- DEtype(results = deg, threshold = 0.05)
write.table(deg,"DEG.csv", sep = ",", row.names = TRUE)
write.table(deglist,"2.csv", sep = ",", row.names = TRUE)

###Volcano Plot  ###Extended Data Fig.6b P0/3_Granulosa
ggplot(data, aes(x = log2FC, y = logP, colour=Status)) +
        geom_point(alpha=0.5, size=1.5) +
        scale_color_manual(values=c("#16982b","#e1a05a", "#DCDCDC"))+
        theme_base()+
        geom_hline(yintercept = -1*log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")

###KEGG
pathway=read.csv("Fig.6f KEGG.csv") 
ggplot(pathway,aes(X,Description))+
  geom_point(aes(size=-log10(qvalue),color=pvalue))+
  scale_color_gradient2(low = "#ff7300",
                        mid = "#d94dff", high = "#4169e1",midpoint =0.45)+
  theme_bw()  
  



######   Metabolome code
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




#####  16s rDNA-seq Code
###PCOA
library(vegan)  
library(ggplot2)
library(plyr)
otu1=read.csv("Fig.4a OTU_PCOA.csv") 
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

###Barplot
rm(list = ls())
phylum <-  read.delim('Fig.4b Phylum_top10.csv', row.names = 1, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
par(xpd = TRUE, mar = par()$mar + c(1, 3, 1, 16))
barplot(as.matrix(100* phylum),
        col = c('#9ab3f5', '#70af85', '#f9e0ae', '#ffabe1', '#c5a880', '#999b84', '#adb36e','#ecb390', '#51adcf', '#79a3b1', 'gray'),
        legend = rownames(phylum),
        cex.axis = 2, cex.names = 2, ylim = c(0, 100), las = 1, width = 0.8, space = 0.6, beside = FALSE,
        args.legend = list(x = 5.5,y=90, 
                           bty = 'n', 
                           inset = -0.18, 
                           cex = 0.5, 
                           y.intersp = 0.8,
                           x.intersp = 0.7, 
                           text.width = 1
        ))
mtext('Relative Abundance(%)', cex = 2, side = 2, line = 4

###heatmap
rm(list = ls())
library(pheatmap)
x=read.csv("Fig.4d Heatmap.csv")
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




#####  ATAC-seq Code
###ATAC-Volcano Plot
library(ggplot2)
data=read.csv("Fig.6b ATAC-Volcano.csv")
head(data)
ggplot(data, aes(x = log2FC, y = logP, colour=Status)) +
        geom_point(alpha=0.5, size=1.5) +
        scale_color_manual(values=c("#16982b","#e1a05a", "#DCDCDC"))+
        theme_base()+
        geom_hline(yintercept = -1*log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")

        
###KEGG
library(ggplot2)
data=read.csv("Fig.6d ATAC-KEGG.csv")
ggplot(pathway,aes(X,Description))+
  geom_point(aes(size=-log10(qvalue),color=pvalue))+
  scale_color_gradient2(low = "#ff7300",
                        mid = "#d94dff", high = "#4169e1",midpoint =0.45)+
  theme_bw()  
