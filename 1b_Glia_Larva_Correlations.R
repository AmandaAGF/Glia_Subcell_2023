library(Seurat)
library(dplyr)
library(ggplot2)


#Correlation between the Larval clusters and the Repo and Elav sorted cells transcriptomes
#Load Larva from Konstatinides 2022
larva <- readRDS("larvaOL.integrated150.rds")

#Load GSE142788_NormalizedBulk.csv.gz but with results combined for triplicates and genes names converted to match the genome annotation used in the scRNAseq
FACS_Bulk <- read.csv("Full_Katarina_data_converted.csv")

rownames(FACS_Bulk) <- FACS_Bulk[,1]

#Select only elav and repo
FACS_Bulk <- FACS_Bulk[,22:23]
colnames(FACS_Bulk) <- c("Elav", "Repo")

Seurat_D_data = CreateSeuratObject(counts = FACS_Bulk)
Seurat_D_data = NormalizeData(Seurat_D_data)

scRNAseq_expression = AverageExpression(larva, assay = 'RNA', return.seurat = T)

larva_markers <- read.table("Markers_Glia_Larva.txt")
larva_markers <- FindAllMarkers(larva, only.pos = T, assay = "RNA", logfc.threshold = 0.5)
top10 <- larva_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
genes = unique(top10$gene)
genes = make.names(genes)


genes_in_common_D = intersect(rownames(FACS_Bulk), genes)

#Make Pearson correlation matrix to compare transcriptome of larval clusters with the bulk seq repo and elav
cor_FACS_Bulk_scRNAseq = cor(as.matrix(scRNAseq_expression@assays$RNA@data)[genes_in_common_D,],
                            as.matrix(Seurat_D_data@assays$RNA@data)[genes_in_common_D,])

cor_FACS_Bulk_scRNAseq <- as.data.frame(cor_FACS_Bulk_scRNAseq)
cor_FACS_Bulk_scRNAseq$label <- rownames(cor_FACS_Bulk_scRNAseq)

cor_FACS_Bulk_scRNAseq <- cor_FACS_Bulk_scRNAseq[order(-cor_FACS_Bulk_scRNAseq$Repo),]

cor_FACS_Bulk_scRNAseq$identity <- "Non-Glial Cluster"
cor_FACS_Bulk_scRNAseq$identity[which(cor_FACS_Bulk_scRNAseq$Repo > 0.2 & cor_FACS_Bulk_scRNAseq$Elav < 0.2)] <- "Glial Cluster"


#Plot
ggplot(data = cor_FACS_Bulk_scRNAseq, mapping = aes(x = Elav, y = Repo, color = identity, label = label)) + 
  geom_point() +  scale_color_manual(values = c("red", "black")) +  
  geom_text(data = cor_FACS_Bulk_scRNAseq[which(cor_FACS_Bulk_scRNAseq$identity == "Glial Cluster"),],aes(label=label),hjust=0, vjust=0, size = 6) +
  theme(panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size = 20, face = "bold"), plot.title = element_text(size = 20), axis.text.y = element_text(size = 15))
