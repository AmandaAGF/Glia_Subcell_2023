library(Seurat)
library(dplyr)

#To manually separate Epithelial and Marginal glia in Adults, we compared them with 
#Davis 2020 FACSorted Epithelial and Marginal glia transcriptomes
glia_adult <- readRDS("glia_adult.rds")
#Converted to match our nomenclature
davis <- read.csv("Full_davis_data_converted.csv")

rownames(davis) <- davis[,1]
davis <- davis[,16:17]

Seurat_D_data = CreateSeuratObject(counts = davis)
Seurat_D_data = NormalizeData(Seurat_D_data)

#Subset Ep/Mg cluster
Ep_Mg <- subset(glia_adult, idents = "Ep/Mg") #or change Ep/Mg for "14" if not annotated yet
Ep_Mg <- FindVariableFeatures(Ep_Mg, selection.method = "vst")
Ep_Mg = ScaleData(Ep_Mg)
Ep_Mg_data <- as.matrix(Ep_Mg@assays$RNA@data)

#Define the genes that are markers of Ep/Mg cluster
markers <- FindMarkers(glia_adult, assay = "RNA", ident.1 = "Ep/Mg") #or change Ep/Mg for "14" if not annotated yet
genes = rownames(markers)
genes = make.names(genes)

#Intersect the genes that are markers with the ones present in Davis data
genes_in_common_D = intersect(rownames(davis), genes)

#Make correlation between each cell and Davis data
cor_Davis_scRNAseq = cor(as.matrix(Ep_Mg_data)[genes_in_common_D,],
                         as.matrix(Seurat_D_data@assays$RNA@data)[genes_in_common_D,])


Epithelial <- rownames(cor_Davis_scRNAseq)[which(cor_Davis_scRNAseq[,1] > cor_Davis_scRNAseq[,2])]
Marginal <- rownames(cor_Davis_scRNAseq)[which(cor_Davis_scRNAseq[,1] < cor_Davis_scRNAseq[,2])]



#Add new IDs to glia_adult
Idents(glia_adult) <- glia_adult$seurat_clusters

glia_adult@meta.data$final.IDs <- ""

for(i in 0:20) {
  names_cells <- WhichCells(subset(glia_adult, idents = i))
  index <- which(colnames(glia_adult) %in% names_cells)
  glia_adult@meta.data$final.IDs[index] <- new.cluster.id[i+1]
}

glia_adult$final.IDs[which(colnames(glia_adult) %in% Epithelial)] <- "Epithelial"
glia_adult$final.IDs[which(colnames(glia_adult) %in% Marginal)] <- "Marginal"
Idents(glia_adult) <- glia_adult$final.IDs
DimPlot(glia_adult, label = T) + NoLegend()



