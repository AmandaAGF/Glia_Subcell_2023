library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(xlsx)

#Read larva data from Konstantinides 2022
larva = readRDS("larvaOL.integrated150.rds")
larva@meta.data$Clusters_Konstantinides2022 <- larva@meta.data$integrated_snn_res.2

#Select clusters identified in 2b_Glia_Larva_Correlations
glias_larva <- subset(larva, idents = c("44", "87", "45", "56", "67", "91", "88", "63", "85", "2", "3"))
glias_larva <- FindVariableFeatures(glias_larva, selection.method = "vst")
glias_larva = ScaleData(glias_larva)
glias_larva <- RunPCA(glias_larva, npcs = 100, verbose = FALSE)
glias_larva <- RunUMAP(glias_larva, reduction = "pca", dims = 1:30, assay = "RNA")
glias_larva <- FindNeighbors(glias_larva, reduction = "pca", dims = 1:30)
glias_larva <- FindClusters(glias_larva, resolution = 0.8)
DimPlot(glias_larva, label = T) + NoLegend()


#Hierarchical Clustering Tree
glias_larva <- BuildClusterTree(glias_larva, assay = "integrated")
PlotClusterTree(glias_larva, no.margin = T, label.offset = 20, direction = "rightwards")

#saveRDS(glias_larva, "glias_larva.rds")
glias_larva <- readRDS("glias_larva.rds")


#Make a new column in the metadata that contains the percentage of All TFs genes extracted from flymine.org GO == "Transcription Factor"
table_All_TFs <- read.table("GO_All_TFs_Drosophila.tsv")
list_table_All_TFs <- table_All_TFs[,2]
all_gene_names <- rownames(glias_larva@assays$RNA@counts)
All_TFs_in_the_dataset <- intersect(list_table_All_TFs, all_gene_names)
glias_larva@meta.data$All_TFs <- PercentageFeatureSet(glias_larva, features = All_TFs_in_the_dataset, assay = "RNA")
glias_larva@meta.data$All_TFs <- glias_larva@meta.data$All_TFs[,1]


#IDs after validation for publication
new.cluster.id <- c("Cortex_Central_Brain", "Cortex", "Cortex_Central_Brain", "Cortex", "Chiasm", "Cortex_Central_Brain", "Medulla_Neuropil", "Cortex", "Perineurial", "Medulla_Neuropil", "Epithelial_Marginal" ,"LobulaC_Neuropil", "Chiasm_1", "Central_Brain", "Perineurial", "Epithelial_Marginal", "Satellite", "Central_Brain", "Cortex_Central_Brain", "Retina_Wrapping", "Cortex", "Subperineurial", "Medulla_Neuropil")
names(new.cluster.id) <- levels(glias_larva)
glias_larva <- RenameIdents(glias_larva, new.cluster.id)
glias_larva@meta.data$new.cluster.id <- glias_larva@active.ident
DimPlot(glias_larva, label = T) + NoLegend()

#Separate Xgo and Xgi
Chiasm <- subset(glias_larva, idents = "Chiasm_1")
Chiasm <- FindVariableFeatures(Chiasm, selection.method = "vst")
Chiasm = ScaleData(Chiasm)
Chiasm <- RunUMAP(Chiasm, reduction = "pca", dims = 1:25)
Chiasm <- FindNeighbors(Chiasm, reduction = "pca", dims = 1:25)
Chiasm <- FindClusters(Chiasm, resolution = 0.2)
DimPlot(Chiasm)
Xgo <- WhichCells(subset(Chiasm, idents = "0"))
Xgi <- WhichCells(subset(Chiasm, idents = "1"))

#Add new IDs to glias_larva
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8

glias_larva@meta.data$Clusters_L3 <- ""

for(i in 0:22) {
  names_cells <- WhichCells(subset(glias_larva, idents = i))
  index <- which(colnames(glias_larva) %in% names_cells)
  glias_larva@meta.data$Clusters_L3[index] <- new.cluster.id[i+1]
}

index_Xgo <- which(colnames(glias_larva) %in% Xgo)
index_Xgi <- which(colnames(glias_larva) %in% Xgi)
glias_larva$Clusters_L3[index_Xgi] <- "Inner_Chiasm"
glias_larva$Clusters_L3[index_Xgo] <- "Outer_Chiasm"
Idents(glias_larva) <- glias_larva$Clusters_L3
DimPlot(glias_larva, label = T) + NoLegend()


#saveRDS(glias_larva, "glias_larva_1.rds")
glias_larva <- readRDS("glias_larva_1.rds")

#Subset to exclude central brain
glias_larva <- subset(glias_larva, idents = c("Cortex", "Chiasm", "Medulla_Neuropil", "Perineurial", "Epithelial_Marginal" ,"LobulaC_Neuropil", "Inner_Chiasm", "Outer_Chiasm","Perineurial", "Satellite", "Retina_Wrapping", "Subperineurial"))

glias_larva <- FindVariableFeatures(glias_larva, selection.method = "vst")
glias_larva = ScaleData(glias_larva)
glias_larva <- RunPCA(glias_larva, npcs = 100, verbose = FALSE)
glias_larva <- RunUMAP(glias_larva, reduction = "pca", dims = 1:15, assay = "RNA")
glias_larva <- FindNeighbors(glias_larva, reduction = "pca", dims = 1:15)
DimPlot(glias_larva, label = T) + NoLegend()


#saveRDS(glias_larva, "glias_larva_2.rds")
glias_larva <- readRDS("glias_larva_2.rds")


#Separate Chiasm and Cortex PC and CB
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8

new.cluster.id_1 <- c("Cortex_PC", "Cortex_PC", "Chiasm_PC", "Medulla_Neuropil", "Cortex_PC", "Perineurial", "Medulla_Neuropil", "Epithelial_Marginal" ,"LobulaC_Neuropil", "Chiasm_CB", "Perineurial", "Epithelial_Marginal", "Satellite","Retina_Wrapping", "Cortex_CB", "Subperineurial", "Medulla_Neuropil")
names(new.cluster.id_1) <- levels(glias_larva)
glias_larva <- RenameIdents(glias_larva, new.cluster.id_1)
glias_larva@meta.data$new.cluster.id_1 <- glias_larva@active.ident
DimPlot(glias_larva, label = T) + NoLegend()

#Get markers of Chiasm_PC and Cortex_PC
Idents(glias_larva) <- glias_larva$new.cluster.id_1
Markers_Chiasm_PC <- FindMarkers(glias_larva, ident.1 = "Chiasm_PC", ident.2 = "Chiasm_CB", only.pos = T, assay = "RNA")
write.table(Markers_Chiasm_PC, "Markers_Chiasm_Processes_L3.txt", sep = "\t", quote = F)

Markers_Cortex_PC <- FindMarkers(glias_larva, ident.1 = "Cortex_PC", ident.2 = "Cortex_CB", only.pos = T, assay = "RNA")
write.table(Markers_Cortex_PC, "Markers_Cortex_Processes_L3.txt", sep = "\t", quote = F)



#Separate Glia vs No-Glia in Larva - Also in Fig. S1
index <- WhichCells(larva, idents = c("44", "87", "45", "56", "67", "91", "88", "63", "85", "2", "3"))
index_1 <- which(colnames(larva) %in% index)
larva$Glia <- "Non-Glial Cells"
larva$Glia[index_1] <- "Glia"
DimPlot(larva, group.by = "Glia")
Idents(larva) <- larva$Glia

#Get markers Glia vs Non-Glia
markers_glia <- FindMarkers(larva, ident.1 = "Glia", ident.2 = "Non-Glial Cells", only.pos = T, assay = "RNA", min.diff.pct = 0.5, logfc.threshold = 1)
write.table(markers_glia, "Markers_Glia_L3.txt", sep = "\t", quote = F)



#Remove processes from dataset
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8
TestAllProcesses(glias_larva)
TestAllProcesses_Correlation(glias_larva)

glias_larva <- subset(glias_larva, idents = c("20", "12", "19", "14", "11", "15", "10", "16", "22", "9", "6", "21"))

glias_larva <- FindVariableFeatures(glias_larva, selection.method = "vst")
glias_larva = ScaleData(glias_larva)
glias_larva <- RunPCA(glias_larva, npcs = 100, verbose = FALSE)
glias_larva <- RunUMAP(glias_larva, reduction = "pca", dims = 1:15, assay = "RNA")
glias_larva <- FindNeighbors(glias_larva, reduction = "pca", dims = 1:15)
Idents(glias_larva) <- glias_larva$Clusters_L3
DimPlot(glias_larva, label = T) + NoLegend()

#saveRDS(glias_larva, "glias_larva_final.rds")
glias_larva <- readRDS("glias_larva_final.rds")



#Get Markers for each cluster
DEG_glias_larva <- FindAllMarkers(glias_larva, only.pos = T, assay = "RNA")
write.table(DEG_glias_larva, "Markers_Glia_Larva.txt", sep = "\t", quote = F)



