library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(ggpubr)


#Read Adult data
adult = readRDS("GSE142787_Adult.rds")

adult@meta.data$Current_Idents = adult@active.ident

#Add annotation from TableS1 from Ozel 2020
TableS1 <- read_excel("ClusterAnnotationTable.xlsx")
TableS1 = TableS1[,1:2]
colnames(TableS1)[2] = "Annotations"
Annotations = match(as.character(Idents(adult)), TableS1$`Cluster number`)
Annotations = TableS1$Annotations[Annotations]
Idents(adult) = Annotations
adult@meta.data$Annotations = adult@active.ident


# Select clusters that were previously identified as Glia and/or Glia/LQ but not LQ
glia_adult <- subset(adult, idents = c("G04", "G08", "G03", "G09", "G06", "G/LQ3", "G02",
                                       "G05", "G14ab", "G/LQ1", "G11", "G/LQ4", "G10", "G07",
                                       "G01", "G/LQ2", "G12", "G13", "G16/LQ"))

glia_adult <- FindVariableFeatures(glia_adult, selection.method = "vst")
glia_adult = ScaleData(glia_adult, assay = "RNA")
glia_adult <- RunPCA(glia_adult, npcs = 100, verbose = FALSE)
glia_adult <- RunUMAP(glia_adult, reduction = "pca", dims = 1:20, assay = "RNA")
glia_adult <- FindNeighbors(glia_adult, reduction = "pca", dims = 1:20, assay = "RNA")
glia_adult <- FindClusters(glia_adult, resolution = 0.6)
DimPlot(glia_adult, label = T) + NoLegend()


markers <- FindAllMarkers(glia_adult, assay = "RNA", logfc.threshold = 1)

glia_adult <- BuildClusterTree(glia_adult, assay = "integrated", 
                               features = unique(markers$gene[which(markers$avg_log2FC > 2)]))
PlotClusterTree(glia_adult, no.margin = T, label.offset = 20, direction = "rightwards")


#saveRDS(glia_adult, "glia_adult.rds")


#Rename for publication
new.cluster.id <- c("Astrocytes_4", "Astrocytes_3", "Perineurial", "Medulla_ALG", "Ensheathing_1", "Astrocytes_1", "Ensheathing_2", "Distal_Satellite", "Inner_Chiasm", "Medulla_EG_1", "Chalice", "Cortex_1", "Pseudocartridge", "Astrocytes_2", "Ep/Mg", "Outer_Chiasm", "Fenestrated", "Cortex_2", "Medulla_EG_2", "Proximal_Satellite", "Subperineurial")
names(new.cluster.id) <- levels(glia_adult)
glia_adult <- RenameIdents(glia_adult, new.cluster.id)
glia_adult@meta.data$new.cluster.id <- glia_adult@active.ident
DimPlot(glia_adult, label = T, group.by= "new.cluster.id")+NoLegend()

#Add All_TFs to the metadata
table_All_TFs <- read.table("GO_Term_Tables/GO_All_TFs_Drosophila.tsv")
list_table_All_TFs <- table_All_TFs[,2]
all_gene_names <- rownames(glia_adult)
All_TFs_in_the_dataset <- intersect(list_table_All_TFs, all_gene_names)
glia_adult@meta.data$All_TFs <- PercentageFeatureSet(glia_adult, features = All_TFs_in_the_dataset, assay = "RNA")
glia_adult@meta.data$All_TFs <- glia_adult@meta.data$All_TFs[,1]


#Separate Epithelial and Marginal
#We determined which one is Epithelial or Marginal based on the script 2b_Glia_Adult_Corr_Epi_Marg


#saveRDS(glia_adult, "glia_adult_1.rds")
glia_adult <- readRDS("glia_adult_1.rds")


#Overcluster to test for processes
glia_adult <- FindClusters(glia_adult, resolution = 2)

#Run function to test for processes
TestAllProcesses(glia_adult)
TestAllProcesses_Correlation(glia_adult)

#Remove processes from dataset
Idents(glia_adult) <- glia_adult$final.IDs
glia_adult <- subset(glia_adult, idents = c("Astrocytes_3", "Perineurial", "Medulla_ALG", "Ensheathing_1", "Ensheathing_2", "Distal_Satellite", "Inner_Chiasm", "Chalice", "Pseudocartridge", "Epithelial", "Marginal", "Outer_Chiasm", "Fenestrated", "Cortex_2", "Medulla_EG_2", "Proximal_Satellite", "Subperineurial"))
Idents(glia_adult) <- glia_adult$integrated_snn_res.2
glia_adult <- subset(glia_adult, idents = c("26", "21", "8", "30", "20", "14", "13", "16","24", "5", "18", "23", "1", "17", "25", "29", "12", "4", "31", "28"))


glia_adult <- FindVariableFeatures(glia_adult, selection.method = "vst")
glia_adult = ScaleData(glia_adult, assay = "RNA")
glia_adult <- RunPCA(glia_adult, npcs = 100, verbose = FALSE)
glia_adult <- RunTSNE(glia_adult, reduction = "pca", dims = 1:15, assay = "RNA")
glia_adult <- FindNeighbors(glia_adult, reduction = "pca", dims = 1:15, assay = "RNA")
Idents(glia_adult) <- glia_adult$final.IDs
DimPlot(glia_adult, label = T, reduction = "tsne") + NoLegend()

#saveRDS(glia_adult, "glia_adult_final.rds")
glia_adult <- readRDS("glia_adult_final.rds")


#Get Markers for each cluster
DEG_glia_adult <- FindAllMarkers(glia_adult, only.pos = T, assay = "RNA")
write.table(DEG_glia_adult, "Markers_Glia_Adult.txt", sep = "\t", quote = F)
