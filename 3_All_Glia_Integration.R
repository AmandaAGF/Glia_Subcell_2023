library(Seurat)

#Load glia larva from script 2_Glia_Larva
glias_larva <- readRDS("glias_larva_2.rds")

#All the identified glial clusters from Ozel 2020 across all stages
b <- c("198","212","210", "201","197","207","203","196","213","219","204","200","202", 
       "185", "186", "187", "188", "189", "190", "199", "205", "206", "208", "209", "211", "216", "217")

#Load P30 data from Ozel 2020
P30_Glia = readRDS("GSE142787_P30.rds")
P30_Glia@active.assay <- "RNA"
P30_Glia@meta.data$Clusters_P30 <- P30_Glia@active.ident
P30_Glia@meta.data$Clusters_ozel2020 <- P30_Glia@active.ident

#Will select the clusters that match the glial clusters
c <- levels(P30_Glia)
d <- intersect(b,c)
P30_Glia <- subset(P30_Glia, idents = d)
P30_Glia <- FindVariableFeatures(P30_Glia, selection.method = "vst")
P30_Glia = ScaleData(P30_Glia)
P30_Glia <- RunPCA(P30_Glia, npcs = 25, verbose = FALSE)
P30_Glia <- RunUMAP(P30_Glia, reduction = "pca", dims = 1:20)
P30_Glia <- FindNeighbors(P30_Glia, reduction = "pca", dims = 1:20)

#Load P15 data from Ozel 2020
P15_Glia = readRDS("GSE142787_P15.rds")
P15_Glia@active.assay <- "RNA"
P15_Glia@meta.data$Clusters_P15 <- P15_Glia@active.ident
P15_Glia@meta.data$Clusters_ozel2020 <- P15_Glia@active.ident

#Selected all the clusters that were annotated as glia in adults that are present in the P15 dataset except 206
e <- b[-22]
c <- levels(P15_Glia)
d <- intersect(e,c)
P15_Glia <- subset(P15_Glia, idents = d)
P15_Glia <- FindVariableFeatures(P15_Glia, selection.method = "vst")
P15_Glia = ScaleData(P15_Glia)
P15_Glia <- RunPCA(P15_Glia, npcs = 25, verbose = FALSE)
P15_Glia <- RunUMAP(P15_Glia, reduction = "pca", dims = 1:15)
P15_Glia <- FindNeighbors(P15_Glia, reduction = "pca", dims = 1:15)

#Load P50 data from Ozel 2020
P50_Glia = readRDS("GSE142787_P50.rds")
P50_Glia@active.assay <- "RNA"
P50_Glia@meta.data$Clusters_P50 <- P50_Glia@active.ident
P50_Glia@meta.data$Clusters_ozel2020 <- P50_Glia@active.ident

#Selected all the clusters that were annotated as glia in adults that are present in the P50 dataset
c <- levels(P50_Glia)
d <- intersect(b,c)
P50_Glia <- subset(P50_Glia, idents = d)
P50_Glia <- FindVariableFeatures(P50_Glia, selection.method = "vst")
P50_Glia = ScaleData(P50_Glia)
P50_Glia <- RunPCA(P50_Glia, npcs = 15, verbose = FALSE)
P50_Glia <- RunUMAP(P50_Glia, reduction = "pca", dims = 1:15)
P50_Glia <- FindNeighbors(P50_Glia, reduction = "pca", dims = 1:15)

#Load P70 data from Ozel 2020
P70_Glia = readRDS("GSE142787_P70.rds")
P70_Glia@active.assay <- "RNA"
P70_Glia@meta.data$Clusters_P70 <- P70_Glia@active.ident
P70_Glia@meta.data$Clusters_ozel2020 <- P70_Glia@active.ident

#Selected all the clusters that were annotated as glia in adults that are present in the P70 dataset
c <- levels(P70_Glia)
d <- intersect(b,c)
P70_Glia <- subset(P70_Glia, idents = d)
P70_Glia <- FindVariableFeatures(P70_Glia, selection.method = "vst")
P70_Glia = ScaleData(P70_Glia)
P70_Glia <- RunPCA(P70_Glia, npcs = 20, verbose = FALSE)
P70_Glia <- RunUMAP(P70_Glia, reduction = "pca", dims = 1:15)
P70_Glia <- FindNeighbors(P70_Glia, reduction = "pca", dims = 1:15)

#Load Adult data from script 2_Glia_Adult
Adult_Glia <- readRDS("glia_adult_1.rds")
Adult_Glia@meta.data$Clusters_Adult <- Adult_Glia@active.ident
Idents(Adult_Glia) <- Adult_Glia$Clustering
Adult_Glia@meta.data$Clusters_ozel2020 <- Adult_Glia@active.ident

##Integrating both L3, P15, P30, P50, P70 and Adults glia datasets
#Looking for anchors for the two datasets
glia.anchors <- FindIntegrationAnchors(object.list = list(glias_larva, P15_Glia, P30_Glia, P50_Glia, P70_Glia, Adult_Glia), dims = 1:30, assay = c("RNA", "RNA", "RNA", "RNA", "RNA", "RNA"))

#Combine the dataset
all_glia <- IntegrateData(glia.anchors, dims = 1:30)
all_glia <- FindVariableFeatures(all_glia, selection.method = "vst")
all_glia = ScaleData(all_glia)
all_glia <- RunPCA(all_glia, npcs = 50, verbose = FALSE)
ElbowPlot(all_glia, ndims = 50) # To determine how many pcs I am going to use in the next steps
all_glia <- RunUMAP(all_glia, reduction = "pca", dims = 1:20, assay = "RNA")
all_glia <- FindNeighbors(all_glia, reduction = "pca", dims = 1:20)
all_glia <- FindClusters(all_glia, resolution = 0.7)
DimPlot(all_glia, label = T) + NoLegend()
DimPlot(all_glia, label = T, group.by = "Clusters_Adult", order = T) + NoLegend()
DimPlot(all_glia, label = T, group.by = "Clusters_L3", order = T) + NoLegend()

#Label each stage
index_L3 <- which(all_glia@meta.data$Clusters_L3 != "NA")
index_P15 <- which(all_glia@meta.data$Clusters_P15 != "NA")
index_P50 <- which(all_glia@meta.data$Clusters_P50 != "NA")
index_P70 <- which(all_glia@meta.data$Clusters_P70 != "NA")
index_adult <- which(all_glia@meta.data$Clusters_Adult != "NA")
all_glia@meta.data$stage <- "P30"
all_glia@meta.data$stage[index_L3] <- "Larva"
all_glia@meta.data$stage[index_P15] <- "P15"
all_glia@meta.data$stage[index_P50] <- "P50"
all_glia@meta.data$stage[index_P70] <- "P70"
all_glia@meta.data$stage[index_adult] <- "Adult"
DimPlot(all_glia, group.by = "stage")

Idents(all_glia) <- all_glia$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
DimPlot(all_glia, cols = colfunc(7))


#exclude cluster "4" which is formed by many different clusters and expresses neuronal markers
before_exclude <- all_glia
Idents(all_glia) <- all_glia$integrated_snn_res.0.7
a <- levels(all_glia)
a <- a[-5]
all_glia <- subset(all_glia, idents = a)
all_glia <- FindVariableFeatures(all_glia, selection.method = "vst")
all_glia = ScaleData(all_glia)
all_glia <- RunPCA(all_glia, npcs = 25, verbose = FALSE)
all_glia <- RunUMAP(all_glia, reduction = "pca", dims = 1:14, assay = "RNA")
all_glia <- FindNeighbors(all_glia, reduction = "pca", dims = 1:14)
#all_glia <- FindClusters(all_glia, resolution = 1.2)
DimPlot(all_glia, label = T) + NoLegend()
DimPlot(all_glia, label = T, group.by = "Clusters_L3", order = T) + NoLegend()
DimPlot(all_glia, label = T, group.by = "Clusters_Adult", order = T) + NoLegend()

#saveRDS(all_glia, file = "all_glia_final.rds")
all_glia <- readRDS("all_glia_final.rds")





#Separating Neuropil glia
#Medulla neuropil
Idents(all_glia) <- all_glia$integrated_snn_res.0.7
neuropil <- subset(all_glia, idents = c("27", "23", "5", "10"))
neuropil <- FindVariableFeatures(neuropil, selection.method = "vst")
neuropil = ScaleData(neuropil)
neuropil <- RunPCA(neuropil, npcs = 25, verbose = FALSE)
ElbowPlot(neuropil, ndims = 100) # To determine how many pcs I am going to use in the next steps
neuropil <- RunUMAP(neuropil, reduction = "pca", dims = 1:14, assay = "RNA")
neuropil <- FindNeighbors(neuropil, reduction = "pca", dims = 1:14)
DimPlot(neuropil, label = T) + NoLegend()

Idents(neuropil) <- neuropil$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
DimPlot(neuropil, cols = colfunc(6)) +NoLegend()
DimPlot(neuropil, group.by = "Clusters_Adult", label= T) + NoLegend()

#Remove label from contaminating cells
Idents(neuropil) <- neuropil$Clusters_Adult
a <- levels(neuropil)
index <- which(neuropil@active.ident != "Medulla_ALG" & neuropil@active.ident != "Medulla_EG_1" & neuropil@active.ident != "Medulla_EG_2")
neuropil$Clusters_Adult[index] <- NA


Idents(neuropil) <- neuropil$Clusters_Adult
DimPlot(neuropil, order = T, cols = c("#00B1D5" ,"#00B1D5" ,"#009DFF" ,"#cccccc"),na.value = "#cccccc") + NoLegend()
#saveRDS(neuropil, file = "medulla_neuropil_allstages.rds")



#Neuropil Glia
Idents(all_glia) <- all_glia$integrated_snn_res.0.7
neuropil_1 <- subset(all_glia, idents = c("0", "24", "2"))
neuropil_1 <- FindVariableFeatures(neuropil_1, selection.method = "vst")
neuropil_1 = ScaleData(neuropil_1)
neuropil_1 <- RunPCA(neuropil_1, npcs = 25, verbose = FALSE)
ElbowPlot(neuropil_1, ndims = 100) # To determine how many pcs I am going to use in the next steps
neuropil_1 <- RunUMAP(neuropil_1, reduction = "pca", dims = 1:15, assay = "RNA")
neuropil_1 <- FindNeighbors(neuropil_1, reduction = "pca", dims = 1:15)
DimPlot(neuropil_1, label = T) + NoLegend()

Idents(neuropil_1) <- neuropil_1$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
DimPlot(neuropil_1, cols = colfunc(6)) +NoLegend()
DimPlot(neuropil_1, group.by = "Clusters_Adult", label= T) + NoLegend()

Idents(neuropil_1) <- neuropil_1$Clusters_Adult
a <- levels(neuropil_1)
index <- which(neuropil_1@active.ident != "Astrocytes_1" & neuropil_1@active.ident != "Astrocytes_4" &
                 neuropil_1@active.ident != "Astrocytes_2" & neuropil_1@active.ident != "Astrocytes_3" &
                 neuropil_1@active.ident != "Ensheathing_1" & neuropil_1@active.ident != "Surface_1")
neuropil_1$Clusters_Adult[index] <- NA
Idents(neuropil_1) <- neuropil_1$Clusters_Adult


DimPlot(neuropil_1, order = T, cols = c("#c1aed4" ,"#8e809c","#524b59", "#d7c5e8","#677a6e","#a3c2ae","#cccccc"),na.value = "#cccccc") + NoLegend()
#saveRDS(neuropil_1, file = "neuropilX_allstages.rds")



#Genes expressed in early neuropil glia
early_markers <- FindMarkers(all_glia, ident.1 = c("10", "14"), only.pos = T, assay = "RNA", 
                             logfc.threshold = 1, min.diff.pct = 0.4)

write.table(early_markers, "Early_Markers_Neuropil_AllGlia.txt", sep = "\t", quote = F)




