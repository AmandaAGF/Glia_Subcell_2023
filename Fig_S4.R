neuropil <- readRDS("medulla_neuropil_allstages.rds")
neuropil_1 <- readRDS("neuropilX_allstages.rds")

#Only adult Medulla Neuropil
Idents(neuropil) <- neuropil$Clusters_Adult
tiff(filename = "Fig_S4/medulla_neuropil_DimPlot_Adult.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(neuropil, order = T, pt.size = 0.75 ,cols = c("#00B9DD", "#00B9DD","#00A5FF", "#cccccc") ,na.value = "#cccccc") + NoLegend()
dev.off()

#Only adult X Neuropil
Idents(neuropil_1) <- neuropil_1$Clusters_Adult
tiff(filename = "Fig_S4/neuropilX_DimPlot_Adult.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(neuropil_1, order = T, pt.size = 0.75,cols = c("#c1aed4" ,"#8e809c","#524b59", "#d7c5e8","#677a6e","#a3c2ae","#e6e3e3"),na.value = "#e6e3e3") + NoLegend()
dev.off()


#Genes expressed in early neuropil glia
all_glia@active.assay = "RNA"
tiff(filename = "Fig_S4/all_glia_FeaturePlot_neuropilmarkers.tiff",width = 16, height = 10,units = "in", res = 400)
FeaturePlot(all_glia, features = c("CR30009", "gcm","E(spl)mgamma-HLH", "CG11670", "E(spl)malpha-BFM", "Traf4"), order = T, ncol = 3)
dev.off()

#Genes expressed in early neuropil glia - plot with only neuropil (except medulla and lamina)
neuropil_1@active.assay = "RNA"
tiff(filename = "Fig_S4/neuropil_1_FeaturePlot_neuropilmarkers.tiff",width = 16, height = 5,units = "in", res = 400)
FeaturePlot(neuropil_1, features = c("CR30009", "gcm","E(spl)malpha-BFM"), order = T, ncol = 3)
dev.off()


