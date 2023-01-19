glias_larva <- readRDS("glias_larva.rds")
glia_adult <- readRDS("glia_adult.rds")

##### Fig. 2
#Larva
tiff(filename = "Fig_2/glias_larva_DimPlot.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glias_larva, pt.size = 0.1, label = T, label.size = 5) + NoLegend()
dev.off()

tiff(filename = "Fig_2/glias_larva_ClusterTree.tiff",width = 15, height = 15,units = "cm", res = 400)
PlotClusterTree(glias_larva, no.margin = T, label.offset = 20, direction = "rightwards", edge.width = 3, font = "arial")
dev.off()

#Adult
tiff(filename = "Fig_2/glia_adult_DimPlot.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glia_adult, label = T, pt.size = 0.1, label.size = 5) + NoLegend()
dev.off()

tiff(filename = "Fig_2/glia_adult_ClusterTree.tiff",width = 15, height = 15,units = "cm", res = 400)
PlotClusterTree(glia_adult, no.margin = T, label.offset = 20, direction = "rightwards", edge.width = 3, font = "arial")
dev.off()
