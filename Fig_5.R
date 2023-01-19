glias_larva <- readRDS("glias_larva_final.rds")
glia_adult <- readRDS("glia_adult_final.rds")

##### Fig. 2
#Larva
tiff(filename = "Fig_5/glias_larva_DimPlot.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glias_larva, cols = c("#00A7FF", "#C77CFF", "#FE73D2", "#64B100", 
                              "#02BE67" ,"#F55442", "#00BFC4", "#E7861B",
                              "#395FFA", "#C59900"), pt.size = 0.1) + NoLegend()
dev.off()


#Adult
tiff(filename = "Fig_5/glia_adult_DimPlot.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glia_adult, cols = c("#524B59", "#AA9E00", "#64AE00", "#00B21B", 
                             "#64AA00", "#A3C2AE", "#677A6E", "#02B767", 
                             "#00B4D8", "#00D3DA", "#009FFF", "#BE7CFF", 
                             "#AC80FF", "#E0801B", "#BC9000", "#F573C9", 
                             "#395FF3"), pt.size = 0.1, reduction = "tsne") + NoLegend()
dev.off()
