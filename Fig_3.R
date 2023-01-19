library(Seurat)
library(ggplot2)

all_glia <- readRDS("all_glia_final.rds")
neuropil <- readRDS("medulla_neuropil_allstages.rds")
neuropil_1 <- readRDS("neuropilX_allstages.rds")

#Fig. 3
#Small Dimplots - After run 3_All_Glia_Integration
tiff(filename = "Fig_3/L3_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(glias_larva) + NoLegend() + NoAxes()
dev.off()
tiff(filename = "Fig_3/P15_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(P15_Glia) + NoLegend() + NoAxes()
dev.off()
tiff(filename = "Fig_3/P30_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(P30_Glia) + NoLegend() + NoAxes()
dev.off()
tiff(filename = "Fig_3/P50_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(P50_Glia) + NoLegend() + NoAxes()
dev.off()
tiff(filename = "Fig_3/P70_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(P70_Glia) + NoLegend() + NoAxes()
dev.off()
tiff(filename = "Fig_3/Adult_miniDimPlot.tiff",width = 2, height = 2,units = "in", res = 400)
DimPlot(Adult_Glia) + NoLegend() + NoAxes()
dev.off()


#Dimplot all stages
Idents(all_glia) <- all_glia$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
tiff(filename = "Fig_3/all_glia_DimPlot_stages.tiff",width = 8, height = 6.5,units = "in", res = 400)
DimPlot(all_glia, cols = colfunc(7)) +theme(legend.text = element_text(size=20, face = "bold"))
dev.off()

##Only Adult
Idents(all_glia) <- all_glia$Clusters_Adult
levels(all_glia) <- levels(glia_adult)
tiff(filename = "Fig_3/all_glia_DimPlot_Adult.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(all_glia, order = T,pt.size = 0.75 ,cols = c("#02BC67", "#c1aed4", "#8e809c", "#524b59", 
                                                     "#ADA100", "#d7c5e8", "#64B100", "#00B81B", 
                                                     "#64B100", "#a3c2ae", "#C57CFF", "#677a6e", 
                                                     "#00BADE", "#00BADE", "#00A5FF", "#00D9E0", 
                                                     "#395FFA", "#B284FF", "#E7861B", "#C39700", 
                                                     "#FC73D0", "#395FFA", "#cccccc") ,na.value = "#cccccc") + NoLegend()
dev.off()

##Only L3
Idents(all_glia) <- all_glia$Clusters_L3
tiff(filename = "Fig_3/all_glia_DimPlot_L3.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(all_glia, order = T, pt.size = 0.75,cols = c("#395FFA", "#00A7FF", "#BB9D00", "#C77CFF", 
                                                     "#FE73D2" ,"#64B100", "#02BE67", "#F55442", 
                                                     "#00BFC4", "#E7861B", "#C59900", "#cccccc") ,na.value = "#cccccc") + NoLegend()
dev.off()

#For Madulla Neuropil Glia
##All stages
Idents(neuropil) <- neuropil$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
tiff(filename = "Fig_3/medulla_neuropil_DimPlot_stages.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(neuropil, cols = colfunc(6), pt.size = 0.75) + NoLegend()
dev.off()


#For Neuropil X Glia
##All stages
Idents(neuropil_1) <- neuropil_1$stage
colfunc <- colorRampPalette(c("#FFE302", "magenta", "darkblue"))
tiff(filename = "Fig_3/neuropilX_DimPlot_stages.tiff",width = 7, height = 7,units = "in", res = 400)
DimPlot(neuropil_1, cols = colfunc(6), pt.size = 0.75) + NoLegend()
dev.off()



