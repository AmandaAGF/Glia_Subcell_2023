glias_larva <- readRDS("glias_larva_2.rds")
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8
levels(glias_larva) <- c("20", "1", "3", "7", "12", "4", "19", "14", "8", "11", "15", "10", "16", "22", "9", "6", "21")
VlnPlot(glias_larva, features = c("nCount_RNA", "All_TFs"), stack = T, flip = T) + NoLegend()

#Change value of 12 cells expressing more than 6 in All_TFs to 6 for better visualization
glias_larva$All_TFs[which(glias_larva$All_TFs > 6)] <- 6

#Change value of 25 cells containing more than 50000 UMIs to 50000 for better visualization
glias_larva$nCount_RNA[which(glias_larva$nCount_RNA > 50000)] <- 50000

p = VlnPlot(glias_larva, features = c("nCount_RNA", "All_TFs"), stack = T, flip = T, fill.by = "ident",
            cols = c("#395FF8", "#395FF8", "#395FF8", "#395FF8", "#C39700", "#C39700", "#C57CFF", 
                     "#02BC67", "#02BC67", "#F35442", "#00BDC2", "#00BDC2","#64AF00", "#00A5FF", "#00A5FF", 
                     "#00A5FF", "#FC73D0")) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14), axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold"))
p$layers[[1]]$aes_params$size = 1
tiff(filename = "Fig_S6/glias_larva_VlnPlot_UMI_TFs.tiff",width = 10, height = 7,units = "in", res = 400)
p
dev.off()


#For adults
glia_adult <- readRDS("glia_adult_1.rds")

Idents(glia_adult) <- glia_adult$integrated_snn_res.0.6
tiff(filename = "Fig_S6/glia_adult_DimPlot.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glia_adult, label = T, pt.size = 0.1, label.size = 5) + NoLegend()
dev.off()

glia_adult <- FindClusters(glia_adult, resolution = 2)

tiff(filename = "Fig_S6/glia_adult_DimPlot_over.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glia_adult, label = T, pt.size = 0.1, label.size = 5) + NoLegend()
dev.off()



levels(glia_adult) <- c("26", "15", "21", "8", "30", "20", "14", "13", "6", "27", "16","24", "5", "18", "23", "1", "17", "25", "7", "11", "19", "0", "2", "22", "29", "9", "12", "3","4", "31", "28", "10")
VlnPlot(glia_adult, features = c("nCount_RNA", "All_TFs"), stack = T, flip = T) + NoLegend()

#Change value of 7 cells expressing more than 4 in All_TFs to 4 for better visualization
glia_adult$All_TFs[which(glia_adult$All_TFs > 4)] <- 4

#Change value of 36 cells containing more than 15000 UMIs to 15000 for better visualization
glia_adult$nCount_RNA[index <- which(glia_adult$nCount_RNA > 15000)] <- 15000

q = VlnPlot(glia_adult, features = c("nCount_RNA", "All_TFs"), stack = T, flip = T, fill.by = "ident",
            cols = c("#395FF8", "#395FF8", "#E5841B", "#C19500", "#FC73D0", "#00BDC2", "#B082FF", 
                     "#02BC67", "#02BC67", "#02BC67", "#64AF00", "#AB9F00", "#00A3FF", "#00A3FF", 
                     "#00A3FF", "#524B59", "#524B59", "#524B59", "#C1AED4", "#C1AED4", "#D7C5E8",
                     "#8E809C", "#8E809C", "#8E809C", "#00B8DC", "#00B8DC", "#A3C2AE", "#A3C2AE", 
                     "#677A6E", "#64AF00", "#C37CFF", "#C37CFF")) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14), axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold"))
q$layers[[1]]$aes_params$size = 1
tiff(filename = "Fig_S6/glia_adult_VlnPlot_UMI_TFs.tiff",width = 10, height = 7,units = "in", res = 400)
q
dev.off()




