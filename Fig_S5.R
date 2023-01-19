##### Fig. S5
glias_larva <- readRDS("glias_larva_2.rds")
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8
#DimPlot of the L3 glia
tiff(filename = "Fig_S5/glias_larva_DimPlot_res08.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(glias_larva, pt.size = 0.1, label = T, label.size = 5, cols = ) + NoLegend()
dev.off()

#Heatmap
DEG_glias_larva <- FindAllMarkers(glias_larva, only.pos = T, assay = "RNA")
top10_glias_larva <- DEG_glias_larva %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
a <- levels(glias_larva)
a <- a[c(-1, -2, -3, -5, -10, -15)]
a <- append(c("12", "4", "20", "7", "1", "3"), a)
levels(glias_larva) <- a
DoHeatmap(glias_larva, features = top10_glias_larva$gene) + NoLegend()
#Error: The following features were omitted as they were not found in the scale.data slot for the integrated assay: l(3)neo38, CG43778, Mdr65, SoxN, crim, Tsp42Ea, alphaTub85E, Tet, CG2611, Hsp27
#Remove those genes
top10_glias_larva <-  top10_glias_larva[-which(top10_glias_larva$gene == "l(3)neo38" | top10_glias_larva$gene == "CG43778" | top10_glias_larva$gene == "Mdr65"|
        top10_glias_larva$gene == "SoxN" | top10_glias_larva$gene == "crim" | top10_glias_larva$gene == "Tsp42Ea" |
        top10_glias_larva$gene == "alphaTub85E" | top10_glias_larva$gene == "Tet" | top10_glias_larva$gene == "CG2611" |
        top10_glias_larva$gene == "Hsp27"),]
tiff(filename = "Fig_S5/glias_larva_Heatmap.tiff",width = 20, height = 30,units = "cm", res = 400)
DoHeatmap(glias_larva, features = top10_glias_larva$gene, group.colors = c("#00BFC4", "#CD9600", "#ED68ED", "#7CAE00", "#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", 
                                                                           "#00C19A", "#00B8E7", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1")) + NoLegend()
dev.off()
