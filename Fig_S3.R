####Fig. S3

glia_adult <- readRDS("glia_adult_1.rds")


#DotPlot for publication - Fig. S2
levels(glia_adult) <- c("Perineurial", "Fenestrated", "Subperineurial", "Chalice", "Pseudocartridge", "Distal_Satellite", "Proximal_Satellite", 
                        "Epithelial", "Marginal", "Outer_Chiasm", "Inner_Chiasm", "Cortex_1", "Cortex_2", "Medulla_EG_1", "Medulla_EG_2", 
                        "Medulla_ALG", "Astrocytes_3", "Ensheathing_1", "Ensheathing_2","Astrocytes_2", "Astrocytes_1", "Astrocytes_4")
tiff(filename = "Fig_S3/glia_adult_DotPlot_validations.tiff",width = 10, height = 6,units = "in", res = 400)
DotPlot(glia_adult, features = rev(c("Tret1-1", "Vmat","rost", "CG3036", "vvl", "Optix","b","GstT4" ,"DAT", "CG31663", "ome", "nkd")), assay = "RNA") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold", size = 14),
                                                                                                                                                             axis.text.y = element_text(vjust = 0.5, hjust=1, face = "bold", size = 18),
                                                                                                                                                             axis.title=element_text(face="bold"))
dev.off()
