glias_larva <- readRDS("glias_larva_1.rds")


##### Fig. S2
#DotPlot for publication - Fig. S2
levels(glias_larva) <- c("Subperineurial", "Perineurial", "Satellite", "Epithelial_Marginal", "Outer_Chiasm", "Inner_Chiasm", "Chiasm", "Medulla_Neuropil", "LobulaC_Neuropil", "Cortex", "Retina_Wrapping", "Cortex_Central_Brain", "Central_Brain")
tiff(filename = "Fig_S2/glias_larva_DotPlot.tiff",width = 10, height = 6,units = "in", res = 400)
DotPlot(glias_larva, features = rev(c("rost", "svp", "Optix", "Dll", "CG14598", "DAT", "hth", 
                                      "Hand", "Cpr51A","CG4288", "Cyp4g15","alphaTub85E", "Npc2b", 
                                      "Adgf-D")), assay = "RNA") + coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold", 
                                   size = 14), axis.text.y = element_text(vjust = 0.5, 
                                                                          hjust=1, face = "bold", size = 18), 
        axis.title=element_text(face="bold"))
dev.off()
