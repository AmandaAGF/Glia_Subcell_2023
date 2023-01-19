glias_larva <- readRDS("glias_larva.rds")
larva = readRDS("larvaOL.integrated150.rds")
glia_adult <- readRDS("glia_adult.rds")
adult = readRDS("GSE142787_Adult.rds")


##### Fig. 1
#Save DotPlot for publication - new markers
a <- CellsByIdentities(glias_larva)
b <- levels(glias_larva)
larva$Fig1 <- "Non-Glial Cells"

for(i in 1:length(a)){
  names_cells <- a[[i]]
  index <- which(colnames(larva) %in% names_cells)
  larva@meta.data$Fig1[index] <- b[i]
}

Idents(larva) <- larva$Fig1
levels(larva) <- c("Non-Glial Cells", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")

tiff(filename = "Fig_1/glias_larva_DotPlot_FourMarkers.tiff",width = 8.5, height = 6,units = "in", res = 400)
DotPlot(larva, features = rev(c("repo", "CG32032", "AnxB9", "GstE12")), assay = "RNA", scale = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14), axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()
dev.off()


#Save DotPlot for publication - new markers - Fig. 1 - ADULT
#Load glia_adult final again
a <- CellsByIdentities(glia_adult)
b <- levels(glia_adult)
adult$Fig1 <- "Non-Glial Cells"

for(i in 1:length(a)){
  names_cells <- a[[i]]
  index <- which(colnames(adult) %in% names_cells)
  adult@meta.data$Fig1[index] <- b[i]
}

Idents(adult) <- adult$Fig1
levels(adult) <- c("Non-Glial Cells", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")

tiff(filename = "Fig_1/glia_adult_DotPlot_FourMarkers.tiff",width = 8.5, height = 6,units = "in", res = 400)
DotPlot(adult, features = rev(c("repo", "CG32032", "AnxB9", "GstE12")), assay = "RNA", scale = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14), axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()
dev.off()
