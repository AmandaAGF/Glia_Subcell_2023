library(Seurat)
library(ggplot2)

#Run script 1b_Glia_Larva_Correlations first
#ggplot for Glia vs Non-Glial Cells
tiff(filename = "Fig_S1/larva_ggplot_Glia_vs_NonGlia.tiff",width = 25, height = 15,units = "cm", res = 400)
ggplot(data = cor_FACS_Bulk_scRNAseq, mapping = aes(x = Elav, y = Repo, color = identity, label = label)) + 
  geom_point() +  scale_color_manual(values = c("red", "black")) +  
  geom_text(data = cor_FACS_Bulk_scRNAseq[which(cor_FACS_Bulk_scRNAseq$identity == "Glial Cluster"),],aes(label=label),hjust=0, vjust=0, size = 6) +
  theme(panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size = 20, face = "bold"), plot.title = element_text(size = 20), axis.text.y = element_text(size = 15))
dev.off()


larva = readRDS("larvaOL.integrated150.rds")

#Separate Glia vs No-Glia in Larva
index <- WhichCells(larva, idents = c("44", "87", "45", "56", "67", "91", "88", "63", "85", "2", "3"))
index_1 <- which(colnames(larva) %in% index)
larva$Glia <- "Non-Glial Cells"
larva$Glia[index_1] <- "Glia"
DimPlot(larva, group.by = "Glia")
Idents(larva) <- larva$Glia



#DimPlot with GLia vs Non-Glial Cells
Idents(larva) <- larva$Glia
tiff(filename = "Fig_S1/larva_DimPlot_Glia_vs_NonGlia.tiff",width = 15, height = 15,units = "cm", res = 400)
DimPlot(larva, cols = c("#4f4e4e","#f50505")) + NoLegend()+ theme(axis.title=element_text(size=14,face="bold"))
dev.off()

#Save DotPlot for publication
#Get metadata$Fig1 from script Fig_1
Idents(larva) <- larva$Fig1
tiff(filename = "Fig_S1/larva_DotPlot_AllMarkers.tiff",width = 8, height = 6,units = "in", res = 400)
DotPlot(larva, features = rev(c("repo", "CG9686", "Ama", "Npc2a", "MtnA", "CG32032", "AnxB9", "CG7433", "CG15209", "GstE12")), assay = "RNA") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14), axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()
dev.off()

###### P96 Optic Lobe Dataset - Kurmangaliyev 2020 GSE156455
#Read Adult data
zipursky_96h = readRDS("Kurmangaliyev96h_classified_from_OzelAdult_Seurat.RDS")

#Study Z data set to find glial clusters - Only clusters that were not annotated as a specific neuron
Idents(zipursky_96h) = zipursky_96h@meta.data$type
DotPlot(zipursky_96h, features = c("repo", "CG32032", "AnxB9", "GstE12"), assay = "RNA" , stack = T, flip = T,idents = c("C112",   "C130",   "C16" ,   "C165",   "C167",   "C192",    "C32" ,   "C46" ,  "C67" ,   "C68" ,   "C72" ,   "C86" ,     "G113",   "G121" ,  "G125"  , "G145" ,  "G157" ,  "G173"   ,"G186",  "G193"  , "G37"   , "G43" ,   "G49"  ,  "G51"  ,  "G57"  ,  "G63"   , "G66"  ,  "G74"  ,  "G77" ,  "G78"  ,  "G91" ,      "N1"   ,  "N100"  , "N101" ,  "N102"  , "N103" ,  "N104"   ,"N106",  "N108" ,  "N110"  , "N111" ,  "N114" ,  "N115" ,  "N116" ,  "N117" ,  "N118" ,  "N119"  , "N120",  "N122" ,  "N123" ,  "N124"  , "N126"  , "N128" ,  "N129",   "N131" ,  "N132"  , "N134" ,  "N135" , "N136" ,  "N137" ,  "N138" ,  "N139"  , "N140" ,  "N142" ,  "N143"  , "N146" ,  "N148",   "N150" , "N151" ,  "N152"  , "N153"  , "N154" ,  "N155",   "N156" ,  "N158"  , "N159" ,  "N160" ,  "N161" , "N162"  , "N163" ,  "N164" ,  "N168"  , "N169"  , "N170",   "N172" ,  "N174"  , "N176"  , "N177"  ,"N178"  , "N179",   "N180" ,  "N181" ,  "N182" ,  "N183"  , "N188" ,  "N19"   , "N190"  , "N194" , "N195"  , "N33"  , "N34" ,   "N35",    "N39"   , "N42"  ,  "N48" ,   "N52" ,   "N54" ,   "N55" ,  "N56" ,   "N58" ,   "N59"  ,  "N61"   , "N65" ,   "N69"  ,  "N70" ,   "N71"  ,  "N73"   , "N75" ,  "N76"  ,  "N79"  ,  "N80"  ,  "N81"    ,"N82"  ,  "N83"  ,  "N84" ,   "N85" ,   "N87"  ,  "N88",   "N89"  ,  "N90"  ,  "N92"  ,  "N93"    ,"N95" ,   "N96"  ,  "N99"))+NoLegend()

#Save graphs for publication - Fig. S1
tiff(filename = "Fig_S1/Zipursky_DotPlot_FourMarkers.tiff",width = 25, height = 5,units = "in", res = 400)
DotPlot(zipursky_96h, features = rev(c("repo", "CG32032", "AnxB9", "GstE12")), assay = "RNA",
        idents = c("C112",   "C130",   "C16" ,   "C165",   "C167",   "C192",    "C32" ,   "C46" ,  
                   "C67" ,   "C68" ,   "C72" ,   "C86" ,     "G113",   "G121" ,  "G125"  , "G145" ,  
                   "G157" ,  "G173"   ,"G186",  "G193"  , "G37"   , "G43" ,   "G49"  ,  "G51"  ,  
                   "G57"  ,  "G63"   , "G66"  ,  "G74"  ,  "G77" ,  "G78"  ,  "G91" ,      "N1"   ,  
                   "N100"  , "N101" ,  "N102"  , "N103" ,  "N104"   ,"N106",  "N108" ,  "N110"  , 
                   "N111" ,  "N114" ,  "N115" ,  "N116" ,  "N117" ,  "N118" ,  "N119"  , "N120",  
                   "N122" ,  "N123" ,  "N124"  , "N126"  , "N128" ,  "N129",   "N131" ,  "N132"  , 
                   "N134" ,  "N135" , "N136" ,  "N137" ,  "N138" ,  "N139"  , "N140" ,  "N142" ,  
                   "N143"  , "N146" ,  "N148",   "N150" , "N151" ,  "N152"  , "N153"  , "N154" ,  
                   "N155",   "N156" ,  "N158"  , "N159" ,  "N160" ,  "N161" , "N162"  , "N163" ,  
                   "N164" ,  "N168"  , "N169"  , "N170",   "N172" ,  "N174"  , "N176"  , "N177"  ,
                   "N178"  , "N179",   "N180" ,  "N181" ,  "N182" ,  "N183"  , "N188" ,  "N19"   , 
                   "N190"  , "N194" , "N195"  , "N33"  , "N34" ,   "N35",    "N39"   , "N42"  ,  
                   "N48" ,   "N52" ,   "N54" ,   "N55" ,  "N56" ,   "N58" ,   "N59"  ,  "N61"   , 
                   "N65" ,   "N69"  ,  "N70" ,   "N71"  ,  "N73"   , "N75" ,  "N76"  ,  "N79"  ,  
                   "N80"  ,  "N81"    ,"N82"  ,  "N83"  ,  "N84" ,   "N85" ,   "N87"  ,  "N88",   
                   "N89"  ,  "N90"  ,  "N92"  ,  "N93"    ,"N95" ,   "N96"  ,  "N99")) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 14,vjust = 0.5), 
        axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()
dev.off()



###### VNC dataset
library(loomR)
VNC_loom <- connect(filename = "VNC/Goodwin_Fly_AdultVNC_elife54074.loom", mode = "r+")

#Extract cell_ID from the loom object
cell_ID <- VNC_loom$col.attrs$CellID
cell_ID <- cell_ID[1:cell_ID$dims]

#Extract gene_ID from the loom object
gene_ID <- VNC_loom$row.attrs$Gene
gene_ID <- gene_ID[1:gene_ID$dims]

#Extract UMI counts from the loom object
VNC_matrix <- VNC_loom[["matrix"]] #I did before with VNC_loom$matrix, but this way stopped working and [[]] works
VNC_matrix <- t(VNC_matrix[1:length(cell_ID), 1:length(gene_ID)])

#Name columns with cell_ID and rows with gene_ID
rownames(VNC_matrix) <- gene_ID
colnames(VNC_matrix) <- cell_ID

#Make Seurat object
VNC <- CreateSeuratObject(VNC_matrix)
VNC <- NormalizeData(VNC)

#Add Clusters in the metadata
VNC_loom$col.attrs$ClusterID
attrs <- c("ClusterName", "ClusterID")
attr.df <- VNC_loom$get.attribute.df(MARGIN = 2, attributes = attrs)
head(attr.df)
VNC$ClusterName <- attr.df[,1]
VNC$ClusterID <- attr.df[,2]

Idents(VNC) <- VNC$ClusterID
a <- levels(VNC)
b <- c("98", "24", "70", "106", "23", "99", "116")
c <- setdiff(a, b)
c <- append(b,c)
levels(VNC) <- c
VlnPlot(VNC, features = c("repo", "CG32032", "AnxB9", "GstE12", "Gabat"), assay = "RNA" , stack = T, flip = T) + NoLegend()

#saveRDS(VNC, "VNC.rds")
VNC <- readRDS("VNC.rds")

#Save graphs for publication - Fig. S1
tiff(filename = "Fig_S1/VNC_DotPlot_FourMarkers.tiff",width = 25, height = 5,units = "in", res = 400)
DotPlot(VNC, features = rev(c("repo", "CG32032", "AnxB9", "GstE12")), assay = "RNA") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 14,vjust = 0.5), 
        axis.text=element_text(size=12), axis.text.y = element_text(face = "bold", size = 14),                                                                                                                                                                                                                                                
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()
dev.off()

