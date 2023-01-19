##### Fig. 4
#Load the previous version of glias_larva and label All_TFs to be able to have Chiasm_PC and Chiasm_CB
glias_larva <- readRDS("glias_larva_2.rds")
Idents(glias_larva) <- glias_larva$integrated_snn_res.0.8

new.cluster.id_1 <- c("Cortex_PC", "Cortex_PC", "Chiasm_PC", "Medulla_Neuropil", "Cortex_PC", "Perineurial", "Medulla_Neuropil", "Epithelial_Marginal" ,"LobulaC_Neuropil", "Chiasm_CB", "Perineurial", "Epithelial_Marginal", "Satellite","Retina_Wrapping", "Cortex_CB", "Subperineurial", "Medulla_Neuropil")
names(new.cluster.id_1) <- levels(glias_larva)
glias_larva <- RenameIdents(glias_larva, new.cluster.id_1)
glias_larva@meta.data$new.cluster.id_1 <- glias_larva@active.ident
DimPlot(glias_larva, label = T) + NoLegend()

#Use this function to generate the diff expt between clusters VlnPlots
#For nCount_RNA of Chiasm
vp_case1 <- function(gene_sig, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(glias_larva, features = signature, assay = "RNA",
            pt.size = 0.1, 
            idents = c("Chiasm_CB", "Chiasm_PC"), #idents that I am going to be comparing
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c("#C49A00", "#00BDD2")
    ) + NoLegend() + stat_compare_means(comparisons = test_sign, label = "p.signif", show.legend = T, method = "wilcox.test") + stat_summary(fun = median, geom='point', size = 30, colour = "blue", shape = 95, show.legend = T) +
      theme(axis.text.x = element_text(vjust = 1, hjust=1, face = "bold", size = 12), axis.text=element_text(size=12),
            axis.title=element_text(size=12))
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_sig) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 6000) )
  }
  cowplot::plot_grid(plotlist = plot_list)
}

gene_sig <- c("nCount_RNA")
comparisons <- list(c("Chiasm_PC", "Chiasm_CB"))
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_nCount_RNA.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

#For All_TFs VlnPlot of Chiasm
#Change expression of the two cells expressing too much TFs to 5 for better visualization
Chiasm_PC <- subset(glias_larva, idents = "Chiasm_PC")
index <- which(Chiasm_PC$All_TFs > 5)
cells_to_change <- colnames(Chiasm_PC)[index]
index <- which(colnames(glias_larva) %in% cells_to_change)
glias_larva$All_TFs[index] <- 5

vp_case1 <- function(gene_sig, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(glias_larva, features = signature, assay = "RNA",
            pt.size = 0.1, 
            idents = c("Chiasm_CB", "Chiasm_PC"), #idents that I am going to be comparing
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c("#C49A00", "#00BDD2")
    ) + NoLegend() + stat_compare_means(comparisons = test_sign, label = "p.signif", show.legend = T, method = "wilcox.test") + stat_summary(fun = median, geom='point', size = 30, colour = "blue", shape = 95, show.legend = T) +
      theme(axis.text.x = element_text(vjust = 1, hjust=1, face = "bold", size = 12), axis.text=element_text(size=12),
            axis.title=element_text(size=12))
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_sig) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
}

gene_sig <- c("All_TFs")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_All_TFs.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

#For the other genes in Chiasm
##repo
gene_sig <- c("repo")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_repo.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

##Act5C
gene_sig <- c("Act5C")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_Act5C.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

##FK506-bp2
gene_sig <- c("FK506-bp2")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_FK506-bp2.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

##CG14598
gene_sig <- c("CG14598")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_CG14598.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()


#For nCount_RNA of Cortex
vp_case1 <- function(gene_sig, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(glias_larva, features = signature, assay = "RNA",
            pt.size = 0.1, 
            idents = c("Cortex_CB", "Cortex_PC"), #idents that I am going to be comparing
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c("#F8766D", "#FB61D7")
    ) + NoLegend() + stat_compare_means(comparisons = test_sign, label = "p.signif", show.legend = T, method = "wilcox.test") + stat_summary(fun = median, geom='point', size = 30, colour = "blue", shape = 95, show.legend = T) +
      theme(axis.text.x = element_text(vjust = 1, hjust=1, face = "bold", size = 12), axis.text=element_text(size=12),
            axis.title=element_text(size=12))
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_sig) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 6000) )
  }
  cowplot::plot_grid(plotlist = plot_list)
}

gene_sig <- c("nCount_RNA")
comparisons <- list(c("Cortex_CB", "Cortex_PC"))
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Cortex_nCount_RNA.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

#For All_TFs VlnPlot of Cortex
vp_case1 <- function(gene_sig, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(glias_larva, features = signature, assay = "RNA",
            pt.size = 0.1, 
            idents = c("Cortex_CB", "Cortex_PC"), #idents that I am going to be comparing
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c("#F8766D", "#FB61D7")
    ) + NoLegend() + stat_compare_means(comparisons = test_sign, label = "p.signif", show.legend = T, method = "wilcox.test") + stat_summary(fun = median, geom='point', size = 30, colour = "blue", shape = 95, show.legend = T) + 
      theme(axis.text.x = element_text(vjust = 1, hjust=1, face = "bold", size = 12), axis.text=element_text(size=12),
            axis.title=element_text(size=12))
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_sig) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
}

gene_sig <- c("All_TFs")
comparisons <- list(c("Cortex_CB", "Cortex_PC"))
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Cortex_All_TFs.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

#For the other genes in Cortex
##repo
gene_sig <- c("repo")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Cortex_repo.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

##fabp
gene_sig <- c("fabp")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Cortex_fabp.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()

##FK506-bp2
gene_sig <- c("FK506-bp2")
vp_case1(gene_sig, test_sign = comparisons)

#Save VlnPlot
tiff(filename = "Fig_4/glias_larva_VlnPlot_Cortex_FK506-bp2.tiff",width = 7, height = 9,units = "cm", res = 400)
vp_case1(gene_sig, test_sign = comparisons)
dev.off()




#Introns
#Load raw counts after cellranger considering introns
GG_PC_I <- Read10X("GG_Cell_Barcodes/bams/GG_PC_all/outs/raw_feature_bc_matrix_introns_new/")
GG_PC_I <- CreateSeuratObject(GG_PC_I)

#Load raw counts after cellranger NOT considering introns
GG_PC_E <- Read10X("GG_Cell_Barcodes/bams/GG_PC_all/outs/raw_feature_bc_matrix_nointrons_new/")
GG_PC_E <- CreateSeuratObject(GG_PC_E)

matrix_GG_PC_I <- as.matrix(GG_PC_I@assays$RNA@counts,"dgCMatrix")
matrix_GG_PC_E <- as.matrix(GG_PC_E@assays$RNA@counts,"dgCMatrix")

sum_counts_PC_I <- as.data.frame(colSums(matrix_GG_PC_I))
sum_counts_PC_E <- as.data.frame(colSums(matrix_GG_PC_E))
total_introns_PC <- ((sum_counts_PC_I-sum_counts_PC_E)/sum_counts_PC_I)*100


#Load raw counts after cellranger considering introns
GG_CB_I <- Read10X("GG_Cell_Barcodes/bams/GG_CB_all/outs/raw_feature_bc_matrix_introns_new/")
GG_CB_I <- CreateSeuratObject(GG_CB_I)

#Load raw counts after cellranger NOT considering introns
GG_CB_E <- Read10X("GG_Cell_Barcodes/bams/GG_CB_all/outs/raw_feature_bc_matrix_nointrons_new/")
GG_CB_E <- CreateSeuratObject(GG_CB_E)

matrix_GG_CB_I <- as.matrix(GG_CB_I@assays$RNA@counts,"dgCMatrix")
matrix_GG_CB_E <- as.matrix(GG_CB_E@assays$RNA@counts,"dgCMatrix")

sum_counts_CB_I <- as.data.frame(colSums(matrix_GG_CB_I))
sum_counts_CB_E <- as.data.frame(colSums(matrix_GG_CB_E))
total_introns_CB <- ((sum_counts_CB_I-sum_counts_CB_E)/sum_counts_CB_I)*100


plot <- as.data.frame(c(total_introns_PC$"colSums(matrix_GG_PC_I)", total_introns_CB$"colSums(matrix_GG_CB_I)"))
colnames(plot) <- "Introns"
plot$Identity <- "Chiasm_PC"
plot$Identity[326:length(plot$Introns)] <- "Chiasm_CB"

plot$Introns[which(plot$Introns > 10)] <- 10

tiff(filename = "Fig_4/glias_larva_VlnPlot_Chiasm_Introns.tiff",width = 7, height = 9,units = "cm", res = 400)
ggplot(plot, aes(x=Identity, y=Introns, fill=Identity)) + geom_violin(adjust =1,trim=TRUE, scale = "width") + 
  stat_compare_means(comparison = list(c("Chiasm_PC", "Chiasm_CB")), label = "p.signif", show.legend = T, method = "wilcox.test") + 
  stat_summary(fun = median, geom='point', size = 30, colour = "blue", shape = 95, show.legend = T) +
  scale_fill_manual(values=c("#00BDD2", "#C49A00")) + ylim(0,12.5) +
  scale_x_discrete(limits=c("Chiasm_PC", "Chiasm_CB")) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 0.3) + theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(vjust = 1, hjust=1, face = "bold", size = 12, angle = 45), axis.text=element_text(size=12),
        axis.title=element_text(size=12))
dev.off()





