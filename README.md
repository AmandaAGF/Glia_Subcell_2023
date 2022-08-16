# Glia_Subcell_2022

This repo contains most of the information needed to generate the Figures of the paper "A Novel Approach to Identify Subcellularly Localized Transcripts in Single Cells".

It contains:

R files for the preparation and annotation of the Adult and Larva datasets:
1_Glia_Adult.R 2_Glia_Larva.R

R file for the separation of Epithelial and Marginal glia on the adult dataset:
XXXXXXXXX

R file for finding which cluster in larva are glia based on correlation with repo- and elav-positive FACSorted cells:
2b_Glia_Larva_Correlations.R

R file





R files that contain the commands used to generate the different Figures. 3_Figure1.R 4_Figure2_S3.R 5_Figure4_5_S7.R 6_Figure6_S10.R 7_FigureS1.R 8_FigureS2.R 9_FigureS11.R

Supporting functions that are used in the above scripts: classifier_utils.R functions.R LayeredFeaturePlot_white_background_pt.brightness.R

Tables that are used in the different scripts: ConversionTable.csv : Table from Flybase that converts Drosophila gene names to FGgn. FlyBase_TFs.csv : All Drosophila TFs, as indicated in Flybase mart_export_human.txt : A table that contains the GO terms for all human genes mart_export.txt : A table that contains the GO terms for all Drosophila genes colors.txt : A vector that color codes the pupal clusters according to their predicted temporal window colors2.txt : A vector that color codes non-medulla clusters used in EDF 9.

Datasets that are used in this project that are available at GEO: GSE118953_raw_count.tsv : mouse cortical scSeq (Telley et al, Science 2019) P15_10Xv3Adjusted.rds : scSeq data from Drosophila optic lobes at P15 stage P30_10Xv3Adjusted.rds : scSeq data from Drosophila optic lobes at P30 stage P40_10Xv3Adjusted.rds : scSeq data from Drosophila optic lobes at P40 stage P50_10Xv3Adjusted.rds : scSeq data from Drosophila optic lobes at P50 stage P70_10Xv3Adjusted.rds : scSeq data from Drosophila optic lobes at P70 stage
