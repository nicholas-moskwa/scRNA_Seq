      #Version 4
    #This version provides all data for 2022 publication
	#These scripts for working with in-vivo E16 SMG dataset
  #Major parameter and statistical decisions made:
#min.cells = 3, min.features = 200
#subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 7
#dims = 1:40
#resolution = 2.0
#Extracts UMAP data, gene differentiation lists between all clusters
#and creates graphs with standard frame sizing
#Labels both broad and specific to add to graphs
  #Following post-processing additions
#1: code to subset out stromal data without re-calculation
#2: code to subset out and re-calculate epithelial and stromal data individually. This data should be used for publications.
#Epithelium: dims = 1:21 
#Epithelium: resolution = 1.1
#Stroma: dims = 1:40 
#Stroma: resolution = 0.8
#3: code to regress out Cell cycling transcriptome from stromal data. This data is not used for publication.
#4: code to integrate datasets from different scRNA seq experiments.  This data is not used for publication.


############################## sessionInfo() ###################################################
################################################################################################

#R version 3.6.3 (2020-02-29)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Debian GNU/Linux 10 (buster)

#Matrix products: default
#BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] stringr_1.4.0   magrittr_1.5    patchwork_1.0.0 pheatmap_1.0.12
#[5] ggplot2_3.3.0   dplyr_0.8.5     cowplot_1.0.0   Seurat_3.1.5   

#loaded via a namespace (and not attached):
#  [1] httr_1.4.1          tidyr_1.0.2         jsonlite_1.6.1     
#[4] viridisLite_0.3.0   splines_3.6.3       lsei_1.2-0         
#[7] leiden_0.3.3        gtools_3.8.2        assertthat_0.2.1   
#[10] ggrepel_0.8.2       globals_0.12.5      pillar_1.4.3       
#[13] lattice_0.20-44     glue_1.4.0          reticulate_1.15    
#[16] digest_0.6.25       RColorBrewer_1.1-2  colorspace_1.4-1   
#[19] htmltools_0.4.0     Matrix_1.2-18       plyr_1.8.6         
#[22] pkgconfig_2.0.3     tsne_0.1-3          listenv_0.8.0      
#[25] purrr_0.3.4         scales_1.1.0        RANN_2.6.1         
#[28] gdata_2.18.0        RSpectra_0.16-0     Rtsne_0.15         
#[31] tibble_3.0.1        farver_2.0.3        ellipsis_0.3.0     
#[34] withr_2.4.2         ROCR_1.0-7          pbapply_1.4-2      
#[37] lazyeval_0.2.2      survival_3.2-7      crayon_1.4.1       
#[40] future_1.17.0       nlme_3.1-152        MASS_7.3-51.5      
#[43] gplots_3.0.3        ica_1.0-2           tools_3.6.3        
#[46] fitdistrplus_1.0-14 data.table_1.12.8   lifecycle_0.2.0    
#[49] plotly_4.9.2.1      munsell_0.5.0       cluster_2.1.0      
#[52] irlba_2.3.3         compiler_3.6.3      rsvd_1.0.3         
#[55] caTools_1.18.0      rlang_0.4.5         grid_3.6.3         
#[58] ggridges_0.5.2      rstudioapi_0.13     RcppAnnoy_0.0.16   
#[61] rappdirs_0.3.1      htmlwidgets_1.5.1   igraph_1.2.5       
#[64] labeling_0.3        bitops_1.0-6        npsurv_0.4-0       
#[67] gtable_0.3.0        codetools_0.2-16    reshape2_1.4.4     
#[70] R6_2.5.0            gridExtra_2.3       zoo_1.8-7          
#[73] uwot_0.1.8          future.apply_1.5.0  KernSmooth_2.23-16 
#[76] ape_5.4             stringi_1.4.6       parallel_3.6.3     
#[79] Rcpp_1.0.4.6        vctrs_0.2.4         sctransform_0.2.1  
#[82] png_0.1-7           tidyselect_1.0.0    lmtest_0.9-37  

############################## Libraries and start up ##########################################
################################################################################################

#Allows fro python UMAP calculations
reticulate::use_condaenv(condaenv = "seurat", conda = "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/miniconda3/bin/conda")

library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(magrittr)
library(stringr)
library(biomaRt)
#Sets the working directory where the single cell data
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04")

#Reads the 10x Generated data that is in matrix form; barcodes, features, matrix, into a large data matrix
agg.E16.data <- Read10X(data.dir =
                          "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/68.195.245.254/results/filtered_feature_bc_matrix")
#Creates Seurat object using the 10x raw data matrix
agg.E16 <- CreateSeuratObject(agg.E16.data, project = "agg.E16", min.cells = 3,
                              min.features = 200)
#Creates new parameter selecting for mitochondrial genes
agg.E16[["percent.mt"]] <- PercentageFeatureSet(agg.E16, pattern = "^mt-")

#Generates Volcano Plots to show parameters of dataset
VlnPlot(agg.E16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("Volcano_plot_Feature_RNAcount_Mitochondira_E16_04.pdf")
VlnPlot(agg.E16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Saves a scatterplot of Features vs mitochondrial genes
plot1 <- FeatureScatter(agg.E16, feature1 = "nCount_RNA", feature2 = "percent.mt")
pdf("Feature_vs_Mt_gene_ScatternPlot_E16_04.pdf")
plot1
dev.off()

#Saves a scatterplot of Features vs Counts
plot2 <- FeatureScatter(agg.E16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
pdf("Feature_vs_RNAcount_ScatterPlot_E16_04.pdf")
plot2
dev.off()

#Displays both plots to determine filtering of next step
plot1 + plot2

#Filter out cells that have unique feature counts more than 200 and less than 9000, and filter out cells that have more than 7% mitochondrial counts because those cells are likely dying
#7% was choosen because if your view the mitochondrial distribution as a normal distribution the lower bound would be 0, the center would be ~3.5%, and the upperbound would be ~7%
agg.E16 <- subset(agg.E16, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 7)

VlnPlot(agg.E16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("Volcano_plot_Feature_RNAcount_Mitochondira_E16_04_AfterMtRemoval.pdf")
VlnPlot(agg.E16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Normalizes dataset
agg.E16 <- NormalizeData(agg.E16)

#Identifies the most highly variable features, and returns 2000 features per dataset as a default
agg.E16 <- FindVariableFeatures(agg.E16, selection.method = "vst", nfeatures = 2000)
#generates list of Top 10 most variable genes
top10 <- head(VariableFeatures(agg.E16), 10)
#Saves plot of selected features
plot1 <- VariableFeaturePlot(agg.E16)
pdf("Varianceplot_E16_04.pdf")
plot1
dev.off()

#Saves plot of selected features and displays top10 most variable genes
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("Varianceplot_top10_E16_04.pdf")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()

#Labels all genes
all.genes <- rownames(agg.E16)
#Applies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.E16 <- ScaleData(agg.E16, feature = all.genes) 
#Runs PCA statistics used for generates a PCA plot of the first PCs
agg.E16 <- RunPCA(agg.E16, features = VariableFeatures(object = agg.E16))

#Generates a heatmap showing which genes are contributing to the hetergeneity that determines PC1, dims = "PC number", shows which genes are contributing most to the cell seperation.
pdf("PCA1_Heatmap_E16_04.pdf")
DimHeatmap(agg.E16, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("PCA2_Heatmap_E16_04.pdf")
DimHeatmap(agg.E16, dims = 2, cells = 500, balanced = TRUE)
dev.off()
pdf("PCA3_Heatmap_E16_04.pdf")
DimHeatmap(agg.E16, dims = 3, cells = 500, balanced = TRUE)
dev.off()

#Saves PCA plot total nomalized dataset based on PC 1 and PC 2
pdf("PCA_Single_cells_E16_04.pdf")
DimPlot(agg.E16, reduction = "pca")
dev.off()

#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
#Around 40 dimensions is where adding another dimentions adds little more extra separation by the Principle Components
pdf("Elbowplot_E16_04.pdf")
ElbowPlot(agg.E16, ndims = 50)
dev.off()

#constructs a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset (first 40 PCs, dims = 40), how many dimentions do you want included for the analysis.
agg.E16 <- FindNeighbors(agg.E16, dims = 1:40) 
#Generates the Clusters, resolution increases the specificity of each cluster, resolution = 2.0, the default is 1.0
agg.E16 <- FindClusters(agg.E16, resolution = 2.0)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.E16 <- RunTSNE(agg.E16, dims = 1:40)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.E16 <- RunUMAP(agg.E16, dims = 1:40)

#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.E16)), file = "Cell_per_Cluster_E16_04.csv")

#Generates a t-SNE plot
pdf("TSNE_E16_04.pdf")
  DimPlot(agg.E16, reduction = "tsne", label = TRUE)
dev.off() 

#Generates a UMAP plot
pdf("UMAP_E16_04.pdf", width = 8, height = 6)
  DimPlot(agg.E16, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
    NoLegend() +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15))
dev.off()

#Saves Formated dataset after all major computational processing
saveRDS(agg.E16, file = "E16_SMG_Raw_(SEURAT_v4)_04.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04")
agg.E16 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/E16_SMG_Raw_(SEURAT_v4)_04.rds")


############################## Gene for cell identification ####################################
################################################################################################

    #All UMAP/Feature gene plot maps used for assigning clusters labels
	#NOTE at the end of these there is a coded For-loop that will do this automatically for any gene lists.
  #LYMPHATICS ENDOTHELIUM
#Shows individual cells on cluster plot, Leucocytes = PTPRC
pdf("LymphaticEC_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Ptprc", "Pdpn", "Prox1", "Lyve1"), min.cutoff = "q9")
dev.off()
pdf("LymphaticEC_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Ptprc", "Pdpn", "Prox1", "Lyve1"))
dev.off()

  #RED BLOOD CELLS
#Shows individual cells on cluster plot, Erythrocytes = Hemoglobins
pdf("Erythrocytes_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Hba-a1", "Hba-a2", "Hbb-bs"), min.cutoff = "q9")
dev.off()
pdf("Erythrocytes_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Hba-a1", "Hba-a2", "Hbb-bs"))
dev.off()

  #ENDOTHELIUM
#Shows individual cells on cluster plot, Endothelial = Kdr and Cdh5
pdf("Endothelial_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Kdr", "Cdh5", "Pecam1", "Tie1", "Vcam1"), min.cutoff = "q9")
dev.off()
pdf("Endothelial_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Kdr", "Cdh5", "Pecam1", "Tie1", "Vcam1"))
dev.off()

  #NEURONS
#Shows individual cells on cluster plot, Neurons = Tubb3
pdf("Neuron_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Tubb3"), min.cutoff = "q9")
dev.off()
pdf("Neuron_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Tubb3"))
dev.off()

  #SCHWANN cells / Oligodendrocytes
pdf("Schwann_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Fabp7", "Mbp", "Plp1", "Ptprz1", "Sox10"), min.cutoff = "q9")
dev.off()
pdf("Schwann_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Fabp7", "Mbp", "Plp1", "Ptprz1", "Sox10"))
dev.off()

  #STROMA
#Shows individual cells on cluster plot, Stroma = PDGFRa
pdf("PDGFRa_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Pdgfra"), min.cutoff = "q9")
dev.off()
pdf("PDGFRa_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Pdgfra"))
dev.off()
#Shows individual cells on cluster plot, Stroma = PDGFRb
pdf("PDGFRb_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Pdgfrb"), min.cutoff = "q9")
dev.off()
pdf("PDGFRb_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Pdgfrb"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Endoglin
pdf("Eng_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Eng"), min.cutoff = "q9")
dev.off()
pdf("Eng_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Eng"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Vimentin
pdf("Vim_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Vim"), min.cutoff = "q9")
dev.off()
pdf("Vim_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Vim"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Thy1
pdf("Thy1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Thy1"), min.cutoff = "q9")
dev.off()
pdf("Thy1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Thy1"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Sca-1
pdf("Sca-1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Ly6a"), min.cutoff = "q9")
dev.off()
pdf("Sca-1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Ly6a"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Nt5e
pdf("Nt5e_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Nt5e"), min.cutoff = "q9")
dev.off()
pdf("Nt5e_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Nt5e"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Col1a1
pdf("Col1a1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Col1a1"), min.cutoff = "q9")
dev.off()
pdf("Col1a1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Col1a1"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Gli1
pdf("Gli1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Gli1"), min.cutoff = "q9")
dev.off()
pdf("Gli1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Gli1"))
dev.off()
#Shows individual cells on cluster plot, Stroma = FGF2
pdf("FGF2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Fgf2"), min.cutoff = "q9")
dev.off()
pdf("FGF2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Fgf2"))
dev.off()
#Shows individual cells on cluster plot, Stroma = FGF10
pdf("FGF10_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Fgf10"), min.cutoff = "q9")
dev.off()
pdf("FGF10_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Fgf10"))
dev.off()
#Shows individual cells on cluster plot, Stroma = BMP4
pdf("BMP4_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Bmp4"), min.cutoff = "q9")
dev.off()
pdf("BMP4_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Bmp4"))
dev.off()
#Shows individual cells on cluster plot, Stroma = BMP7
pdf("BMP7_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Bmp7"), min.cutoff = "q9")
dev.off()
pdf("BMP7_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Bmp7"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Kit Ligand
pdf("KitL_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Kitl"), min.cutoff = "q9")
dev.off()
pdf("Kitl_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Kitl"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Tgfb1
pdf("Tgfb1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Tgfb1"), min.cutoff = "q9")
dev.off()
pdf("Tgfb1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Tgfb1"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Tgfb2
pdf("Tgfb2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Tgfb2"), min.cutoff = "q9")
dev.off()
pdf("Tgfb2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Tgfb2"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Tgfbr3
pdf("Tgfbr3_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Tgfbr3"), min.cutoff = "q9")
dev.off()
pdf("Tgfbr3_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Tgfbr3"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Wnt2
pdf("Wnt2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Wnt2"), min.cutoff = "q9")
dev.off()
pdf("Wnt2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Wnt2"))
dev.off()
#Shows individual cells on cluster plot, Stroma = Wnt4
pdf("Wnt4_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Wnt4"), min.cutoff = "q9")
dev.off()
pdf("Wnt4_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Wnt4"))
dev.off()

#Shows individual cells on cluster plot, Col4a1
pdf("Col4a1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Col4a1"), min.cutoff = "q9")
dev.off()
pdf("Col4a1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Col4a1"))
dev.off()

  #EPITHELIUM
#Shows individual cells on cluster plot, Acinar Epi = Aqp5
pdf("Aqp5_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Aqp5"), min.cutoff = "q9")
dev.off()
pdf("Aqp5_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Aqp5"))
dev.off()
#Shows individual cells on cluster plot, Proacinar Epi = Bhlha15
pdf("Bhlha15_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Bhlha15"), min.cutoff = "q9")
dev.off()
pdf("Bhlha15_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Bhlha15"))
dev.off()
#Shows individual cells on cluster plot, Proacinar Epi = Bpifa2
pdf("PSP_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Bpifa2"), min.cutoff = "q9")
dev.off()
pdf("PSP_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Bpifa2"))
dev.off()
#Shows individual cells on cluster plot, Proacinar Epi = Muc19
pdf("Muc19_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Muc19"), min.cutoff = "q9")
dev.off()
pdf("Muc19_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Muc19"))
dev.off()
#Shows individual cells on cluster plot, Proacinar Epi = Smgc
pdf("Smgc_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Smgc"), min.cutoff = "q9")
dev.off()
pdf("Smgc_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Smgc"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Cdh1
pdf("Ecad_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Cdh1"), min.cutoff = "q9")
dev.off()
pdf("Ecad_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Cdh1"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Epcam
pdf("Epcam_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Epcam"), min.cutoff = "q9")
dev.off()
pdf("Epcam_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Epcam"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Cytokeratin 5
pdf("Krt5_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Krt5"), min.cutoff = "q9")
dev.off()
pdf("Krt5_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Krt5"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Cytokeratin 14
pdf("Krt14_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Krt14"), min.cutoff = "q9")
dev.off()
pdf("Krt14_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Krt14"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Cytokeratin 19
pdf("Krt19_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Krt19"), min.cutoff = "q9")
dev.off()
pdf("Krt19_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Krt19"))
dev.off()
#Shows individual cells on cluster plot, Epithelial = Smooth muscle actin
pdf("Acta2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Acta2"), min.cutoff = "q9")
dev.off()
pdf("Acta2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Acta2")) + NoLegend()
dev.off()
#Shows individual cells on cluster plot, Epithelial = Calponin 1
pdf("Cnn1_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Cnn1"), min.cutoff = "q9")
dev.off()
pdf("Cnn1_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Cnn1"))
dev.off()
#Shows individual cells on cluster plot, Epithelium = Sox10
pdf("Sox10_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Sox10"), min.cutoff = "q9")
dev.off()
pdf("Sox10_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Sox10"))
dev.off()
#Shows individual cells on cluster plot, Epithelium = Sox2
pdf("Sox2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Sox2"), min.cutoff = "q9")
dev.off()
pdf("Sox2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Sox2"))
dev.off()

  #IMMUNE
#Shows individual cells on cluster plot, Immune Memory CD4+ = Il7r
pdf("Il7r_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Il7r"), min.cutoff = "q9")
dev.off()
pdf("Il7r_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Il7r"))
dev.off()
#Shows individual cells on cluster plot, Immune Memory CD4+ = S100A4
pdf("S100a4_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("S100a4"), min.cutoff = "q9")
dev.off()
pdf("S100a4_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("S100a4"))
dev.off()
#Shows individual cells on cluster plot, Immune Mast cell = Ms4a2
pdf("Ms4a2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Ms4a2"), min.cutoff = "q9")
dev.off()
pdf("Ms4a2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Ms4a2"))
dev.off()
#Shows individual cells on cluster plot, Immune Mast cell = Hdc
pdf("Hdc_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Hdc"), min.cutoff = "q9")
dev.off()
pdf("Hdc_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Hdc"))
dev.off()

  #CELL CYCLE / Mitosis
#Shows individual cells on cluster plot, Mitosis G2 = "Cenpa", "Cenpe", "Cenpf", "Mki67"
pdf("Mitosis_G2_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Cenpa", "Cenpe", "Cenpf", "Mki67"), min.cutoff = "q9")
dev.off()
pdf("Mitosis_G2_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Cenpa", "Cenpe", "Cenpf", "Mki67"))
dev.off()
#Shows individual cells on cluster plot, Mitosis G1/S = "Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae"
pdf("Mitosis_G1_S_Clusters_E16_04.pdf")
FeaturePlot(agg.E16, features = c("Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae"), min.cutoff = "q9")
dev.off()
pdf("Mitosis_G1_S_ViolinPlot_E16_04.pdf")
VlnPlot(agg.E16, features = c("Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae"))
dev.off()

  #loop to make UMAP/Feature gene graphs and violin plots for assigning cluster labels
#List of genes to run in loop
Genes_of_Interest <- c("Ptprc", "Pdpn", "Prox1", "Lyve1","Hba-a1", "Hba-a2",
                      "Hbb-bs", "Kdr", "Cdh5", "Pecam1", "Tie1", "Vcam1", "Tubb3",
                      "Fabp7", "Mbp", "Plp1", "Ptprz1", "Pdgfra", "Pdgfrb",
                      "Eng", "Vim", "Thy1", "Ly6a", "Nt5e", "Col1a1", "Gli1", "Fgf2",
                      "Fgf10", "Bmp4", "Bmp7", "Kitl", "Tgfb1", "Tgfb2", "Tgfbr3",
                      "Wnt2", "Wnt4", "Col4a1", "Aqp5", "Bhlha15", "Bpifa2", "Muc19",
                      "Smgc", "Cdh1", "Epcam", "Krt5", "Krt14", "Krt19", "Acta2",
                      "Cnn1", "Sox10", "Sox2", "Il7r", "S100a4", "Ms4a2", "Hdc",
                      "Cenpa", "Cenpe", "Cenpf", "Mki67", "Hist1h1a", "Hist1h1b",
                      "Top2a", "Hist1h2ae", "Cd74", "Lsp1", "H2-Aa", "H2-Ab1")
Stromal_Markers <- c("Thy1", "Acta2","Pdgfrb", "Pdgfra", "Vim", "Col1a1")
ECM_Genes_of_Interest <- c("Col1a1", "Col1a2", "Col3a1", "Col4a1","Col5a1", "Col6a1",
                          "Col7a1", "Col8a1", "Col9a1", "Col11a1","Col12a1", "Col15a1",
                          "Eln", "Fn1", "Fmod", "Lama1", "Lama2", "Lama3",
                          "Lama4", "Hspg2")
Mmp_Genes_of_Interest <- c("Mmp2", "Mmp3", "Mmp7", "Mmp8", "Mmp9", 
                           "Mmp11", "Mmp12", "Mmp13", "Mmp14",
                           "Mmp15", "Mmp16", "Mmp17", "Mmp19", 
                           "Mmp23", "Mmp24", "Mmp25", "Mmp28")
Stress_Genes_of_Interest <- c("Atf3", "Jun", "Btg2", "Cdkn1a", "Bad", "Bax")
TGFb_super_Genes_of_Interest <- c("Bmpr1a", "Bmpr1b", "Bmp2", "Bmp4", "Bmp7",
                                  "Tgfb1", "Tgfb2", "Tgfb3",
                                  "Tgfbr1", "Tgfbr2", "Tgfbr3")
Cell_Cycle_Genes_of_Interest <- c("Cenpa", "Cenpe", "Cenpf", "Mki67",
                                  "Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae")
Secreted_Factors_AA <- c("Angpt2", "Bmp2", "Bmp5", "Bmp6", "Ccl2", "Ccl3",
                         "Ccl6", "Ccl9", "Col14a1", "Col4a1", "Cxcl12", "Cxcl13",
                         "Cxcl2", "Dkk2", "Dkkl1", "Edn1", "Egfl7", "Eln",
                         "Esm1", "Fgf1", "Fgf7", "Fgf9", "Fn1", "Gja1",
                         "Gja4", "Gjb2", "Gjc1", "Igf1", "Il1a", "Il1f9",
                         "Il33", "Il6", "Jag1", "Kitl", "Lama4", "Lamb1",
                         "Lamb2", "Lamc1", "Mmp13", "Mmp8", "Mmp9", "Pdgfa",
                         "Pdgfb", "Pdgfc", "Pdgfd", "Tgfb1", "Tgfb2", "Tgfb3",
                         "Tnf", "Vegfa", "Vegfb", "Vegfc", "Vegfb", "Wnt2",
                         "Wnt5a", "Wnt8a")
Surface_Markers_AA <- c("Cd34", "Cd36", "CD37", "CD44", "Cdh5", "Csf1r",
                        "Cxcr4", "Cxcr7", "Dll1", "Dll4", "Dpp4", "Eng",
                        "Flt1", "Fzd4", "Fzd6", "Icam1", "Kdr", "Kit",
                        "Lyve1", "Mcam", "Notch2", "Prominin1", "Robo4",
                        "Sele", "Sell", "Tek", "Tspan7", "Vcam1", "Vwf")
Top100Up_AA  <- c("Acta2", "Mpeg1", "C3", "Cybb", "Axl", "Ctss",
                  "Trp53i11", "Spp1", "Tgfbi", "Tpm2", "Tnc",
                  "Krt14", "Slc7a2", "C1qb", "Mmp14", "Irs2",
                  "Col16a1", "Ms4a7", "Dapk2", "Cxcl2", "Lhfp",
                  "Kcnj2", "Rasgef1b", "Cilp", "Sfrp1", "Pmaip1",
                  "C1qa", "Cdh3", "Mmp12", "Lilrb4", "C1qc", "Clec7a",
                  "Gpr153", "Pld4", "Krt17", "Myl9", "Aspn",
                  "Cx3cr1", "Itga8", "Col17a1", "Cnn1", "Fst",
                  "Postn", "Trf", "Cyp1b1", "Col28a1", "Gas1",
                  "Ncf2", "Lmod1", "Dkk3", "Lamb3", "Kcne4",
                  "Adamtsl1", "Ncf1", "Tenm3", "Tnf", "Fcgr3",
                  "Sdc2", "Pdzrn3", "Igfbp5", "Antxr1", "Abcc3",
                  "Dgkb", "Col1a1", "Il1b", "Pcdh18", "Cxcl14",
                  "Ntrk3", "Abca9", "Igsf10", "Fibin", "Pdgfrl",
                  "Fcgr2b", "Enpp1", "Kcnj8", "Mark1", "Ssc5d",
                  "Myom1", "Scn4b", "Igfbp2", "Inmt", "Trp63",
                  "Fxyd1", "Olfml3", "Cacna1c", "Thbs4", "Cacna2d1",
                  "Icos", "Plxnc1", "Cpxm1", "Col27a1", "Dchs2",
                  "Trim29", "Nr2f1", "Olfr78", "Lum", "Jph2",
                  "Adamts5", "Adcy5", "Il34")
Top100Down_AA  <- c("Bhlha15", "Lpo", "Aqp5", "A630073D07Rik", "Car6", "Derl3",
                    "Prol1", "Oit1", "Gjb1", "5330417C22Rik", "Ppp1r1b", "Prom2",
                    "Agt", "Rap1gap", "Spdef", "Elf5", "Pycr1", "Creb3l1", "Lrrc26",
                    "Aldh3b2", "Scgb2b27", "Galnt3", "Slc2a5", "Scgb1b27", "Ncald",
                    "Wfdc12", "Pip", "Cldn10", "Fbp2", "9130230L23Rik", "Mucl2",
                    "Slc6a14", "Gm47865", "Arfgef3", "AC163018.1", "Scgb2b26",
                    "Arhgef19", "Dnah11", "Slc25a34", "Tex15", "Folr1", "Esp8",
                    "Slc25a48", "4631405J19Rik", "Syne4", "Tpd52l1", "Fkbp11",
                    "Barx2", "Gm45644", "Lman1l", "Gucy2c", "A2ml1", "Fndc5",
                    "Kcne1", "Pex11a", "Muc19", "Slc5a1", "Gbp10", "Ces1e", "Gpd1",
                    "Apol7a", "Gm33586", "Gm43534", "Aldh1a7", "Egf", "Ccl28",
                    "Capn11", "Esp18", "Large2", "Tmem56", "Azgp1", "Ttc25",
                    "Inpp5j", "Gnmt", "Prlr", "Mup6", "Klk1b21", "Slc1a3", "Cyp2s1",
                    "B3galt5", "Tekt5", "A830018L16Rik", "Trpv6", "Brip1", "Smgc",
                    "Kcnk1", "Ttc39aos1", "BC049987", "Mettl21b", "Marveld3",
                    "Klk1b9", "Klk1b22", "Frmpd1", "Pigr", "Klk1b5", "Cgref1",
                    "2010016I18Rik", "Klk1b26", "C2cd4d", "Klk1b8")
Organoid_S_Microarry_Top11_Up <- c("Bpifa2", "Csn3", "Lpo", "Chil1", "Smgc", "Muc19",
                                "Kcnn4", "Crispld2", "Dlk1", "Pnlip", "Aqp5")
Organoid_S_Microarry_Top11_Down <- c("Klk1b9", "Ncam1", "Ntng1", "Vsnl1", "Klk1",
                                     "Klk1b1", "Dsc3", "Htra1", "Il33", "Klk1b11", "Krt5")
Organoid_D_Microarry_Top25_Up <- c("Cpa3", "Slc14a1", "Glp1r", "Cma1", "Tpsb2",
                                   "Dkk2", "2010007H06Rik", "Car12", "Ccdc3", "Serpina3n",
                                   "Lcn2", "Frzb", "Hhip", "Ppbp", "Ano3",
                                   "Nefl", "Mcpt4", "Cd34", "Fmo2", "Csn1s1",
                                   "Tm4sf1", "Bmp7", "Ereg", "Crip2", "Gm2115")
Organoid_D_Microarry_Top25_Down <- c("Abi3bp", "Tgfb3", "Npr3", "Serpine1", "Adamts6", 
                                     "Srpx", "Adamts5","Fat4", "Lox", "Itga11",
                                     "Ndufa4l2", "Podn", "Ncam1","Mfap5", "Cnn1", 
                                     "Ccn5", "Fbn2", "Pappa", "Avpr1a","Epha3",
                                     "Rspo3", "Col12a1", "Ccn4", "Vegfd", "Fmod")
Organoid_D_Microarry_SecretedSignals_Up <- c("Areg", "Bmp2", "Bmp7", "Chil1", "Deptor", "Ereg",
                                             "Frzb", "Grem1", "Hhip", "Igfbp5", "Kitl", "Pgf",
                                             "Sema3a", "Sema3b", "Sema6a", "Sema6d", "Spry2",
                                             "Spry4", "Tgfbi")
Organoid_D_Microarry_SecretedSignals_Down <- c("Bmp3", "Bmpr1b", "Fgf18", "Fgf2", "Fgf5",
                                               "Fzd1", "Fzd4", "Nrg1", "Rspo3", "Sema3c",
                                               "Slit2", "Slit3", "Tgfb1i1", "Tgfb2", "Tgfb3",
                                               "Vegfa", "Vegfc", "Vegfd", "Ccn4", "Ccn5", "Wls")
All_Up_AA  <- c("Acta2", "Mpeg1", "C3", "Cybb", "Axl", "Ctss", "Trp53i11",
                "Spp1", "Tgfbi", "Tpm2", "Tnc", "Krt14", "Slc7a2", "C1qb",
                "Mmp14", "Irs2", "Col16a1", "Ms4a7", "Dapk2", "Cxcl2", "Lhfp",
                "Kcnj2", "Rasgef1b", "Cilp", "Sfrp1", "Pmaip1", "C1qa", "Cdh3",
                "Mmp12", "Lilrb4", "C1qc", "Clec7a", "Gpr153", "Pld4", "Krt17",
                "Myl9", "Aspn", "Cx3cr1", "Itga8", "Col17a1", "Cnn1", "Fst",
                "Postn", "Trf", "Cyp1b1", "Col28a1", "Gas1", "Ncf2", "Lmod1",
                "Dkk3", "Lamb3", "Kcne4", "Adamtsl1", "Ncf1", "Tenm3", "Tnf",
                "Fcgr3", "Sdc2", "Pdzrn3", "Igfbp5", "Antxr1", "Abcc3", "Dgkb", 
                "Col1a1", "Il1b", "Pcdh18", "Cxcl14", "Ntrk3", "Abca9", "Igsf10",
                "Fibin", "Pdgfrl", "Fcgr2b", "Enpp1", "Kcnj8", "Mark1", "Ssc5d",
                "Myom1", "Scn4b", "Igfbp2", "Inmt", "Trp63", "Fxyd1", "Olfml3",
                "Cacna1c", "Thbs4", "Cacna2d1", "Icos", "Plxnc1", "Cpxm1",
                "Col27a1", "Dchs2", "Trim29", "Nr2f1", "Olfr78", "Lum", "Jph2",
                "Adamts5", "Adcy5", "Il34", "Tgfb3", "Capn6", "Cd33", "Col1a2",
                "Slc11a1", "Pi15", "Fgf7", "Pou2f2", "Tmem119", "Ltbp2", "Col6a3", 
                "Wisp2", "Mrvi1", "Crispld2", "Hpgds", "Chil1", "Pgf", "Fbln1",
                "Svep1", "Loxl1", "Podn", "Des", "Gpc6", "Barx1", "Pdpn", "Csf3r",
                "Basp1", "Il17a", "Col3a1", "Ogn", "Il33", "Fzd2", "Emilin2",
                "Itga7", "Dclk1", "Tspan11", "Ctsk", "C1s1", "Cacna1g", "Oaf",
                "Fcgr1", "Lacc1", "Slc4a8", "Olfr558", "Tacstd2", "5830411N06Rik",
                "Robo1", "Srpx2", "Naalad2", "Olfml1", "Pde1b", "Wisp1", "Mfap4",
                "Kcnk2", "Vgll3", "Scml4", "F3", "Col6a2", "Mfap5", "Il1a", "Nlrp3",
                "Trac", "Heyl", "Aldh1a1", "Clec12a", "Adap2", "Mrc2", "Rnf150",
                "Gm6377", "Col5a1", "Il1rl2", "Ccbe1", "Sult5a1", "Gpx3", "Rab7b",
                "Akr1b8", "Stc2", "Kit", "Siglec1", "Efemp1", "Rbpms2", "Susd5",
                "Myh11", "Lpl", "Cacna1h", "Cdh11", "Nrip2", "Trem2", "Rasl11b",
                "Col4a6", "Col6a6", "Col6a1", "Cd8a", "Osm", "Rasl12", "Pdgfra",
                "Clec11a", "Lrp1", "Itgbl1", "Cd300c2", "Sphk1", "Alx4", "Twist1",
                "Map3k7cl", "Cdh19", "Synpo2", "Abca6", "Atp1a2", "Lamc3", "Plxdc2",
                "Lpar1", "Tlr13", "Lcn2", "Havcr2", "Gdf6", "C3ar1", "Inhba",
                "Ptges", "Serpina3n", "A630072L19Rik", "Sod3", "Osr1", "Hmgcs2",
                "Egflam", "Rarres2", "Npy1r", "Flnc", "Smoc2", "Kcnab1", "Col14a1",
                "Robo3", "Islr", "Ms4a6d", "Lama2", "Scn7a", "Ptgs2", "Loxl3",
                "Maf", "Dnm1", "Cdh6", "Fbln7", "Prelp", "Fcrls", "Apod", "Ildr2",
                "Col8a2", "Rbm24", "Lgi2", "Spon1", "Lsamp", "Omd", "Ggt5", "Hgf",
                "Dcn", "Cd28", "Cd3d", "Tbx3os1", "Ms4a14", "Astn2", "Hcar2",
                "Serping1", "Lyz2", "Mmp19", "Col4a5", "Trpc6", "Xirp1", "Cygb",
                "Zfhx4", "Nox4", "F2rl3", "Cd6", "Dio2", "Abca8b", "Clmp",
                "Ndufa4l2", "Nkd1", "Ednra", "Ephb3", "Wnt9a", "Ror1", "Cd83",
                "Clec4n", "Col5a3", "Rgs4", "Pla2g4a", "Prss12", "Disp2", "Mmp2",
                "Fhl2", "Pamr1", "Wnt5b", "Col6a5", "C5ar1", "Prrx1", "Kcnq4",
                "Cfh", "Cd276", "Scn2b", "Sorcs2", "Sobp", "Gsn", "Fbn2", "Cacng7",
                "Medag", "A830082K12Rik", "Clec5a", "Plin4", "Cap2", "Colec12",
                "Adamts8", "Lrrc15", "Cd8b1", "Slc6a17", "Dusp10", "Ptgis",
                "Gabra3", "Higd1b", "Epha3", "Siglece", "Dzip1l", "P2ry14",
                "D630003M21Rik", "Tbx15", "Frem1", "Adamts12", "Pla1a", "Lin7a",
                "Avpr1a", "Pdgfrb", "Cacna1d", "Htra3", "Gpnmb", "Mt2", "Abca8a",
                "Fgf10", "Wif1", "Lrrc25", "Hc", "Abcd2", "Serpinf1", "Cd40lg",
                "Wdr66", "C1qtnf7", "Gpr20", "P2ry1", "Pde1a", "Fam169b", "Blk",
                "Mrgprf", "Tcrg-C2", "Gm42793", "Col5a2", "Ebf2", "Itga11", "Entpd2",
                "Lrrn2", "Shc4", "Tacr1", "Pla2g7", "Fam180a", "Fmod", "Cacnb2",
                "Pth1r", "Il1r2", "Ccl12", "Gm45714", "Gucy1b1", "Tnip3", "Clec3b",
                "B3gnt9", "Igsf9b", "Tmem200a", "Siglecg", "Ccl2", "Tmem132c",
                "Fcrl1", "1500015O10Rik", "Mmp13", "Pygo1", "Slc36a2", "Spsb1",
                "Sectm1b", "Foxs1", "Art4", "Sfrp2", "Cacna1e", "Dpep1", "Scara5",
                "Rai2", "Aff3", "Gm11639", "Cd79b", "Slamf6", "Ldb3", "Bmper",
                "Cd5", "Rubcnl", "BC067074", "Irf4", "Slc38a11", "Ntng1", "Ror2",
                "Tnfaip6", "Gpm6b", "Iglon5", "Klhl29", "Spon2", "Kcnmb1", "Osr2",
                "Srpx", "Ripor3", "Trnp1", "Serpine1", "Krt6a", "Anxa8", "Zim1",
                "Zfp469", "Abi3bp", "Gpr176", "Clec4a2", "Dnase1l3", "Sema3a",
                "Myocd", "Has2", "Gdf3", "Il1rn", "C2", "Cpz", "Ccr6", "Krt19",
                "Ccdc80", "Wdr86", "I830077J02Rik", "Actg2", "Nfasc", "Nrxn2",
                "Cspg4", "Fgf18", "Slfn1", "Abcc9", "Dio3", "Gpc3", "Cercam",
                "Gm48878", "Msr1", "Pnmal2", "Tspoap1", "Trdv4", "Cxcl5", "Crabp1",
                "Adamts17", "Grem2", "H2-M2", "5430437J10Rik", "Cys1", "Nlrp1c-ps",
                "Lgi4", "Tpbg", "Ncs1", "Has1", "Hmcn2", "Hrc", "Bgn", "Lef1",
                "Slc8a2", "Inhbb", "Tfap2c", "Csdc2", "Pou2af1", "Rgs6", "S100b",
                "Fgf14", "Kcp", "Rad51b", "Timp1", "Slitrk6", "Cd300lf", "Bnc2",
                "Adcyap1r1", "Cacna1i", "4833422C13Rik", "AC121574.1", "Enox1",
                "Pnck", "Treml2", "Mmp3", "Sh3rf3", "Serpinb5", "Clec4e", "Igfbp6",
                "Itih2", "Galnt17", "Pcp4", "Ccr4", "Astn1", "Dlx2", "Gldn", "Npl",
                "Sh2d5", "Ccr7", "Gm13861", "Arhgap22", "Wnt7b", "Ighd", "Prkg2",
                "Sdk1", "Trpc4", "B4galnt4", "Acod1", "Brinp3", "Npas4", "Msc",
                "Gp2", "Sit1", "Fam19a5", "3632451O06Rik", "Fam26e", "Jakmip1",
                "Trbv13-3", "Sell", "4933424M12Rik", "Cxcr5", "Dennd2a", "Fcer2a",
                "Zfp536", "Gm34829", "B430306N03Rik", "Nxph3", "Dpt", "Eya4",
                "Pdcd1lg2", "Lipc", "Steap2", "Lingo1", "Ky", "Gm5086", "Il23r",
                "4933431K23Rik", "Srpk3", "Lrrc17", "Clca3a2", "Gcnt4", "Gm9176",
                "Cd22", "Trhde", "Cd19", "Pax5", "Grp", "Gm5834", "Gm43145",
                "Izumo1r", "Gli2", "Gm8369", "Cxcr2", "Ntm", "Btla", "Bank1",
                "Iglc2", "E030013I19Rik", "Dppa3", "Tfec", "1500009L16Rik", "Otogl",
                "Gucy2g", "Tmem132e", "Dmp1", "Kcnc1", "Rasgrf1", "Acss3", "Cr2",
                "Crhr2", "Vmn2r97", "D030007L05Rik", "Gm48727", "Gm34354", "Il13ra2",
                "Gm14023", "Fam26f", "Adgrv1", "Fcmr", "Frzb", "Begain", "Ephb2",
                "Clec4d", "Cd1d2", "Serpina1b", "Itih3", "Arsi", "Drd1", "Aff2",
                "Col26a1", "Tcf24", "Cemip", "S100a8", "Arl11", "Gm15684", "Xpnpep2",
                "Stra6", "C230012O17Rik", "Prrx2", "Tmem45a", "Ryr2", "Casq2",
                "Tnfrsf8", "Cd79a", "Pcyt1b", "Gm47015", "9530034E10Rik",
                "4932435O22Rik", "Muc4", "Olfm2", "Gm13881", "H2-DMb2", "Mab21l2",
                "Trbv16", "Lrrc10b")
All_Down_AA  <- c("Bhlha15", "Lpo", "Aqp5", "A630073D07Rik", "Car6", "Derl3",
                  "Prol1", "Oit1", "Gjb1", "5330417C22Rik", "Ppp1r1b", "Prom2",
                  "Agt", "Rap1gap", "Spdef", "Elf5", "Pycr1", "Creb3l1",
                  "Lrrc26", "Aldh3b2", "Scgb2b27", "Galnt3", "Slc2a5",
                  "Scgb1b27", "Ncald", "Wfdc12", "Pip", "Cldn10", "Fbp2",
                  "9130230L23Rik", "Mucl2", "Slc6a14", "Gm47865", "Arfgef3",
                  "AC163018.1", "Scgb2b26", "Arhgef19", "Dnah11", "Slc25a34",
                  "Tex15", "Folr1", "Esp8", "Slc25a48", "4631405J19Rik", "Syne4",
                  "Tpd52l1", "Fkbp11", "Barx2", "Gm45644", "Lman1l", "Gucy2c",
                  "A2ml1", "Fndc5", "Kcne1", "Pex11a", "Muc19", "Slc5a1",
                  "Gbp10", "Ces1e", "Gpd1", "Apol7a", "Gm33586", "Gm43534",
                  "Aldh1a7", "Egf", "Ccl28", "Capn11", "Esp18", "Large2",
                  "Tmem56", "Azgp1", "Ttc25", "Inpp5j", "Gnmt", "Prlr", "Mup6",
                  "Klk1b21", "Slc1a3", "Cyp2s1", "B3galt5", "Tekt5",
                  "A830018L16Rik", "Trpv6", "Brip1", "Smgc", "Kcnk1",
                  "Ttc39aos1", "BC049987", "Mettl21b", "Marveld3", "Klk1b9",
                  "Klk1b22", "Frmpd1", "Pigr", "Klk1b5", "Cgref1", "2010016I18Rik",
                  "Klk1b26", "C2cd4d", "Klk1b8", "4833423E24Rik", "Adtrp", "Klk14",
                  "Rab26", "Adora1", "Frmpd1os", "Tulp1", "Slc22a3", "Kndc1",
                  "Klk1b4", "Gabra4", "Asgr2", "Klk10", "Abo", "Sval2", "Gm5886",
                  "Lmx1b", "Gdf5", "Jchain", "Vip", "Gm5546", "Sypl2", "D7Ertd443e",
                  "1810019D21Rik", "2810442N19Rik", "Igsf5", "Erbb4", "Ddo",
                  "BC024139", "Hsbp1l1", "Dcpp3", "Klk1b11", "Dcpp2", "Bpifa2",
                  "AI463170", "1700019D03Rik", "1810062G17Rik", "Dbh", "Gm5737",
                  "Ccdc162", "Tmcc3os", "8430419K02Rik", "Muc13", "Apof", "Ptprq",
                  "Folh1", "Gm10863", "Dcpp1", "Mup4", "Gdpd2", "Smim22", "Gm38244",
                  "Klk1b16", "Klk1b27", "9530002B09Rik", "Ugt8a", "Gm10790", "Zan",
                  "Cd164l2", "Gm14133", "Scgb2b24", "Gm13421", "Uox", "Tmem139",
                  "Mup5", "Gm20275", "Nxf7", "Gm12888", "Esp4", "Lgals12", "Itgb2l",
                  "Hapln4", "Atp13a4", "2810030D12Rik", "1700003E16Rik",
                  "B230206L02Rik", "C130021I20Rik", "Ly6g6e", "Cspg5", "Gabrb3",
                  "Htr3a", "Gm16025", "Gm41177", "Gm37737", "Cyp26b1", "Hhatl",
                  "Slc7a14", "Crisp3", "Klk1b1", "AI463229", "4930517G19Rik",
                  "Neu2", "Gm13648", "1700003D09Rik", "Klk1b3", "Klk1", "Ceacam10",
                  "Grin1os", "Gm20554", "Cyp1a1")
Stromal_Transcription_Factors <- c("Atf2",
                                  "Pelp1",
                                  "Nfkb1",
                                  "Ets1",
                                  "Dmtf1",
                                  "Egr1",
                                  "Myocd",
                                  "Gtf2i",
                                  "Jund",
                                  "Srebf2",
                                  "Erg",
                                  "E2f1",
                                  "Rbl2",
                                  "Trp53",
                                  "Stat3",
                                  "Sp1")
#https://science.sciencemag.org/content/sci/371/6534/eabc3172.full.pdf
Growth_Factors_Receptors_Zepp2021 <- c("Pdgfa",
                                       "Pdgfb",
                                       "Wnt10a",
                                       "Wnt6",
                                       "Wnt5b",
                                       "Wnt4",
                                       "Fzd1",
                                       "Fzd2",
                                       "Fzd4",
                                       "Lrp1",
                                       "Lrp6",
                                       "Gli3",
                                       "Smo",
                                       "Vegfb",
                                       "Mfap5")
Genes_Mascharak2021 <- c("Msn",
                         "Dlk1",
                         "Tead1",
                         "Tead2",
                         "Tead3",
                         "Tead4",
                         "Yap1",
                         "Wwtr1",
                         "Fat4",
                         "Dchs1",
                         "Dchs2")
#https://advances.sciencemag.org/content/7/24/eabg6005?elqTrackId=5ca8a69ac8ff4f3080e47a4a072a3d66&elq=6c66d2cfe7a44810b72937392ed0cbc5&elqaid=31499&elqat=1&elqCampaignId=10598
Genes_Xie2021 <- c("Eln", "Cxcl14", "Npnt", "Fn1", "Hspb1", "Xist", "Nebl", "Tagln", "Cd74", "Cxcl1", "Cxcl2", "Cxcl5", "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl16")
FGFR_Genes <- c("Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4", "Fgfrl1")
#https://www.nature.com/articles/s41586-020-2941-1
Myofibroblast_Genes <- c("Notch3", "Rgs5", "Meg3", "Colec11", "Cxcl12")

#Loop to create Feature plots and violin plots for Gene lists
for(i in Myofibroblast_Genes[1:5]){
  print(i)
  plot3 <- FeaturePlot(agg.E16, features = i, min.cutoff = "q9")
  plot4 <- VlnPlot(agg.E16, features = i) + NoLegend()
  Dot_Plot_Title <- str_c("UMAP_Plot_", i, "_E16_04.pdf")
  Violin_Plot_Title <- str_c("Violin_Plot_", i, "_E16_04.pdf")
  pdf(file = Dot_Plot_Title)
      print(plot3)
  dev.off()  
  pdf(file = Violin_Plot_Title)
      print(plot4)
  dev.off()  
}

#Make double feature plots for markers that shows overlap of Gli1 with Pdgfra
pdf("UMAP_E16_Pdgfra_Gli1_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16, features = c("Pdgfra", "Gli1"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Green","Red"))
dev.off()
pdf("UMAP_E16_Pdgfrb_Gli1_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16, features = c("Pdgfrb", "Gli1"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Blue","Red"))
dev.off()
pdf("UMAP_E16_Pdgfra_Gli3_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16, features = c("Pdgfra", "Gli3"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Green","Red"))
dev.off()
pdf("UMAP_E16_Pdgfrb_Gli3_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16, features = c("Pdgfrb", "Gli3"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Blue","Red"))
dev.off()

############################## Differential gene list creation #################################
################################################################################################

  #Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.E16.markers.All <- FindAllMarkers(agg.E16, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.All, file = "All_Gene_Averages_perCluster_E16_04.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos <- FindAllMarkers(agg.E16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.Pos, file = "Positive_Gene_Averages_perCluster_E16_04.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos.Top10 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.E16.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perCluster_E16_04.csv")
#Reduces to top 5 postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos.Top5 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.csv(agg.E16.markers.Pos.Top5, file = "Positive_Top5_Gene_Averages_perCluster_E16_04.csv")
#Finds all postive marker genes for cluster #13 compared to all other clusters
agg.E16.markers.Cluster13 <- FindMarkers(agg.E16, ident.1 = 13, min.pct = 0.25)
write.csv(agg.E16.markers.Cluster13, file = "All_Gene_Averages_Cluster13_E16_04.csv")

#Attempt to find specific markers in each cluster
agg.E16.markers.All <- FindAllMarkers(agg.E16.2, min.pct = 0.25, min.diff.pct = 0.25, logfc.threshold = log(1.5))
write.csv(agg.E16.markers.All, file = "All_Specific_Genes_Averages_perCluster_E16_04.csv")
agg.E16.markers.All.Top30 <- agg.E16.markers.All %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


    #Naming clusters
#Makes a copy of original Total Seurat data object
agg.E16.2 <- agg.E16
#Old list with cluster numbers
old.cluster.ids <- c("0", "1", "2", "3", "4",
                     "5", "6", "7", "8", "9",
                     "10", "11", "12", "13", "14",
                     "15", "16", "17", "18", "19",
                     "20", "21", "22", "23", "24",
                     "25", "26", "27", "28", "29",
                     "30", "31", "32", "33", "34",
                     "35", "36", "37", "38")
  #Making Broad Cellular labels
#New list with simple cluster names
new.cluster.ids <- c(
                    "Stroma",
                    "Stroma",
                    "Stroma",
                    "Stroma",
                    "Endothelium",
                    "Endothelium",
                    "Endothelium",
                    "Stroma",
                    "Stroma",
                    "Progenitor",
                    "Epithelium",
                    "Immune",
                    "Stroma",
                    "Lymphatic_Endothelium",
                    "Endothelium",
                    "Progenitor",
                    "Stroma",
                    "Epithelium",
                    "Immune",
                    "Epithelium",
                    "Epithelium",
                    "Erythrocytes",
                    "Stroma",
                    "Schwann_Cells",
                    "Epithelium",
                    "Immune",
                    "Immune",
                    "Stroma",
                    "Endothelium",
                    "Lymphatic_Endothelium",
                    "Erythrocytes",
                    "Endothelium",
                    "Epithelium",
                    "Stroma",
                    "Epithelium",
                    "Immune",
                    "Neurons",
                    "Stroma",
                    "Lymphatic_Endothelium")
#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.E16.2)
agg.E16.2 <- RenameIdents(agg.E16.2, new.cluster.ids)
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Broad_Names_04.csv")
#Draws a new UMAP plot with written labels from new.cluster.ids
pdf("UMAP_E16_labeled_Simple_04.pdf", width = 8, height = 6)
  DimPlot(agg.E16.2, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.5, repel = TRUE) +
  NoLegend() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()
#Genes that show up in broad cell types
Broad_Label_Genes  <- c("Col1a1",
                        "Dlk1",
                        "Kdr",
                        "Cdh5",
                        "Epcam",
                        "Sox9",
                        "Ptprc",
                        "Cd37",
                        "C1qa",
                        "C1qb",
                        "Hba-a1",
                        "Hba-a2",
                        "Plp1",
                        'Ptprz1',
                        "Prph",
                        "Tubb3")
#Draws a Dot plot for genes of interest for Broad labels
pdf("DotPlot_GeneIDs_Simple_Labeled_E16_04.pdf", width = 8, height = 4)
    DotPlot(agg.E16.2, features = Broad_Label_Genes, cols = "Spectral", col.min = 0, dot.scale = 5) +
    RotatedAxis()
dev.off()

#Remakes the copy of original Total Seurat data object
agg.E16.2 <- agg.E16
  #New list with Specific cluster names
new.cluster.ids <- c("Stroma_Pdgfra_Pdgfrb_01", "Stroma_Pdgfra_Pdgfrb_02", "Stroma_Pdgfra_Pdgfrb_03", "Stroma_Pdgfra_Pdgfrb_Proliferative_01", "Endothelium_Proliferative_01",
                     "Stroma_Proliferative_01", "Endothelium_01", "Stroma_Pdgfrb_01", "Stroma_Pdgfrb_Proliferative_01", "Unknown_01",
                     "Myoepithelium_01", "Immune_CD45+_01", "Stroma_01", "Lymphatics_Endothelium_01", "Stroma_02",
                     "Unknown_02", "Stroma_Thy1_01", "Myoepithelium_02", "Immune_Macrophage_01", "Acinar_Smgc+_01",
                     "Acinar_Psp+_01", "Erythrocytes_01", "Stroma_Pdgfra_01", "Schwann_Cells_01", "Ducts_K19+_01",
                     "Immune_Memory_CD4+_01", "Immune_Mast_Cells_01", "Stroma_Pdgfrb_Proliferative_02", "Endothelium_02", "Lymphatics_Endothelium_Proliferative_01",
                     "Erythrocytes_02", "Endothelium_03", "Ducts_K5+_01", "Stroma_Pdgfrb_aSMA_01", "Myoepithelium_Proliferative_01",
                     "Immune_B_Cells_01", "Neurons_01", "Stroma_Pdgfrb_Proliferative_03", "Lymphatic_Endothelium_02")
new.cluster.ids <- c("01_Stroma_Pa_Pb",
                     "02_Stroma_Pa_Pb",
                     "03_Stroma_Pa_Pb",
                     "04_Stroma_Pa_Pb_CC",
                     "05_Endo_CC",
                     "06_Endo_CC",
                     "07_Endo",
                     "08_Stroma_Pb",
                     "09_Stroma_CC",
                     "10_Stromal Progenitor",
                     "11_Myoepi",
                     "12_Immune_CD45",
                     "13_Stroma",
                     "14_Lym_Endo",
                     "15_Endo",
                     "16_Stromal_Progenitor",
                     "17_Stroma_Thy1",
                     "18_Myoepi",
                     "19_Immune_Macro",
                     "20_Proacinar_Smgc",
                     "21_Proacinar_PSP",
                     "22_Erythro",
                     "23_Stroma_Pa",
                     "24_Schwann_Cells",
                     "25_Duct_K19",
                     "26_Immune_CD4",
                     "27_Immune_Mast",
                     "28_Stroma_Pb",
                     "29_Endo",
                     "30_Lym_Endo_CC",
                     "31_Erythro",
                     "32_Endo",
                     "33_Duct_K5",
                     "34_Stroma_Pb_SMA",
                     "35_Myoepi_CC",
                     "36_Immune_B_cells",
                     "37_Neurons",
                     "38_Stroma_Pb_CC",
                     "39_Lym_Endo")
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Names_04.csv")
#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.E16.2)
agg.E16.2 <- RenameIdents(agg.E16.2, new.cluster.ids)
#Draws a new UMAP plot with written labels from new.cluster.ids
pdf("UMAP_E16_labeled_04.pdf", width = 12, height = 8)
  DimPlot(agg.E16.2, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()
pdf("UMAP_E16_labeled_key_04.pdf", width = 16, height = 8)
  DimPlot(agg.E16.2, reduction = "umap", label = TRUE, label.size = 1.5, pt.size = 0.5) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()

  #Dot plot of cell-type defining markers
#List of genes to graph
Genes_of_Interest <- c("Col1a1", "Vim", "Pdgfra", "Pdgfrb", "Thy1", "Fgf2", "Acta2", "Epcam", "Krt5", "Krt14", "Krt19", "Aqp5", "Bpifa2", "Smgc", "Cdh5", "Kdr", "Pdpn", "Lyve1", "Tubb3", "Hba-a1", "Hba-a2", "Il7r", "Ms4a2", "Ptprc", "Mki67", "Top2a")
#Draws a Dot plot for genes of interest for broad labels
pdf("DotPlot_GeneIDs_Broad_Labeled_E16_All_04.pdf", width = 9, height = 5)
  DotPlot(agg.E16.2, features = Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Specific_Labeled_E16_All_04.pdf", width = 18, height = 10)
  DotPlot(agg.E16.2, features = Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for ECM genes of interest for Specific labels
pdf("DotPlot_GeneIDs_ECM_Specific_Labeled_E16_All_04.pdf", width = 12, height = 10)
  DotPlot(agg.E16.2, features = ECM_Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Mmp genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Mmp_Specific_Labeled_E16_All_04.pdf", width = 10, height = 10)
  DotPlot(agg.E16.2, features = Mmp_Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Tgfb genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Tgfb_All_Specific_Labeled_E16_All_04.pdf", width = 8, height = 10)
  DotPlot(agg.E16.2, features = TGFb_super_Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA Secreted Factors genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_SecretedFactors_All_Specific_Labeled_E16_All_04.pdf", width = 15, height = 10)
  DotPlot(ag.E16.2, features = Secreted_Factors_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA Secreted Factors genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_SurfaceMarkers_All_Specific_Labeled_E16_All_04.pdf", width = 10, height = 10)
  DotPlot(agg.E16.2, features = Surface_Markers_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA Top 100 upregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_Top100Up_All_Specific_Labeled_E16_All_04.pdf", width = 25, height = 10)
  DotPlot(agg.E16.2, features = Top100Up_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA All upregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_AllUp_All_Specific_Labeled_E16_All_04.pdf", width = 125, height = 10)
  DotPlot(agg.E16.2, features = All_Up_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA Top 100 Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_Top100Down_All_Specific_Labeled_E16_All_04.pdf", width = 25, height = 10)
  DotPlot(agg.E16.2, features = Top100Down_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for AA All Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_AA_AllDown_All_Specific_Labeled_E16_All_04.pdf", width = 125, height = 10)
  DotPlot(agg.E16.2, features = All_Down_AA, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Epitheial S microarray Top 11 upregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_S_Microarray_Top11Up_All_Specific_Labeled_E16_All_04.pdf", width = 8, height = 10)
  DotPlot(agg.E16.2, features = Organoid_S_Microarry_Top11_Up, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Epitheial S microarray Top 11 Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_S_Microarray_Top11Down_All_Specific_Labeled_E16_All_04.pdf", width = 8, height = 10)
  DotPlot(agg.E16.2, features = Organoid_S_Microarry_Top11_Down, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Stroma D microarray Top 25 upregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_D_Microarray_Top25Up_All_Specific_Labeled_E16_All_04.pdf", width = 11, height = 10)
  DotPlot(agg.E16.2, features = Organoid_D_Microarry_Top25_Up, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Stroma D microarray Top 25 Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_D_Microarray_Top25Down_All_Specific_Labeled_E16_All_04.pdf", width = 11, height = 10)
  DotPlot(agg.E16.2, features = Organoid_D_Microarry_Top25_Down, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Stroma D microarray Secreted Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_D_Microarray_SecretedSignalsUp_All_Specific_Labeled_E16_All_04.pdf", width = 10, height = 10)
  DotPlot(agg.E16.2, features = Organoid_D_Microarry_SecretedSignals_Up, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for Organoid Stroma D microarray Secreted Downregulated genes of interest for Specific labels
pdf("DotPlot_GeneIDs_Oragnoid_D_Microarray_SecretedSignalsDown_All_Specific_Labeled_E16_All_04.pdf", width = 10, height = 10)
  DotPlot(agg.E16.2, features = Organoid_D_Microarry_SecretedSignals_Down, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
#Draws a Dot plot for genes that define stroma in lung based on Zepp 2021 with  Specific labels
pdf("DotPlot_GeneIDs_Zepp2021_Specific_Labeled_E16_All_04.pdf", width = 10, height = 10)
DotPlot(agg.E16.2, features = Growth_Factors_Receptors_Zepp2021, cols = "Spectral", col.min = 0, dot.scale = 5) +
  RotatedAxis()
dev.off()

#Heatmap plot of Top 10 genes defining cell clusters
pdf("Heatmap_Top10_GeneIDs_E16_04.pdf", width = 60, height = 30)
    DoHeatmap(agg.E16.2, features = agg.E16.markers.Pos.Top10$gene)
dev.off()
#Heatmap plot of Top 5 genes defining cell clusters
pdf("Heatmap_Top5_GeneIDs_E16_04.pdf", width = 40, height = 30)
    DoHeatmap(agg.E16.2, features = agg.E16.markers.Pos.Top5$gene)
dev.off()


  #Stroma subset analysis Without re-calculating the data
  #Pulls out the stromal data alone for seperate analyses
#Creates a subset of agg.E16 including all stromal cells
agg.E16.Stroma <- subset(agg.E16.2, idents = c("Stroma_Pdgfra_Pdgfrb_01",
                                               "Stroma_Pdgfra_Pdgfrb_02",
                                               "Stroma_Pdgfra_Pdgfrb_03",
                                               "Stroma_Pdgfra_Pdgfrb_Proliferative_01",
                                               "Stroma_Proliferative_01",
                                               "Stroma_Pdgfrb_01",
                                               "Stroma_Pdgfrb_Proliferative_01",
                                               "Stroma_01",
                                               "Stroma_02",
                                               "Stroma_Thy1_01",
                                               "Stroma_Pdgfra_01",
                                               "Stroma_Pdgfrb_Proliferative_02",
                                               "Stroma_Pdgfrb_aSMA_01",
                                               "Stroma_Pdgfrb_Proliferative_03"))
#UMAP plot for only the stromal clusters
pdf("UMAP_E16_Stroma_Labeled_key_04.pdf", width = 12, height = 8)
  DimPlot(agg.E16.Stroma, reduction = "umap", label = FALSE, label.size = 6, pt.size = 0.5, repel = TRUE) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()

#Make double feature plots for markers that separate stromal clusters
pdf("UMAP_E16_Stroma_Thy1_Pdgfra_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16.Stroma, features = c("Thy1", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Magenta", "Green"))
dev.off()
pdf("UMAP_E16_Stroma_Pdgfrb_Pdgfra_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16.Stroma, features = c("Pdgfrb", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","Green"))
dev.off()
pdf("UMAP_E16_Stroma_Pdgfrb_aSMA_04.pdf", width = 40, height = 8)
  FeaturePlot(agg.E16.Stroma, features = c("Pdgfrb", "Acta2"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
dev.off()

#Draws a Dot plot for genes of interest for Specific labels
pdf("DotPlot_Stroma_Pdgfra_FGF2_Specific_Labeled_E16_All_04.pdf", width = 7, height = 4)
  DotPlot(agg.E16.Stroma, features = c("Fgf2", "Acta2", "Pdgfrb", "Pdgfra", "Thy1"), cols = "Spectral", col.min = 0, dot.scale = 5) +
    RotatedAxis()
dev.off()

#Finds all marker genes for every cluster compared to all other clusters
agg.E16.markers.All <- FindAllMarkers(agg.E16.Stroma, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.All, file = "All_Gene_Averages_perStromalCluster_E16_04.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos <- FindAllMarkers(agg.E16.Stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.Pos, file = "Positive_Gene_Averages_perStromalCluster_E16_04.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos.Top10 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.E16.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perStromalCluster_E16_04.csv")

#Heatmap plot of Top 10 genes defining only the stromal clusters
pdf("Heatmap_Top10_Stroma_E16_04.pdf", width = 60, height = 30)
  DoHeatmap(agg.E16.Stroma, features = agg.E16.markers.Pos.Top10$gene)
dev.off()


############################## Epithelium subset analysis ######################################
################################################################################################
 
  #Pulls out the epithelial data and re-normalized the cells to only epithelial cells
#Creates a subset of agg.E16 based on Epcam expression, epithelial cells
agg.E16.Epi <- subset(agg.E16.2, idents = c("Epithelium"))
# Check that the subset took the right subpopulations
DimPlot(agg.E16.Epi, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
#Labels all genes
all.genes <- rownames(agg.E16.Epi)
#Re-applies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.E16.Epi <- ScaleData(agg.E16.Epi, feature = all.genes) 
#Re-runs PCA statistics used for generates a PCA plot of the first PCs
agg.E16.Epi <- RunPCA(agg.E16.Epi, features = VariableFeatures(object = agg.E16))
#generates list of Top 10 most variable genes
top10.Epi <- head(VariableFeatures(agg.E16.Epi), 10)
#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
pdf("Elbowplot_E16_Epi_04.pdf")
ElbowPlot(agg.E16.Epi, ndims = 50)
dev.off()
#Reruns for constructing a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset (first 40 PCs, dims = 40), how many dimentions do you want included for the analysis.
agg.E16.Epi <- FindNeighbors(agg.E16.Epi, dims = 1:21) 
#Reruns generation of the Clusters, resolution increases the specificity of each cluster, resolution = 2.0, the default is 1.0
agg.E16.Epi <- FindClusters(agg.E16.Epi, resolution = 1.1)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.E16.Epi <- RunTSNE(agg.E16.Epi, dims = 1:21)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.E16.Epi <- RunUMAP(agg.E16.Epi, dims = 1:21)

#Saves Formated dataset after all major computational processing
saveRDS(agg.E16.Epi, file = "E16_SMG_Raw_Epi_(SEURAT_v4)_04.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Epi_only")
agg.E16.Epi <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Epi_only/E16_SMG_Raw_Epi_(SEURAT_v4)_04.rds")

#Creates UMAP with cluster numbers
pdf("UMAP_E16_Epi_ClusterNumbers_04.pdf", width = 6.25, height = 5)
  DimPlot(agg.E16.Epi, reduction = "umap", label = TRUE, pt.size = 1.5, group.by = "seurat_clusters") + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()
#Cluster names for Epithelium
new.cluster.ids <- c("Myo.Epi_01",
                     "Myo.Epi_02",
                     "Myo.Epi_CC",
                     "Duct_K19",
                     "Myo.Epi_02",                     
                     "Proacinar_PSP",
                     "Proacinar_Smgc_CC",
                     "Proacinar_Smgc_01",
                     "Proacinar_Smgc_02",
                     "Duct_K5",
                     "Myo.Epi_CC")
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Names_Epi_04.csv")
#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.E16.Epi)
agg.E16.Epi <- RenameIdents(agg.E16.Epi, new.cluster.ids)
#Creates UMAP with cluster labels
pdf("UMAP_E16_Epi_Labels_04.pdf", width = 8, height = 5)
    DimPlot(agg.E16.Epi, reduction = "umap", label = TRUE, label.size = 6, pt.size = 1.5, repel = TRUE)  + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()
#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.E16.Epi)), file = "Cell_per_Cluster_E16_Epi_04.csv")
#List of Epithelial genes to find clusters
Epi_Genes <- c("Cnn1","Acta2",
               "Krt5","Krt14","Krt19",
               "Smgc","Bpifa2","Aqp5","Bhlha15",
               "Epcam","Cdh1")
Proliferation_Genes <- c("Cenpa","Cenpe","Cenpf","Mki67",
                         "Hist1h1a","Hist1h1b","Hist1h2ae","Top2a")
#Loop to create Feature plots and violin plots for Geneslists
for(i in Proliferation_Genes[1:8]){
  print(i)
  plot3 <- FeaturePlot(agg.E16.Epi, features = i, pt.size = 3, min.cutoff = "q9") + theme(title = element_text(size = 30), axis.title = element_text(size = 20), axis.text = element_text(size = 18))
  plot4 <- VlnPlot(agg.E16.Epi, features = i) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
  Dot_Plot_Title <- str_c("UMAP_Plot_Epi_", i, "_E16_04.pdf")
  Violin_Plot_Title <- str_c("Violin_Plot_Epi_", i, "_E16_04.pdf")
  pdf(file = Dot_Plot_Title,  width = 8, height = 5)
      print(plot3)
  dev.off()  
  pdf(file = Violin_Plot_Title, width = 8, height = 4.2)
      print(plot4)
  dev.off()  
}
#Draws a Dot plot for genes of interest for Epithelium
pdf("DotPlot_GeneIDs_Epi_ClusterNumber_E16_All_04.pdf", width = 7, height = 5)
    DotPlot(agg.E16.Epi, features = Epi_Genes, cols = "Spectral", col.min = 0, dot.scale = 5) +
            RotatedAxis()
dev.off()
#Finds all marker genes for every cluster compared to all other clusters
agg.E16.markers.All <- FindAllMarkers(agg.E16.Epi, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.All, file = "All_Gene_Averages_perEpiCluster_E16_04.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos <- FindAllMarkers(agg.E16.Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.Pos, file = "Positive_Gene_Averages_perEpiCluster_E16_04.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos.Top10 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.E16.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perEpiCluster_E16_04.csv")
#Draws a Dot plot for 5 Positive genes defining each stroma cluster
pdf("DotPlot_Epi_Top10genes_Labeled_E16_04.pdf", width = 20, height = 10)
    DotPlot(object = agg.E16.Epi, features = unique(agg.E16.markers.Pos.Top10$gene), cols = "Spectral", col.min = 0, dot.scale = 5) +
            RotatedAxis()
dev.off()

############################## Stromal subset analysis ######################################
################################################################################################

  #Re-calculations of Stroma compared to other stroma
#Creates a subset of agg.E16 based on "Stroma" label
agg.E16.Stroma <- subset(agg.E16.2, idents = c("Stroma"))
#Rerun to reidentifies the most highly variable features, and returns 2000 features per dataset as a default
agg.E16.Stroma <- FindVariableFeatures(agg.E16.Stroma, selection.method = "vst", nfeatures = 2000)
#Rerun to reapplies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.E16.Stroma <- ScaleData(agg.E16.Stroma, feature = all.genes) 

#generates list of Top 10 most variable genes
top10.Stroma <- head(VariableFeatures(agg.E16.Stroma), 10)

#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
pdf("Elbowplot_E16_Stroma_04.pdf")
  ElbowPlot(agg.E16.Stroma, ndims = 50)
dev.off()

#Reruns PCA statistics used for generates a PCA plot of the first PCs
agg.E16.Stroma <- RunPCA(agg.E16.Stroma, features = VariableFeatures(object = agg.E16.Stroma))
#Reruns for constructing a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset (first 40 PCs, dims = 40), how many dimentions do you want included for the analysis.
agg.E16.Stroma <- FindNeighbors(agg.E16.Stroma, dims = 1:40) 
#Reruns generation of the Clusters, resolution increases the specificity of each cluster, resolution = 2.0, the default is 1.0
agg.E16.Stroma <- FindClusters(agg.E16.Stroma, resolution = 0.8)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.E16.Stroma <- RunTSNE(agg.E16.Stroma, dims = 1:40)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.E16.Stroma <- RunUMAP(agg.E16.Stroma, dims = 1:40)

#Saves Formated dataset after all major computational processing
saveRDS(agg.E16.Stroma, file = "E16_SMG_Raw_Str_(SEURAT_v4)_04.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Stroma_only")
agg.E16.Stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Stroma_only/E16_SMG_Raw_Str_(SEURAT_v4)_04.rds")

pdf("UMAP_E16_Stroma_ClusterNumbers_04_2.pdf", width = 6.25, height = 5)
  DimPlot(agg.E16.Stroma, reduction = "umap", label = TRUE, pt.size = 1.5) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()

#New list with simple cluster names
new.cluster.ids <- c("Pa+Pb_00",
                     "Pa+Pb_01",
                     "Pa+Thy1_02",
                     "Pa+Pb_CC_03",
                     "Pa+Pb_CC_04",
                     "Pa+Pb_05",
                     "Pb_06",
                     "Pa+Pb_CC_07",
                     "Pb_CC_08",
                     "Pa_09",
                     "Pb+Acta2_10",
                     "Acta2_11",
                     "Pa+Pb_12"
                     )
#New list with GSEA enrichment labels
new.cluster.ids <- c("Vehicle_enriched",
                     "FGF2_enriched",
                     "Vehicle_enriched",
                     "FGF2_enriched",
                     "FGF2_enriched",
                     "FGF2_enriched",
                     "Vehicle_enriched",
                     "FGF2_enriched",
                     "FGF2_enriched",
                     "Vehicle_enriched",
                     "Vehicle_enriched",
                     "Vehicle_enriched",
                     "FGF2_enriched"
                      )
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Names_Stroma_04.csv")
#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.E16.Stroma)), file = "Cell_per_Cluster_E16_Stroma_04.csv")
#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.E16.Stroma)
agg.E16.Stroma <- RenameIdents(agg.E16.Stroma, new.cluster.ids)

#Generates a t-SNE plot with labels
pdf("TSNE_E16_Stroma_labels_04.pdf")
  DimPlot(agg.E16.Stroma, reduction = "tsne", label = TRUE)
dev.off() 
#Generates a UMAP plot with and without labels
pdf("UMAP_E16_Stroma_labels_04.pdf", width = 8, height = 5)
  DimPlot(agg.E16.Stroma, reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, pt.size = 1.5) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
dev.off()

#Finds all marker genes for every cluster compared to all other clusters
agg.E16.markers.All <- FindAllMarkers(agg.E16.Stroma, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.All, file = "All_Gene_Averages_perStromaCluster_E16_04.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.E16.markers.Pos <- FindAllMarkers(agg.E16.Stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.E16.markers.Pos, file = "Positive_Gene_Averages_perStromalClustersAlone_E16_04.csv")
#Reduces to top 20 postive marker genes for every cluster compared to all other clusters
agg.E16.Stroma.markers.Pos.Top20 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(agg.E16.Stroma.markers.Pos.Top20, file = "Positive_Top20_Stroma_Gene_Averages_perCluster_E16_04.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.E16.Stroma.markers.Pos.Top10 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.E16.Stroma.markers.Pos.Top10, file = "Positive_Top10_Stroma_Gene_Averages_perCluster_E16_04.csv")
#Reduces to top 5 postive marker genes for every cluster compared to all other clusters
agg.E16.Stroma.markers.Pos.Top5 <- agg.E16.markers.Pos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.csv(agg.E16.Stroma.markers.Pos.Top5, file = "Positive_Top5_Stroma_Gene_Averages_perCluster_E16_04.csv")
#Draws a Dot plot for 5 Positive genes defining each stroma cluster
pdf("DotPlot_Stroma_Top5genes_Labeled_E16_04.pdf", width = 20, height = 10)
    DotPlot(object = agg.E16.Stroma, features = unique(agg.E16.Stroma.markers.Pos.Top5$gene), cols = "Spectral", col.min = 0, dot.scale = 5) +
      RotatedAxis()
dev.off()
#Heatmap plot of Top 10 genes defining only the stromal clusters
pdf("Heatmap_Top20_Stroma_E16_04.pdf", width = 40, height = 20)
  DoHeatmap(agg.E16.Stroma, features = agg.E16.Stroma.markers.Pos.Top20$gene)
dev.off()
#Heatmap plot of Top 10 genes defining only the stromal clusters
pdf("Heatmap_Top10_Stroma_E16_04.pdf", width = 30, height = 15)
    DoHeatmap(agg.E16.Stroma, features = agg.E16.Stroma.markers.Pos.Top10$gene)
dev.off()
#Heatmap plot of Top 5 genes defining only the stromal clusters
pdf("Heatmap_Top5_Stroma_E16_04.pdf", width = 20, height = 10)
    DoHeatmap(agg.E16.Stroma, features = agg.E16.Stroma.markers.Pos.Top5$gene)
dev.off()

#Plot for Pdgfra in UMAP orientation.
pdf("UMAP_E16_Stroma_Labels_Pdgfra_04.pdf")
  FeaturePlot(agg.E16.Stroma, features = c("Pdgfra"), min.cutoff = "q9")
dev.off()

#loop to make top gene graphs for differential expression
# List of genes to run in loop
Genes_of_Interest <- c("Ptprc", "Pdpn", "Prox1", "Lyve1","Hba-a1", "Hba-a2",
                       "Hbb-bs", "Kdr", "Cdh5", "Pecam1", "Tie1", "Vcam1", "Tubb3",
                       "Fabp7", "Mbp", "Plp1", "Ptprz1", "Pdgfra", "Pdgfrb",
                       "Eng", "Vim", "Thy1", "Ly6a", "Nt5e", "Col1a1", "Gli1", "Fgf2",
                       "Fgf10", "Bmp2", "Bmp4", "Bmp7", "Kitl", "Tgfb1", "Tgfb2", "Tgfb3",
                       "Tgfbr1", "Tgfbr2", "Tgfbr3", "Ccn2", "Kitl", "Egf", "Areg", "Ereg",
                       "Wnt2", "Wnt4", "Col4a1", "Acta2", "Cnn1", "Il7r", "S100a4",
                       "Ms4a2", "Hdc", "Cenpa", "Cenpe", "Cenpf", "Mki67", "Hist1h1a", 
                       "Hist1h1b", "Top2a", "Hist1h2ae")
for(i in Genes_of_Interest[1:3]){
  print(i)
  plot3 <- FeaturePlot(agg.E16.Stroma, features = i, pt.size = 3, min.cutoff = "q9") + theme(title = element_text(size = 30), axis.title = element_text(size = 20), axis.text = element_text(size = 18))
  plot4 <- VlnPlot(agg.E16.Stroma, features = i) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
  Dot_Plot_Title <- str_c("UMAP_Plot_Stroma_", i, "_E16_04.pdf")
  Violin_Plot_Title <- str_c("Violin_Plot_Stroma_", i, "_E16_04.pdf")
  pdf(file = Dot_Plot_Title, width = 8, height = 5)
      print(plot3)
  dev.off()  
  pdf(file = Violin_Plot_Title, width = 8, height = 4)
      print(plot4)
  dev.off()  
}

Genes_of_Interest <- c("Acta2",
                       "Cnn1",
                       "Ccn2",
                       "Kdr",
                       "Cdh5",
                       "Vcam1",
                       "Pdgfra",
                       "Pdgfrb",
                       "Eng",
                       "Vim",
                       "Thy1",
                       "Gli1", 
                       "Fgf2",
                       "Fgf10",
                       "Bmp2",
                       "Bmp4",
                       "Bmp7",
                       "Kitl", 
                       "Tgfb1",
                       "Tgfb2",
                       "Tgfb3",
                       "Col1a1",
                       "Col4a1"
                       ) 
#Draws a Dot plot for genes of interest
pdf("DotPlot_Stroma_GeneIDs_Labeled_E16_04.pdf", width = 18, height = 10)
  DotPlot(agg.E16.Stroma, features = Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
pdf("DotPlot_Stroma_PDGFRa_FGF2_Labeled_E16_04.pdf", width = 7, height = 5)
  DotPlot(agg.E16.Stroma, features = c("Thy1","Acta2","Pdgfrb","Pdgfra","Fgf2","Fgf10"), cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
dev.off()
pdf("DotPlot_Stroma_PDGFRa_BMPs_EGFs_Labeled_E16_04.pdf", width = 7, height = 5)
  DotPlot(agg.E16.Stroma, features = c("Thy1","Acta2","Pdgfrb","Pdgfra","Bmp2","Bmp4","Bmp7","Areg","Ereg","Egf"), cols = "Spectral", col.min = 0, dot.scale = 5) +
    RotatedAxis()
dev.off()
#Make double feature plots for markers that separate stromal clusters
pdf("UMAP_E16_Stroma_Thy1_Pdgfra_04.pdf", width = 40, height = 8)
    FeaturePlot(agg.E16.Stroma, features = c("Thy1", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Magenta", "Green"))
dev.off()
pdf("UMAP_E16_Stroma_Pdgfrb_Pdgfra_04.pdf", width = 40, height = 8)
    FeaturePlot(agg.E16.Stroma, features = c("Pdgfrb", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","Green"))
dev.off()
pdf("UMAP_E16_Stroma_Pdgfrb_aSMA_04.pdf", width = 40, height = 8)
    FeaturePlot(agg.E16.Stroma, features = c("Pdgfrb", "Acta2"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
dev.off()

  #First method for counting cells
# Function that can calculate proportion of cells expressing a gene
# https://github.com/satijalab/seurat/issues/371
# updated 1/31/2020 to accommodate V3.1
# updated 2/4/2020 to output "NA" for genes not detected in certain subgroups
# calculates total cells expressing a gene (raw counts > 0) by metadata groups
# can be grouped by different samples types or cluster_# based on metadata
# 'ncells' counts to total number of cells, can be passed to have percentages in calc_helper
# you can adjust the threshold for RNA count to select for cells with more 'higher' expression
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}
# for total cells:
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)
  }else{return(NA)}
}
# for percentage of cells use this function:
#calc_helper <- function(object,genes){
#  counts = object[['RNA']]@counts
#  ncells = ncol(counts)
#  if(genes %in% row.names(counts)){
#    sum(counts[genes,]>0)/ncells
#  }else{return(NA)}
#}

# Finds and saves a .csv of cell numbers expressing certain genes
#performs the counts based on the cluster assignment
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04")
agg.E16 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/E16_SMG_Raw_(SEURAT_v4)_04.rds")
NumberCellperGene <- PrctCellExpringGene(agg.E16,genes = c("Pdgfra","Pdgfrb"), group.by = "seurat_clusters")
write.csv(NumberCellperGene, file = "NumberOfCellsPerGene_Pdgfra_Pdgfrb_04.csv")

  #Second method for counting cells, can do co-positive cells
#https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FetchData
#This way for scripting
#FetchData() can pull out cell, gene expression data, and any other data from Seurat objects
CellData <- FetchData(agg.E16, vars = c("Pdgfra","Pdgfrb","Gli1","Acta2"), slot = "data")
#Calculates the sum of cells that express gene1
Pdgfra_Pos <- sum(CellData$Pdgfra>0)
#Calculates the sum of cells that express gene2
Pdgfrb_Pos <- sum(CellData$Pdgfrb>0)
#Calculates the sum of cells that express gene3
Gli1_Pos <- sum(CellData$Gli1>0)
#Calculates the sum of cells that express gene4
Acta2_Pos <- sum(CellData$Acta2>0)
#Subsets the data for Postive for gene1 and then positive for gene2 then positive for gene3
Co_Postitive_cell_ID <- CellData %>% subset(Pdgfra>0) %>% subset(Pdgfrb>0) %>% subset(Acta2>0)  
#Number of Co_positive cells
Co_Postitive_cell_Total <- nrow(Co_Postitive_cell_ID)

#This is a function that does the same thing
# for total co-positive cells:
Co_Positive <- function(object,gene1,gene2){
  CellData <- FetchData(object, vars = c(gene1,gene2), slot = "data")
  print(head(CellData))
  Gene1_total_cells <- sum(CellData[,gene1]>0)
  cat("Total cells positive for", gene1,Gene1_total_cells," ")
  Gene2_total_cells <- sum(CellData[,gene2]>0)
  cat("Total cells positive for", gene2,Gene2_total_cells)
  Co_Postitive_cell_ID <- subset(CellData, CellData[,gene1]>0)
  Co_Postitive_cell_ID <- subset(Co_Postitive_cell_ID, Co_Postitive_cell_ID[,gene2]>0)
  print(head(Co_Postitive_cell_ID))
  Co_Postitive_cell_Total <- nrow(Co_Postitive_cell_ID)
  cat("Total cells positive for",gene1," & ",gene2," ", Co_Postitive_cell_Total)
}
