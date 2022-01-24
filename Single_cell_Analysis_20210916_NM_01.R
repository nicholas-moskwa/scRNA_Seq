#Version 1
#This version is using data that is deconvoluted
#Major parameter and statistical decisions made:
#noFGF2_organoids
	#min.cells = 3, min.features = 200
	#subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & 
  #noFGF2_Organoid percent.mt < 5%
  #noFGF2_Organoid = 35 dimensions
  #resolution = 1.0
#Extracts UMAP data and graphs with standard sizing
#Labels both broad and specific added to graphs
#Has code to subset out stromal data
#Analyzes data, using violin and UMAP plots
#Creates gene expression lists based on clusters
#Integrates noFGF2 organoid stroma with E16 SMG stroma
#Includes function for finding number/percentage of cells expressing a specific gene


######## sessionInfo()
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
#Sets the working directory where the single cell data
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01")

############################## Functions #######################################################
################################################################################################


#Funtion "Make_pdf" makes a pdf file for plot with "Title", quick a dirty image generator/save
Make_pdf <- function(Title, plot){
                                  pdf(Title)
                                    print(plot)
                                  dev.off()
                                  }

############################## Format and compute datpset ######################################
################################################################################################


#Reads the 10x Generated data that is in matrix form; barcodes, features, matrix, into a large data matrix
agg.data.noFGF2 <- Read10X(data.dir =
                        "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/ML-2/outs/filtered_feature_bc_matrix")
#Creates Seurat object using the 10x raw data matrix
agg.data.noFGF2 <- CreateSeuratObject(agg.data.noFGF2, project = "agg.data.noFGF2", min.cells = 3,
                              min.features = 200)
#Creates new parameter selecting for mitochondrial genes
agg.data.noFGF2[["percent.mt"]] <- PercentageFeatureSet(agg.data.noFGF2, pattern = "^mt-")
#Generates Volcano Plots to show parameters of dataset
#Useful QC measurements to establish data quality and cutoffs
plot1 <- VlnPlot(agg.data.noFGF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Make_pdf("Volcano_plot_Feature_RNAcount_Mitochondira_noFGF2_01.pdf", plot1)

#Saves a scatterplot of Features vs mitochondrial genes
plot1 <- FeatureScatter(agg.data.noFGF2, feature1 = "nCount_RNA", feature2 = "percent.mt")
Make_pdf("Feature_vs_Mt_gene_ScatternPlot_noFGF2_01.pdf", plot1)

#Saves a scatterplot of Features vs Counts
plot1 <- FeatureScatter(agg.data.noFGF2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
Make_pdf("Feature_vs_RNAcount_ScatterPlot_noFGF2_01.pdf", plot1)

#Filter out cells that have unique feature counts more than 200 and less than 9000
#idea is to remove fake cells that are either, free RNA in a droplet, <200
#or doublets/triplets, two cells in a droplet, >9000
#noFGF2 organoid Mitochondrial cutoff for near gaussian distribution = 5%
#More mitochondrial RNA reflects apoptotic stress that is mitochondrially driven
agg.data.noFGF2 <- subset(agg.data.noFGF2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

plot1 <- VlnPlot(agg.data.noFGF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Make_pdf("Volcano_plot_Feature_RNAcount_Mitochondira_Truncation_noFGF2_01.pdf", plot1)

#Normalizes dataset
agg.data.noFGF2 <- NormalizeData(agg.data.noFGF2)

#Identifies the most highly variable features, and returns 2000 features per dataset as a default
agg.data.noFGF2 <- FindVariableFeatures(agg.data.noFGF2, selection.method = "vst", nfeatures = 2000)

#generates list of Top 10 most variable genes
top10.noFGF2 <- head(VariableFeatures(agg.data.noFGF2), 10)

#Saves plot of selected features
plot1 <- VariableFeaturePlot(agg.data.noFGF2)

#Saves plot of selected features and displays top10 most variable genes
plot1 <- LabelPoints(plot = plot1, points = top10.noFGF2, repel = TRUE)
Make_pdf("Varianceplot_top10_noFGF2_01.pdf", plot1)

#Labels all genes
all.genes.noFGF2 <- rownames(agg.data.noFGF2)
#Applies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.data.noFGF2 <- ScaleData(agg.data.noFGF2, feature = all.genes.noFGF2)
#Runs PCA statistics used for generates a PCA plot of the first PCs
agg.data.noFGF2 <- RunPCA(agg.data.noFGF2, features = VariableFeatures(object = agg.data.noFGF2))

#Generates a heatmap showing which genes are contributing to the hetergeneity that determines PC1, dims = "PC number", shows which genes are contributing most to the cell seperation.
#SOMETHING IS WRONG WITH THIS "Dimheatmap" FUNCTION, can't be saved in RAM, be careful about physical memory needed to graph these, you might need to just save it without viewing in R.
plot1 <- DimHeatmap(agg.data.noFGF2, dims = 1, cells = 500, balanced = TRUE)
Make_pdf("PCA1_Heatmap_noFGF2_01.pdf", plot1)
plot1 <- DimHeatmap(agg.data.noFGF2, dims = 2, cells = 500, balanced = TRUE)
Make_pdf("PCA2_Heatmap_noFGF2_01.pdf", plot1)
plot1 <- DimHeatmap(agg.data.noFGF2, dims = 3, cells = 500, balanced = TRUE)
Make_pdf("PCA3_Heatmap_noFGF2_01.pdf", plot1)

#Saves PCA plot total nomalized dataset based on PC 1 and PC 2
plot1 <- DimPlot(agg.data.noFGF2, reduction = "pca")
Make_pdf("PCA_1v2_noFGF2_01.pdf", plot1)
#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
#Around 30 dimensions is where adding another dimentions adds little more extra separation by the Principle Components

plot1 <- ElbowPlot(agg.data.noFGF2, ndims = 50)
Make_pdf("Elbowplot_noFGF2_01.pdf", plot1)

#constructs a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset, how many dimentions do you want included for the analysis.
#noFGF2 Oragnoid = 30 dimensions
agg.data.noFGF2 <- FindNeighbors(agg.data.noFGF2, dims = 1:35) 
#Generates the Clusters, resolution increases the specificity of each cluster, resolution = 2.0, the default is 1.0
agg.data.noFGF2 <- FindClusters(agg.data.noFGF2, resolution = 1.0)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.data.noFGF2 <- RunTSNE(agg.data.noFGF2, dims = 1:35)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.data.noFGF2 <- RunUMAP(agg.data.noFGF2, dims = 1:35)

#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.data.noFGF2)), file = "Cell_per_Cluster_Organoid_noFGF2_01.csv")

#Generates a t-SNE plot
plot1 <-  DimPlot(agg.data.noFGF2, reduction = "tsne", label = TRUE)
Make_pdf("TSNE_noFGF2_01.pdf", plot1)

#Generates a UMAP plot
plot1 <-  DimPlot(agg.data.noFGF2, reduction = "umap", label = TRUE)
Make_pdf("UMAP_noFGF2_01.pdf", plot1)

#Saves Formated dataset after all major computational processing; Total data from no FGF2 organoids
saveRDS(agg.data.noFGF2, file = "Organoid_noFGF2_Raw_(SEURAT_v4)_01.rds")

#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01")
agg.data.noFGF2 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01/Organoid_noFGF2_Raw_(SEURAT_v4)_01.rds")

############################## Subsets Stroma only data ########################################
################################################################################################

#Creates a subset of agg.data.noFGF2 based analysis
#This is based of UMAP and gene data that is established in next section
agg.noFGF2.Stroma <- subset(agg.data.noFGF2, idents = c("0","1","2","3","4","5","6","7","8"))
#Saves Formated dataset after all major computational processing; Stroma only data
saveRDS(agg.noFGF2.Stroma, file = "Organoid_noFGF2_Stroma_(SEURAT_v3)_01.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01")
agg.noFGF2.Stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01/Organoid_noFGF2_Stroma_(SEURAT_v3)_01.rds")

    #All UMAP/Feature gene plot maps used for assigning clusters labels
	#NOTE at the end of these there is a coded For-loop that will do this automatically for any gene lists.
  #LYMPHATICS ENDOTHELIUM
#Shows individual cells on cluster plot, Leucocytes = PTPRC
# "Ptprc", "Pdpn", "Prox1", "Lyve1"
  #RED BLOOD CELLS
#Shows individual cells on cluster plot, Erythrocytes = Hemoglobins
# "Hba-a1", "Hba-a2", "Hbb-bs"
  #ENDOTHELIUM
# Endothelial = "Kdr", "Cdh5", "Pecam1", "Tie1", "Vcam1"
  #NEURONS
#Neurons = Tubb3
  #SCHWANN cells / Oligodendrocytes
# "Fabp7", "Mbp", "Plp1", "Ptprz1", "Sox10"
  #STROMAL genes
# Pdgfra
# Pdgfrb
# Eng or "CD105"
# Vim Vimentin
# Thy1 or "CD90"
# Ly6a or "Sca-1"
# Nt5e
# Col1a1
# Gli1
# Fgf2
# Fgf10
# Bmp4
# Bmp7
# Kit Ligand
# Tgfb1
# Tgfb2
# Tgfbr3
# Wnt2
# Wnt4
# Col4a1
  #EPITHELIUM
#Acinar Epi = Aqp5
#Proacinar Epi = Bhlha15 or "Mist1"
#Proacinar Epi = Bpifa2 or "PSP"
#Proacinar Epi = Muc19
#Proacinar Epi = Smgc
#Epithelial = Cdh1
#Epithelial = Epcam
#Epithelial = Krt5
#Epithelial = Krt14
#Epithelial = Krt19
#Epithelial = Smooth muscle actin, Acta2
#Epithelial = Calponin 1, Cnn1
#Epithelium = Sox10
#Epithelium = Sox2
  #IMMUNE
#Immune Memory CD4+ = Il7r
#Immune Memory CD4+ = S100A4
#Immune Mast cell = Ms4a2
#Immune Mast cell = Hdc
  #CELL CYCLE / Mitosis
#Mitosis G2 = "Cenpa", "Cenpe", "Cenpf", "Mki67"
#Mitosis G1/S = "Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae"

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
Stromal_Markers <- c("Thy1", "Cnn1", "Acta2","Pdgfrb", "Pdgfra", "Vim", "Col1a1")
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
#Loop to create Feature plots and violin plots for Geneslists, can pass any gene-list to generate plots
for(i in Genes_of_Interest[1:67]){
  print(i)
  plot3 <- FeaturePlot(agg.data.noFGF2, features = i, min.cutoff = "q9")
  plot4 <- VlnPlot(agg.data.noFGF2, features = i) + NoLegend()
  Dot_Plot_Title <- str_c("UMAP_Plot_", i, "_noFGF2_01.pdf")
  Violin_Plot_Title <- str_c("Violin_Plot_", i, "_noFGF2_01.pdf")
  pdf(file = Dot_Plot_Title)
  print(plot3)
  dev.off()  
  pdf(file = Violin_Plot_Title)
  print(plot4)
  dev.off()  
}

#Draws a Dot plot for genes of interest
pdf("DotPlot_GeneIDs_Broad_Labeled_E16_01.pdf", width = 8, height = 4)
plot1 <- DotPlot(agg.data.noFGF2, features = Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
                  RotatedAxis()
dev.off()


  #Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.noFGF2.markers.All <- FindAllMarkers(agg.data.noFGF2, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.noFGF2.markers.All, file = "All_Gene_Averages_perCluster_noFGF2_01.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.noFGF2.markers.Pos <- FindAllMarkers(agg.data.noFGF2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.noFGF2.markers.Pos, file = "Positive_Gene_Averages_perCluster_noFGF2_01.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.noFGF2.markers.Pos.Top10 <- agg.noFGF2.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.noFGF2.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perCluster_noFGF2_01.csv")

    #Naming clusters
#Old list with cluster numbers
old.cluster.ids <- c("0", "1", "2", "3", "4",
                      "5", "6", "7", "8", "9",
                      "10", "11", "12", "13", "14"
                      )
  #Making Broad Cellular labels
new.cluster.ids <- c(
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Stroma",
  "Epithelium",
  "Neuron",
  "Immune",
  "Neuron",
  "Epithelium",
  "Endothelium"
)
  #New list with specific cluster names
new.cluster.ids <- c(
  "Str_Acta2_Pdgfra_Pdgfrb_00",
  "Str_Pdgfra_Pdgfrb_Thy1_01",
  "Str_Pdgfra_Pdgfrb_02",
  "Str_Acta2_03",
  "Str_Acta2_04",
  "Str_Acta2_Cnn1_05",
  "Str_Pdgfra_Pdgfrb_06",
  "Str_Acta2_CC_07",
  "Str_Acta2_CC_08",
  "Epi_Duct_09",
  "Neuron_10",
  "Immune_11",
  "Neuron_12",
  "Epi_Duct_13",
  "Endothelium_14"
  )

#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.data.noFGF2)
agg.data.noFGF2 <- RenameIdents(agg.data.noFGF2, new.cluster.ids)
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Broad_Names_Organoids_noFGF2_01.csv")
write.csv(new.cluster.ids, file = "Cluster#_Specific_Names_Organoids_noFGF2_01.csv")
#Saves/catalogs cluster names to cell numbers
write.csv(table(Idents(agg.data.noFGF2)), file = "Cell_per_Cluster_Broad_Labels_Organoids_noFGF2_01.csv")
write.csv(table(Idents(agg.data.noFGF2)), file = "Cell_per_Cluster_Specific_Labels_Organoids_noFGF2_01.csv")
#Draws a new UMAP plot with written labels from new.cluster.ids
plot1 <- DimPlot(agg.data.noFGF2, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5, repel = TRUE) +
                NoLegend() +
                theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("UMAP_Organoids_noFGF2_Labeled_Broad_01.pdf", plot1)
Make_pdf("UMAP_Organoids_noFGF2_Labeled_Specific_01.pdf", plot1)

#Genes that show up in broad cell types
Broad_Label_Genes  <- c("Col1a1",
                        "Acta2",
                        "Epcam",
                        "Sox9",
                        "Gap43",
                        "Gas7",
                        "Ptprc",
                        "Ccl6",
                        "Kdr",
                        "Cdh5"
                        )
#Draws a Dot plot for genes of interest for Broad labels
pdf("DotPlot_GeneIDs_Broad_Labeled_E16_01.pdf", width = 8, height = 4)
plot1 <- DotPlot(agg.data.noFGF2, features = Broad_Label_Genes, cols = "Spectral", col.min = 0, dot.scale = 5) +
        RotatedAxis()
plot1
dev.off()
