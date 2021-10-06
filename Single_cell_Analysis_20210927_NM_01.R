#Version 1
#Uses Seurat packages for analysis: https://satijalab.org/seurat/index.html
#This version is using data that is deconvoluted using CellRanger
#Includes function for finding number/percentage of cells expressing a specific gene
#Major parameter and statistical decisions made:
 #yesFGF2_Organoid
	#min.cells = 3, min.features = 200
	#subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & 
  #percent.mt < 5%
  #39 dimensions
  #resolution = 1.0
#Extracts UMAP data and graphs with standard sizing
#Labels both broad and specific added to graphs
#Subsets out stromal data and re-analyzes the data
#Analyzes data, using violin and UMAP plots for lists of genes
#Creates gene expression lists based on clustering
#Integrates yesFGF2 organoid stroma, noFGF2 organoid stroma with E16 SMG stroma
 #noFGF2_Organoid
  #min.cells = 3, min.features = 200
  #subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & 
  #percent.mt < 5%
  #35 dimensions
  #resolution = 1.0
 #E16_InVivo
  #min.cells = 3, min.features = 200
  #subset = nFeature_RNA > 200 & nFeature_RNA < 9000 &
  #percent.mt < 7%
  #40 dimensions
  #resolution = 2.0
 #Integration uses
  #20 dimensions
  #resolution = 0.6
############################## Libraries and start up ##########################################
################################################################################################

#Allows for python UMAP calculations using miniconda environment "seurat"
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
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01")

############################## Custom Functions ###################################################################
###################################################################################################################

#Funtion "Make_pdf" makes a pdf file for plot with "Title", quick a dirty image generator/save
#Title is Title, plot is the plot defined as an object, width and height are the .pdf's dimensions
Make_pdf <- function(Title, plot, width, height){
    pdf(Title, width = width, height = height)
    print(plot)
    dev.off()
}
								  
#Fucntion for quickly making UMAP and violin plots for all genes from a gene list
#Loop to create Feature plots and violin plots for GenesLists, can pass any gene-list to generate plots
Gene_Graphs <- function(GeneList, Object, Name){
	for(i in GeneList[1:100]){
	print(i)
	plot3 <- FeaturePlot(Object, features = i, min.cutoff = "q9")
	plot4 <- VlnPlot(Object, features = i, pt.size = 0.1) + ylim(0,NA) + NoLegend()
	Dot_Plot_Title <- str_c("UMAP_Plot_", i, "_", Name, "_01.pdf")
	Violin_Plot_Title <- str_c("Violin_Plot_", i, "_", Name, "_01.pdf")
	pdf(file = Dot_Plot_Title)
	print(plot3)
	dev.off()  
	pdf(file = Violin_Plot_Title)
	print(plot4)
	dev.off()  
	}
}

# Function that can calculate proportion of cells expressing a gene
# https://github.com/satijalab/seurat/issues/371
# updated 1/31/2020 to accommodate V3.1
# updated 2/4/2020 to output "NA" for genes not detected in certain subgroups
# calculates total cells expressing a gene (raw counts > 0) by metadata groups
# can be grouped by different samples types or cluster_# based on metadata
# 'ncells' counts to total number of cells, can be passed to have percentages in calc_helper
# you can adjust the threshold for RNA count to select for cells with more 'higher' expression
TotalCellExpringGene <- function(object, genes, group.by = "all"){
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

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper.2, object=object))
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
# for percentage of cells use this function:
calc_helper.2 <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
 #Example for finding all cells expressing a gene:
#Finds and saves a .csv of cell numbers expressing certain genes
#performs the counts based on the cluster assignment
NumberCellperGene <- PrctCellExpringGene(agg.E16.combined,genes = c("Vim","Pdgfra","Pdgfrb","Acta2","Cnn1"), group.by = "seurat_clusters")
write.csv(NumberCellperGene, file = "NumberOfCellsPerGene_Vim_Pdgfra_Pdgfrb_Acta2_Cnn1_Integrated_01.csv")
#performs the counts based on the sample assignment
NumberCellperGene <- PrctCellExpringGene(agg.E16.combined,genes = c("Vim","Pdgfra","Pdgfrb","Acta2","Cnn1"), group.by = "protocol")
write.csv(NumberCellperGene, file = "NumberOfCellsPerGenePerSample_Vim_Pdgfra_Pdgfrb_Acta2_Cnn1_Integrated_01.csv")

  #Second method for counting cells, can do co-positive cells
 #Example:
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

#This is a function that does the same thing as finding nubmer of cells expressing a gene, but for two genes.
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

############################## Data Analysis yesFGF2 Sample #######################################################
###################################################################################################################

#Reads the 10x Generated data that is in matrix form; barcodes, features, matrix, into a large data matrix
agg.data.yesFGF2 <- Read10X(data.dir =
                        "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/235FGF2-7/outs/filtered_feature_bc_matrix")
#Creates Seurat object using the 10x raw data matrix
agg.data.yesFGF2 <- CreateSeuratObject(agg.data.yesFGF2, project = "agg.data.yesFGF2", min.cells = 3,
                              min.features = 200)
#Creates new parameter selecting for mitochondrial genes
agg.data.yesFGF2[["percent.mt"]] <- PercentageFeatureSet(agg.data.yesFGF2, pattern = "^mt-")
#Generates Volcano Plots to show parameters of dataset
#Useful QC measurements to establish data quality and cutoffs
plot1 <- VlnPlot(agg.data.yesFGF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Make_pdf("Volcano_plot_Feature_RNAcount_Mitochondira_yesFGF2_01.pdf", plot1, 4, 3)

#Saves a scatterplot of Features vs mitochondrial genes
plot1 <- FeatureScatter(agg.data.yesFGF2, feature1 = "nCount_RNA", feature2 = "percent.mt")
Make_pdf("Feature_vs_Mt_gene_ScatternPlot_yesFGF2_01.pdf", plot1, 4, 3)

#Saves a scatterplot of Features vs Counts
plot1 <- FeatureScatter(agg.data.yesFGF2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
Make_pdf("Feature_vs_RNAcount_ScatterPlot_yesFGF2_01.pdf", plot1, 4, 3)

#Filter out cells that have unique feature counts more than 200 and less than 9000
#idea is to remove fake cells that are either, free RNA in a droplet, <200
#or doublets/triplets, two cells in a droplet, >9000
#noFGF2 organoid Mitochondrial cutoff for near gaussian distribution = 5%
#More mitochondrial RNA reflects apoptotic stress that is mitochondrially driven
agg.data.yesFGF2 <- subset(agg.data.yesFGF2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

plot1 <- VlnPlot(agg.data.yesFGF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Make_pdf("Volcano_plot_Feature_RNAcount_Mitochondira_Truncation_yesFGF2_01.pdf", plot1, 4, 3)

#Normalizes dataset
#Normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
#Normalized values are stored in object[["RNA"]]@data
agg.data.yesFGF2 <- NormalizeData(agg.data.yesFGF2)

#Identifies the most highly variable features, and returns 2000 features per dataset as a default
agg.data.yesFGF2 <- FindVariableFeatures(agg.data.yesFGF2, selection.method = "vst", nfeatures = 2000)

#generates list of Top 10 most variable genes
top10.yesFGF2 <- head(VariableFeatures(agg.data.yesFGF2), 10)

#Saves plot of selected features
plot1 <- VariableFeaturePlot(agg.data.yesFGF2)

#Saves plot of selected features and displays top10 most variable genes
plot1 <- LabelPoints(plot = plot1, points = top10.yesFGF2, repel = TRUE)
Make_pdf("Varianceplot_top10_yesFGF2_01.pdf", plot1, 4, 3)

#Labels all genes
all.genes.yesFGF2 <- rownames(agg.data.yesFGF2)
#Applies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.data.yesFGF2 <- ScaleData(agg.data.yesFGF2, feature = all.genes.yesFGF2)
#Runs PCA statistics used for generates a PCA plot of the first PCs
agg.data.yesFGF2 <- RunPCA(agg.data.yesFGF2, features = VariableFeatures(object = agg.data.yesFGF2))

#Generates a heatmap showing which genes are contributing to the hetergeneity that determines PC1, dims = "PC number", shows which genes are contributing most to the cell seperation.
#SOMETHING IS WRONG WITH THIS "Dimheatmap" FUNCTION, can't be saved in RAM, be careful about physical memory needed to graph these, you might need to just save it without viewing in R.
plot1 <- DimHeatmap(agg.data.yesFGF2, dims = 1, cells = 500, balanced = TRUE)
Make_pdf("PCA1_Heatmap_yesFGF2_01.pdf", plot1, 4, 3)
plot1 <- DimHeatmap(agg.data.yesFGF2, dims = 2, cells = 500, balanced = TRUE)
Make_pdf("PCA2_Heatmap_yesFGF2_01.pdf", plot1, 4, 3)
plot1 <- DimHeatmap(agg.data.yesFGF2, dims = 3, cells = 500, balanced = TRUE)
Make_pdf("PCA3_Heatmap_yesFGF2_01.pdf", plot1, 4, 3)

#Saves PCA plot total nomalized dataset based on PC 1 and PC 2
plot1 <- DimPlot(agg.data.yesFGF2, reduction = "pca")
Make_pdf("PCA_1v2_yesFGF2_01.pdf", plot1, 4, 3)
#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
#Around 30 dimensions is where adding another dimentions adds little more extra separation by the Principle Components

plot1 <- ElbowPlot(agg.data.yesFGF2, ndims = 50)
Make_pdf("Elbowplot_yesFGF2_01.pdf", plot1, 4, 3)

#constructs a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset, how many dimentions do you want included for the analysis.
#noFGF2 Oragnoid = 30 dimensions
agg.data.yesFGF2 <- FindNeighbors(agg.data.yesFGF2, dims = 1:39) 
#Generates the Clusters, resolution increases the specificity of each cluster, resolution = 1.0, the default is 1.0
agg.data.yesFGF2 <- FindClusters(agg.data.yesFGF2, resolution = 1.0)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.data.yesFGF2 <- RunTSNE(agg.data.yesFGF2, dims = 1:39)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.data.yesFGF2 <- RunUMAP(agg.data.yesFGF2, dims = 1:39)

#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.data.yesFGF2)), file = "Cell_per_Cluster_Organoid_yesFGF2_01.csv")

#Generates a t-SNE plot
plot1 <-  DimPlot(agg.data.yesFGF2, reduction = "tsne", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
          NoLegend() +
          theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15))
Make_pdf("TSNE_yesFGF2_01.pdf", plot1, 4, 3)

#Generates a UMAP plot
plot1 <-  DimPlot(agg.data.yesFGF2, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
          NoLegend() +
          theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15))
Make_pdf("UMAP_yesFGF2_01.pdf", plot1, 4, 3)

#Saves Formated dataset after all major computational processing; Total data from no FGF2 organoids
saveRDS(agg.data.yesFGF2, file = "Organoid_yesFGF2_Raw_(SEURAT_v3)_01.rds")

#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01")
agg.data.yesFGF2 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Organoid_yesFGF2_Raw_(SEURAT_v3)_01.rds")

############################## Gene lists and Cluster Identification ##############################################
###################################################################################################################

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
#Creates Feature plots and violin plots for Geneslists, can pass any gene-list to generate plots
Gene_Graphs(Genes_of_Interest, agg.data.yesFGF2, "yesFGF2")

#Draws a Dot plot for genes of interest
pdf("DotPlot_GeneIDs_Broad_NoLabels_yesFGF2_01.pdf", width = 16, height = 4)
  plot1 <- DotPlot(agg.data.yesFGF2, features = Genes_of_Interest, cols = "Spectral", col.min = 0, dot.scale = 5) +
                    RotatedAxis()
  plot1
dev.off()

  #Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.yesFGF2.markers.All <- FindAllMarkers(agg.data.yesFGF2, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.yesFGF2.markers.All, file = "All_Gene_Averages_perCluster_yesFGF2_01.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.yesFGF2.markers.Pos <- FindAllMarkers(agg.data.yesFGF2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Alternative method for finding all Postive marker genes, saves time.
#agg.yesFGF2.markers.Pos.2 <- subset(agg.yesFGF2.markers.All, subset = avg_logFC > 0)
write.csv(agg.yesFGF2.markers.Pos, file = "Positive_Gene_Averages_perCluster_yesFGF2_01.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.yesFGF2.markers.Pos.Top10 <- agg.yesFGF2.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.yesFGF2.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perCluster_yesFGF2_01.csv")

    #Naming clusters
#Old list with cluster numbers
old.cluster.ids <- c("0", "1", "2", "3", "4",
                      "5", "6", "7", "8", "9",
                      "10", "11", "12", "13", "14",
                      "15", "16"
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
  "Stroma", #Cell Cycling
  "Stroma",
  "Epithelium",
  "Stroma",
  "Stroma", #True smooth muscle?
  "Immune", #Macrophage
  "Immune", #Mast Cell
  "Nerve"
)
  #New list with specific cluster names
new.cluster.ids <- c(
  "Str_Pdgfra_Thy1_00",
  "Str_Pdgfra_Thy1_01",
  "Str_Pdgfra_Pdgfrb_Thy1_02",
  "Str_Pdgfra_Thy1_03",
  "Str_Pdgfra_Thy1_Acta2_04",
  "Str_Pdgfra_Pdgfrb_Thy1_05",
  "Str_Pdgfra_Thy1_06",
  "Str_Pdgfra_Pdgfrb_07",
  "Str_Pdgfra_Pdgfrb_08",
  "Str_Pdgfra_Pdgfrb_Thy1_Acta2_CC_09",
  "Str_Pdgfra_Pdgfrb_10",
  "Epi_11",
  "Str_Pdgfra_Pdgfrb_12",
  "Str_Acta2_13",
  "Macrophage_14",
  "Mast_Cell_15",
  "Neuron_16"
)

#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.data.yesFGF2)
agg.data.yesFGF2 <- RenameIdents(agg.data.yesFGF2, new.cluster.ids)
#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Broad_Names_Organoids_yesFGF2_01.csv")
write.csv(new.cluster.ids, file = "Cluster#_Specific_Names_Organoids_yesFGF2_01.csv")
#Saves/catalogs cluster names to cell numbers
write.csv(table(Idents(agg.data.yesFGF2)), file = "Cell_per_Cluster_Broad_Labels_Organoids_yesFGF2_01.csv")
write.csv(table(Idents(agg.data.yesFGF2)), file = "Cell_per_Cluster_Specific_Labels_Organoids_yesFGF2_01.csv")
#Draws a new UMAP plot with written labels from new.cluster.ids
plot1 <- DimPlot(agg.data.yesFGF2, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5, repel = TRUE) +
                NoLegend() +
                theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("UMAP_Organoids_yesFGF2_Labeled_Broad_01.pdf", plot1, 4, 3)
Make_pdf("UMAP_Organoids_yesFGF2_Labeled_Specific_01.pdf", plot1, 10, 6)

#Genes that show up in broad cell types
Broad_Label_Genes  <- c("Col1a1",
                        "Vim",
                        "Epcam",
                        "Sox9",
                        "Wfdc17",
                        "Lyz2",
                        "Mcpt4",
                        "Ms4a2",
                        "Kcna1",
                        "Tubb3"
                        )
#Draws a Dot plot for genes of interest for Broad labels
plot1 <- DotPlot(agg.data.yesFGF2, features = Broad_Label_Genes, cols = "Spectral", col.min = 0, dot.scale = 5) +
          RotatedAxis()
Make_pdf("DotPlot_GeneIDs_Broad_Labeled_yesFGF2_01.pdf", plot1, 6, 3)

############################## Stroma data subset ##############################################
################################################################################################

#Subsets Stroma only data
#Creates a subset of agg.data.yesFGF2 based analysis
#This is based of UMAP and gene data that is established in previous section
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Stroma_Only")
agg.yesFGF2.Stroma <- subset(agg.data.yesFGF2, idents = c("0","1","2","3","4","5","6","7","8","9","10","12","13"))
#Saves Formated dataset after all major computational processing; Stroma only data
saveRDS(agg.yesFGF2.Stroma, file = "Organoid_yesFGF2_Stroma_(SEURAT_v3)_01.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Stroma_Only")
agg.yesFGF2.Stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Organoid_yesFGF2_Stroma_(SEURAT_v3)_01.rds")

	#Stormal subset analysis with re-calculations of Stroma compared to other stroma
  #This method is used in the paper
#Rerun to reidentifies the most highly variable features, and returns 2000 features per dataset as a default
agg.yesFGF2.Stroma <- FindVariableFeatures(agg.yesFGF2.Stroma, selection.method = "vst", nfeatures = 2000)
#Labels all genes
all.genes.yesFGF2.Stroma <- rownames(agg.yesFGF2.Stroma)
#Rerun to reapplies a linear transformation (â scalingâ ) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
agg.yesFGF2.Stroma <- ScaleData(agg.yesFGF2.Stroma, feature = all.genes.yesFGF2.Stroma) 

#generates list of Top 10 most variable genes
top10.Stroma <- head(VariableFeatures(agg.yesFGF2.Stroma), 10)

#Generates an elbow plot:a ranking of principle components based on the percentage of variance explained by each one.
plot1 <- ElbowPlot(agg.yesFGF2.Stroma, ndims = 50)
Make_pdf("Elbowplot_yesFGF2_Stroma_01.pdf",plot1,4,3)

#Reruns PCA statistics used for generates a PCA plot of the first PCs
agg.yesFGF2.Stroma <- RunPCA(agg.yesFGF2.Stroma, features = VariableFeatures(object = agg.yesFGF2.Stroma))
#Reruns for constructing a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset (first 40 PCs, dims = 40), how many dimentions do you want included for the analysis.
agg.yesFGF2.Stroma <- FindNeighbors(agg.yesFGF2.Stroma, dims = 1:36) 
#Reruns generation of the Clusters, resolution increases the specificity of each cluster, resolution = 2.0, the default is 1.0
agg.yesFGF2.Stroma <- FindClusters(agg.yesFGF2.Stroma, resolution = 0.8)
#Generates a t-SNE dataset, use the same dimensionality you used for FindNeighbors
agg.yesFGF2.Stroma <- RunTSNE(agg.yesFGF2.Stroma, dims = 1:36)
#Generates a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.yesFGF2.Stroma <- RunUMAP(agg.yesFGF2.Stroma, dims = 1:36)

#Saves Formated dataset after all major computational processing
saveRDS(agg.yesFGF2.Stroma, file = "Organoid_yesFGF2_Stroma_(SEURAT_v3)_01.rds")
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Stroma_Only")
agg.yesFGF2.Stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Stroma_Only/Organoid_yesFGF2_Stroma_(SEURAT_v3)_01.rds")
#Makes TSNE of stromal only data
plot1 <- DimPlot(agg.yesFGF2.Stroma, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("TSNE_yesFGF2_Stroma_ClusterNumbers_01.pdf",plot1,8,6)
#Makes UMAP of stromal only data
plot1 <- DimPlot(agg.yesFGF2.Stroma, reduction = "umap", label = TRUE, pt.size = 1.5) + NoLegend() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("UMAP_yesFGF2_Stroma_ClusterNumbers_01.pdf",plot1,8,6)

#Generates UMAP and violin plots for Stromal only data
Gene_Graphs(TGFb_super_Genes_of_Interest,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(Stromal_Markers,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(Organoid_D_Microarry_Top25_Up,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(Organoid_D_Microarry_Top25_Down,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(Secreted_Factors_AA,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(ECM_Genes_of_Interest,agg.yesFGF2.Stroma,"yesFGF2")
Gene_Graphs(Mmp_Genes_of_Interest,agg.yesFGF2.Stroma,"yesFGF2")

#Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.yesFGF2.Stroma.markers.All <- FindAllMarkers(agg.yesFGF2.Stroma, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.yesFGF2.Stroma.markers.All, file = "All_Gene_Averages_perCluster_yesFGF2_01.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.yesFGF2.Stroma.markers.Pos <- FindAllMarkers(agg.yesFGF2.Stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Alternative method for finding all Postive marker genes, saves time.
agg.yesFGF2.Stroma.markers.Pos <- subset(agg.yesFGF2.Stroma.markers.All, subset = avg_logFC > 0)
write.csv(agg.yesFGF2.Stroma.markers.Pos, file = "Positive_Gene_Averages_perCluster_yesFGF2_01.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.yesFGF2.Stroma.markers.Pos.Top10 <- agg.yesFGF2.Stroma.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.yesFGF2.Stroma.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_perCluster_yesFGF2_01.csv")

#Draws a Dot plot for stromal genes in stromal only data
plot1 <- DotPlot(object = agg.yesFGF2.Stroma, features = Stromal_Markers, cols = "Spectral", col.min = 0, dot.scale = 5) +
  RotatedAxis()
Make_pdf("DotPlot_StromalGenes_ClusterNumber_yesFGF2_01.pdf",plot1,6,4)

#New list with simple cluster names
new.cluster.ids <- c("Pa_00",
                     "Thy1_01",
                     "Pa+Pb_02",
                     "Pa_03",
                     "Pa+Pb+Thy1_04",
                     "Acta2_05",
                     "Pb_06",
                     "Pb+Thy1_07",
                     "Acta2_08",
                     "Pa_CC_09",
                     "Pa_10"
)

#Saves/catalogs cluster names to numbers
write.csv(new.cluster.ids, file = "Cluster#_Names_yesFGF2_Stroma_01.csv")
#Sets new layer with cluster names
names(new.cluster.ids) <- levels(agg.yesFGF2.Stroma)
agg.yesFGF2.Stroma <- RenameIdents(agg.yesFGF2.Stroma, new.cluster.ids)
#Creates a table showing number of cells per cluster
write.csv(table(Idents(agg.yesFGF2.Stroma)), file = "Cell_per_Cluster_yesFGF2_Stroma_01.csv")

#Generates a t-SNE plot with labels
plot1 <- DimPlot(agg.yesFGF2.Stroma, reduction = "tsne", label = TRUE, label.size = 6, repel = TRUE, pt.size = 1.5) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("TSNE_yesFGF2_Stroma_labels_01.pdf",plot1,8,5)

#Generates a UMAP plot with and without labels
plot1 <- DimPlot(agg.yesFGF2.Stroma, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5, repel = TRUE) +
          NoLegend() +
          theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
Make_pdf("UMAP_yesFGF2_Stroma_labels_01.pdf",plot1,8,6)

#Make double feature plots for markers that separate stromal clusters
plot1 <- FeaturePlot(agg.yesFGF2.Stroma, features = c("Thy1", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "Magenta", "Green"))
Make_pdf("UMAP_Plot_Thy1_Pdgfra_yesFGF2_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.yesFGF2.Stroma, features = c("Pdgfrb", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","Green"))
Make_pdf("UMAP_Plot1_Pdgfrb_Pdgfra_yesFGF2_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.yesFGF2.Stroma, features = c("Pdgfrb", "Acta2"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
Make_pdf("UMAP_Plot_Pdgfrb_aSMA_yesFGF2_01.pdf",plot1,40,8)

####################### InVivo and Oragnoid Integration, All data ###########################################
#############################################################################################################

  #DATA integration between E16 and yesFGF2 and noFGF2 organoids, All data
#base on https://rpubs.com/mathetal/integratedanalysis
#which is based on https://satijalab.org/seurat/articles/integration_rpca.html
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_All")
#Loads E16 data that is only the stomal subset
agg.E16 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/E16_SMG_Raw_(SEURAT_v4)_04.rds")
#Loads noFGF2 organoid data that is only the stomal subset
agg.noFGF2 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01/Organoid_noFGF2_Raw_(SEURAT_v4)_01.rds")
#Loads yesFGF2 organoid data that is only the stomal subset
agg.yesFGF2 <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Organoid_yesFGF2_Raw_(SEURAT_v3)_01.rds")
#Sets metadata header and value to each dataset/Seurat Object for integration function
agg.E16@meta.data[, "protocol"] <- "InVivo"
agg.noFGF2@meta.data[, "protocol"] <- "InVitro_noFGF2"
agg.yesFGF2@meta.data[, "protocol"] <- "InVitro_yesFGF2"
#integration---------------------------------------------------------------
#Merges the three datasets based on the metadata header where each cell is being described by the Protocol ID/header
agg.combined = merge(agg.E16, y = c(agg.noFGF2,agg.yesFGF2), add.cell.ids = c("InVivo","InVitro_noFGF2","InVitro_yesFGF2"), project = "protocol")
#Creates a new list of all Seurat objects based on the metadata
run.list <- SplitObject(agg.combined, split.by = "protocol")
reference.list <- run.list[c("InVitro_yesFGF2","InVitro_noFGF2","InVivo")]
# Calculates "anchors", genes that are stable that can be used to determine cell proximity to each other
# The anchor.feature default is 2000, but I found 3614 is the maximum anchors in this dataset and used as many anchors as by incrementally increasing this number until I found all possible features.
run.anchors <- FindIntegrationAnchors(object.list = reference.list,anchor.features = 3614, dims = 1:40)
# Creates an 'integrated' data assay of both datasets based on the anchors
agg.combined <- IntegrateData(anchorset = run.anchors)
# Specify that we will perform downstream analysis on the corrected data note that the original
# Sets a new default level for integrated data to be stored, Unmodified data still resides in the 'RNA' assay
DefaultAssay(agg.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
agg.combined <- ScaleData(agg.combined, verbose = FALSE)
agg.combined <- RunPCA(agg.combined, npcs = 40, verbose = FALSE)
#Used ElbowPlot to find how many PCs/dimensions to use, using 20.
ElbowPlot(agg.combined, ndims = 30, reduction = "pca")
agg.combined <- FindNeighbors(agg.combined, reduction = "pca", dims = 1:20)
agg.combined <- FindClusters(agg.combined, resolution = 0.6)
agg.combined <- RunUMAP(agg.combined, reduction = "pca", dims = 1:20)
# Saves integrated stromal data
saveRDS(agg.combined, file = "Exp_209_231_235_ALL_(SEURAT_v3)_01.rds")
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_All")
# Reads integrated stromal data
agg.combined <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_All/Exp_209_231_235_(SEURAT_v3)_01.rds")
#Graphs new data
plot1 <- DimPlot(agg.combined, reduction = "umap", group.by = "protocol")
Make_pdf("Integrated_UMAP_Sample_Label_01.pdf",plot1,6.25,5)
plot1 <- DimPlot(agg.combined, reduction = "umap", group.by = "seurat_clusters")
Make_pdf("Integrated_UMAP_Cluster_Label_01.pdf",plot1,6.25,5)
plot1 <- DimPlot(agg.combined, reduction = "umap", group.by = "seurat_clusters", split.by = "protocol")
Make_pdf("Integrated_UMAP_Cluster_Label_SampleBreakdown_01.pdf",plot1,13,5)

#Loop to create UMAP and violin plots for integrated data
Gene_Graphs(Genes_of_Interest,agg.combined,"integrated")
Gene_Graphs(Stromal_Markers,agg.combined,"integrated")
Gene_Graphs(Mmp_Genes_of_Interest,agg.combined,"integrated")
Gene_Graphs(TGFb_super_Genes_of_Interest,agg.combined,"integrated")
Gene_Graphs(ECM_Genes_of_Interest,agg.combined,"integrated")

#Saves/catalogs cluster names to cell numbers
write.csv(table(Idents(agg.combined)), file = "Cell_per_Cluster_Integrated_01.csv")
#Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.combined.markers.All <- FindAllMarkers(agg.combined, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.combined.markers.All, file = "All_Gene_Averages_Integrated_01.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.combined.markers.Pos <- FindAllMarkers(agg.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Alternative method for finding all Postive marker genes, saves time.
agg.combined.markers.Pos <- subset(agg.combined.markers.All, subset = avg_logFC > 0)
write.csv(agg.combined.markers.Pos, file = "Positive_Gene_Averages_Integrated_01.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.combined.markers.Pos.Top10 <- agg.combined.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.combined.markers.Pos.Top10, file = "Positive_Top10_Gene_Averages_Integrated_01.csv")

#Make double feature plots for markers that separate stromal clusters
plot1 <-  FeaturePlot(agg.combined, features = c("Thy1", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
Make_pdf("UMAP_plot_Pdgfra_Thy1_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.combined, features = c("Pdgfrb", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
Make_pdf("UMAP_plot_Pdgfra_Pdgfra_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.combined, features = c("Pdgfra", "Acta2"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "red","blue"))
Make_pdf("UMAP_plot_Pdgfra_aSMA_01.pdf",plot1,40,8)

############## InVivo and Oragnoid Integration, Stromal only data ###########################################
#############################################################################################################

#DATA integration between E16 and yesFGF2 and noFGF2 organoids, STROMA ONLY
#base on https://rpubs.com/mathetal/integratedanalysis
#which is based on https://satijalab.org/seurat/articles/integration_rpca.html
#Loads saved dataset
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_Stroma")
#Loads E16 data that is only the stomal subset
agg.E16.stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Stroma_only/E16_SMG_Raw_Str_(SEURAT_v4)_04.rds")
#Loads yesFGF2 organoid data that is only the stomal subset
agg.yesFGF2.stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Stroma_Only/Organoid_yesFGF2_Stroma_(SEURAT_v3)_01.rds")
#Loads noFGF2 organoid data that is only the stomal subset
agg.noFGF2.stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_231_LigatedMock_noFGF2organoid/Organoid_noFGF2/R_Results_01/Organoid_noFGF2_Stroma_(SEURAT_v3)_01.rds")
#Sets metadata header and value to each dataset/Seurat Object for integration function
agg.E16.stroma@meta.data[, "protocol"] <- "InVivo"
agg.yesFGF2.stroma@meta.data[, "protocol"] <- "InVitro_yesFGF2"
agg.noFGF2.stroma@meta.data[, "protocol"] <- "InVitro_noFGF2"
#integration---------------------------------------------------------------
#Merges the three datasets based on the metadata header where each cell is being described by the Protocol ID/header
agg.combined.stroma = merge(agg.E16.stroma, y = c(agg.noFGF2.stroma,agg.yesFGF2.stroma), add.cell.ids = c("InVivo","InVitro_noFGF2","InVitro_yesFGF2"), project = "protocol")
#Creates a new list of all Seurat objects based on the metadata
run.list.stroma <- SplitObject(agg.combined.stroma, split.by = "protocol")
reference.list.stroma <- run.list.stroma[c("InVitro_yesFGF2","InVitro_noFGF2","InVivo")]
# Calculates "anchors", genes that are stable that can be used to determine cell proximity to each other
run.anchors.stroma <- FindIntegrationAnchors(object.list = reference.list.stroma, anchor.features = 3590, dims = 1:40)
# Creates an 'integrated' data assay of both datasets based on the anchors
agg.combined.stroma <- IntegrateData(anchorset = run.anchors.stroma)
# Specify that we will perform downstream analysis on the corrected data note that the original
# Sets a new default level for integrated data to be stored, Unmodified data still resides in the 'RNA' assay
DefaultAssay(agg.combined.stroma) <- "integrated"
# Run the standard workflow for visualization and clustering
agg.combined.stroma <- ScaleData(agg.combined.stroma, verbose = FALSE)
agg.combined.stroma <- RunPCA(agg.combined.stroma, npcs = 40, verbose = FALSE)
#Used ElbowPlot to find how many PCs/dimensions to use, using 21.
ElbowPlot(agg.combined.stroma, ndims = 40, reduction = "pca")
agg.combined.stroma <- FindNeighbors(agg.combined.stroma, reduction = "pca", dims = 1:21)
agg.combined.stroma <- FindClusters(agg.combined.stroma, resolution = 0.7)
agg.combined.stroma <- RunUMAP(agg.combined.stroma, reduction = "pca", dims = 1:21)
# Saves integrated stromal data
saveRDS(agg.combined.stroma, file = "Exp_209_231_235_Stroma_(SEURAT_v3)_01.rds")
setwd("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_Stroma")
# Reads integrated stromal data
agg.combined.stroma <- readRDS("/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_235_FGF2_Organoid/R_Results_01/Integrated_E16_YesNo_FGF2_All/Exp_209_231_235_(SEURAT_v3)_01.rds")
#Graphs new data
plot1 <- DimPlot(agg.combined.stroma, reduction = "umap", group.by = "protocol")
Make_pdf("Integrated_UMAP_Sample_Label_01.pdf",plot1,6.25,5)
plot1 <- DimPlot(agg.combined.stroma, reduction = "umap", group.by = "seurat_clusters")
Make_pdf("Integrated_UMAP_Cluster_Label_01.pdf",plot1,6.25,5)
plot1 <- DimPlot(agg.combined.stroma, reduction = "umap", group.by = "seurat_clusters", split.by = "protocol")
Make_pdf("Integrated_UMAP_Cluster_Label_SampleBreakdown_01.pdf",plot1,13,5)

#Loop to create UMAP and violin plots for integrated data
Gene_Graphs(Genes_of_Interest,agg.combined.stroma,"integrated")
Gene_Graphs(Stromal_Markers,agg.combined.stroma,"integrated")
Gene_Graphs(Mmp_Genes_of_Interest,agg.combined.stroma,"integrated")
Gene_Graphs(TGFb_super_Genes_of_Interest,agg.combined.stroma,"integrated")
Gene_Graphs(ECM_Genes_of_Interest,agg.combined.stroma,"integrated")
Gene_Graphs(Cell_Cycle_Genes_of_Interest,agg.combined.stroma,"integrated")

#Saves/catalogs cluster names to cell numbers
write.csv(table(Idents(agg.combined.stroma)), file = "Cells_per_Cluster_Integrated_01.csv")
#Generating Lists of differential gene expressions
#Finds all marker genes for every cluster compared to all other clusters
agg.combined.stroma.markers.All <- FindAllMarkers(agg.combined.stroma, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agg.combined.stroma.markers.All, file = "All_DE_Gene_Averages_Integrated_01.csv")
#Finds all postive marker genes for every cluster compared to all other clusters
agg.combined.markers.Pos <- FindAllMarkers(agg.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Alternative method for finding all Postive marker genes, saves time.
agg.combined.stroma.markers.Pos <- subset(agg.combined.stroma.markers.All, subset = avg_logFC > 0)
write.csv(agg.combined.stroma.markers.Pos, file = "Positive_DE_Gene_Averages_Integrated_01.csv")
#Reduces to top 10 postive marker genes for every cluster compared to all other clusters
agg.combined.stroma.markers.Pos.Top10 <- agg.combined.stroma.markers.Pos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(agg.combined.stroma.markers.Pos.Top10, file = "Positive_DE_Top10_Gene_Averages_Integrated_01.csv")

#Gene differeneces in specific clusters enriched in certain samples
# Comparing FGF2 dominate cell cluster 0 with Control dominate cell cluster 2
yesFGF2.0_noFGF2.2_Markers <- FindMarkers(agg.combined.stroma, ident.1 = "0", ident.2 = "2")
write.csv(yesFGF2.0_noFGF2.2_Markers, file = "All_DE_Genes_Cluster0v2_01.csv")
# Comparing E16 dominate cell cluster 4 with Control dominate cell cluster 2
InVivo.4_noFGF2.2_Markers <- FindMarkers(agg.combined.stroma, ident.1 = "4", ident.2 = "2")
write.csv(yesFGF2.0_noFGF2.2_Markers, file = "All_DE_Genes_Cluster4v2_01.csv")
# Comparing E16 dominate cell cluster 4 with FGF2 dominate cell cluster 0
InVivo.4_yesFGF2.0_Markers <- FindMarkers(agg.combined.stroma, ident.1 = "4", ident.2 = "0")
write.csv(InVivo.4_yesFGF2.0_Markers, file = "All_DE_Genes_Cluster4v0_01.csv")

#Make double feature plots for markers that separate stromal clusters
plot1 <-  FeaturePlot(agg.combined.stroma, features = c("Thy1", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
Make_pdf("UMAP_plot_Pdgfra_Thy1_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.combined.stroma, features = c("Pdgfrb", "Pdgfra"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "blue","red"))
Make_pdf("UMAP_plot_Pdgfra_Pdgfra_01.pdf",plot1,40,8)
plot1 <- FeaturePlot(agg.combined.stroma, features = c("Pdgfra", "Acta2"), blend = TRUE, pt.size = 2.0, cols = c("lightgrey", "red","blue"))
Make_pdf("UMAP_plot_Pdgfra_aSMA_01.pdf",plot1,40,8)

#Extracts all integrated gene data from all cells
#Grabs all genes/anchors used in the integration
all.genes <- rownames(agg.combined.stroma)
#Pulls all the data associated with anchored genes for all cells.
All.genes.All.cells <- FetchData(object = agg.combined.stroma, vars = c("protocol","ident",all.genes), slot = "data")
write.csv(All.genes.All.cells, file = "All_Genes_All_Cells_01.csv")


#New list with simple cluster names
new.cluster.ids <- c("Pa_Pb_Thy1_00",
                     "Pa_Pb_Thy1_01",
                     "Pa_Pb_Thy1_Acta2_02",
                     "Pa_Pb_Thy1_03",
                     "Pa_Pb_Thy1_04",
                     "Acta2_05",
                     "Pb_Thy1_CC_06",
                     "Pb_Thy1_07",
                     "Pa_Pb_Thy1_Acta2_CC_08",
                     "Pa_Pb_Thy1_CC_09",
                     "Pb_Thy1_Acta2_10",
                     "Pa_Thy1_11",
                     "Pa_Pb_Thy1_12",
                     "Pb_Thy1_Acta2_CC_13",
                     "CC_14"
                     )
