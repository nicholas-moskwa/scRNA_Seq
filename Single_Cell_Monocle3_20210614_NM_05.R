#Version 4
#This verion is meant to look at Exp_209 Stromal subset data integrated together and where proliferating cells are regressed and re-scaled to diminish cell cycling genes contributions to clustering.
#scRNA data is  calculated and subseted using Seurat
#Trajectory done using Monocle3 function to map cell differentiation trajectories
#Need Seurat v4.0 or higher
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(magrittr)
library(dplyr)
library(rlang)
library(bslib)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(SeuratData)

#Downloads data from Single_cell_Analysis_20210216_NM_04
#Loads saved dataset
#Exp_209 only No cell cycle regression
agg.E16.Stroma <- readRDS("Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Stroma_only/E16_SMG_Raw_Str_(SEURAT_v4)_04.rds")
#Exp_209 only Yes cell cycle regression
agg.E16.Stroma <- readRDS("Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Stroma_only/E16_SMG_Raw_Str_CC_Regressed(SEURAT_v4)_04.rds")
#Exp_209 and Exp_223 only Yes cell cycle regression
agg.E16.Stroma <- readRDS("Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/R_results_04/Integrated_Exp_209_223/E16_SMG_Integrated_CC_Regressed_Stroma(SEURAT_v4)_04.rds")

#Sets the working directory
setwd("Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_209_E16_Str_Isolation_scRNA_202001026/Monocle3_05")
#Plot showing data looks good.
DimPlot(agg.E16.Stroma, reduction = "umap", label = TRUE, pt.size = 1.5) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
#Creates a new object for Monocle3 from Seurat object
cds.2 <- as.cell_data_set(agg.E16.Stroma)
#Calculates cell clusters using "community detection" method based on Leiden algorithm
cds.2 <- cluster_cells(cds.2, reduction_method = "UMAP")
#Calculate size factors using built-in function in monocle3
cds.2 <- estimate_size_factors(cds.2) # check this one https://github.com/satijalab/seurat-wrappers/issues/54
#Add gene names into CDS
cds.2@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(agg.E16.Stroma[["RNA"]])
#Add gene names into CDS, needs to be "integrated" if using combined datasets
cds.2@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(agg.E16.Stroma[["integrated"]])
# Generates 2-axix plot using UMAP dimentions
plot_cells(cds.2)
# Generates 2-axix plot using UMAP dimentions and colorizes based on genes
plot_cells(cds.2, genes=c("Acta2"))
#Creates pseudotime trajectories
cds.2 <- learn_graph(cds.2)
# Generates 2-axix plot using UMAP dimentions but with psuedotime tracks
pdf("UMAP_Pseudotime_No_Labels_E16_All_02.pdf", width = 10, height = 10)
    plot_cells(cds.2)
dev.off()
#Find root nodes for pseudotime prediction
cds.2 <- order_cells(cds.2, reduction_method = "UMAP")
#Draws pseudotime prediction plot based on previous selection
plot1 <- plot_cells(cds.2, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, show_trajectory_graph = TRUE, graph_label_size = 4, cell_size = 1) +
                    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
pdf("UMAP_Pseudotime_Stroma_E16_Exp209-223_Thy1-Start_CC-Regressed_05.pdf", width = 6, height = 4)
    plot1
dev.off()
#Stromal marker list
Stromal_Markers <- c("Thy1", "Acta2","Pdgfrb", "Pdgfra", "Vim", "Col1a1")
#Loop to create Feature plots for Gene lists
for(i in Stromal_Markers[1:6]){
  print(i)
  plot3 <- plot_cells(cds.2, genes = i, graph_label_size = 4, cell_size = 1) +
                      theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18))
  Dot_Plot_Title <- str_c("UMAP_Plot_", i, "_E16_Exp209-223_CC-Regressed_05.pdf")
  pdf(file = Dot_Plot_Title, width = 6, height = 4)
  print(plot3)
  dev.off()  
}
#Calculates the top genes that define each Monocle3 cluster
marker_test_res <- top_markers(cds.2,
                               group_cells_by="seurat_clusters",
                               reference_cells=1000,
                               cores=3)
#Creates a subset of unique markers, used for making figures
top_specific_markers <- marker_test_res %>%
                        filter(fraction_expressing >= 0.10) %>%
                        group_by(cell_group) %>%
                        top_n(4, pseudo_R2)
#Pulls the geneID for each unique marker
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
#Creates a dotplot showing the most highly changing unique genes for each Monocle3 cluster
plot1 <- plot_genes_by_group(cds.2,
                             top_specific_marker_ids,
                             group_cells_by="seurat_clusters",
                             ordering_type="maximal_on_diag",
                             max.size=5)
#Saves dotplot of top unqiue genes per cluster
pdf("DotPlot_Top4_E16_Exp209-223_CC-Regressed_05.pdf", width = 8, height = 10)
    plot1
dev.off()
#Saves dotplot of genes per cluster using Stromal_Markers from Seurat script
plot1 <- plot_genes_by_group(cds.2,
                            Stromal_Markers,
                            group_cells_by="seurat_clusters",
                            ordering_type="cluster_row_col",
                            max.size=3)
pdf("DotPlot_GenesofInterest_E16_Exp209-223_CC-Regressed_05.pdf", width = 8, height = 10)
    plot1
dev.off()
