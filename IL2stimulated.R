library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc <- readRDS("GSE144191 2")

pbmc <- UpdateSeuratObject(pbmc)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#TSNE
# If you haven't installed TSNE, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunTSNE(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "tsne")
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)    



# purple plots
FeaturePlot(pbmc, features = c("KLRC1", "KLRK1", "FCGR3A", "CD1A", "IL2", "IL2RA", "GZMB", "IFNG",
                               "TNF", "NCAM1", "PRF1"))
FeaturePlot(pbmc, features = c("D1DR", "DRD2", "DRD3", "DRD4", "DRD5"))
FeaturePlot(pbmc, features = c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4"))
FeaturePlot(pbmc, features = c("HTR5A", "HTR5BP", "HTR6", "HTR7", "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E"))
FeaturePlot(pbmc, features = c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3")) 
FeaturePlot(pbmc, features = c("NR3C1"))
FeaturePlot(pbmc, features = c("BDKBRB1", "BDKBRB2", "GALR1", "NMUR1", "NPY1R", "OPRM1", "TACR1"))
FeaturePlot(pbmc, features = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B"))
FeaturePlot(pbmc, features = c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"))
FeaturePlot(pbmc, features = c("CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10"))
FeaturePlot(pbmc, features = c("CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", "CHRND", "CHRNE", "CHRNG"))
FeaturePlot(pbmc, features = c("GRM1", "GRM2","GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8"))
FeaturePlot(pbmc, features = c("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRID1", "GRID2"))

#TSNE
new.cluster.ids <- c("KLRC1", "KLRK1", "FCGR3A", "CD1A", "IL2", "IL2RA", "GZMB", "IFNG",
                     "TNF", "NCAM1", "PRF1")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


# violin plots
VlnPlot(pbmc, features = c("KLRC1", "KLRK1", "FCGR3A", "CD1A", "IL2", "IL2RA", "GZMB", "IFNG",
                           "TNF", "NCAM1", "PRF1"))
VlnPlot(pbmc, features = c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3")) 
VlnPlot(pbmc, features = c("D1DR", "DRD2", "DRD3", "DRD4", "DRD5"))
VlnPlot(pbmc, features = c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4"))
VlnPlot(pbmc, features = c("HTR5A", "HTR5BP", "HTR6", "HTR7", "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E"))
VlnPlot(pbmc, features = c("KLRC1", "KLRK1", "FCGR3A", "CD1A", "IL2", "IL2RA", "GZMB", "IFNG",
                           "TNF", "NCAM1"))
VlnPlot(pbmc, features = c("NR3C1"))
VlnPlot(pbmc, features = c("BDKBRB1", "BDKBRB2", "GALR1", "NMUR1", "NPY1R", "OPRM1", "TACR1"))
VlnPlot(pbmc, features = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B"))
VlnPlot(pbmc, features = c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"))
VlnPlot(pbmc, features = c("CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10"))
VlnPlot(pbmc, features = c("CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", "CHRND", "CHRNE", "CHRNG"))
VlnPlot(pbmc, features = c("GRM1", "GRM2","GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8"))
VlnPlot(pbmc, features = c("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRID1", "GRID2"))


#Length
pbmc.adrb2pos <- subset(pbmc, subset = ADRB2 > 0)
pbmc.nr3c1pos <- subset(pbmc, subset = NR3C1 > 0)
pbmc.nmur1pos <- subset(pbmc, subset = NMUR1 > 0)
pbmc.chrnb1pos <- subset(pbmc, subset = CHRNB1 > 0)

#number
length(pbmc.adrb2pos@assays$RNA@counts@p)-1
length(pbmc.nr3c1pos@assays$RNA@counts@p)-1
length(pbmc.nmur1pos@assays$RNA@counts@p)-1
length(pbmc.chrnb1pos@assays$RNA@counts@p)-1

#percentage
(length(pbmc.adrb2pos@assays$RNA@counts@p)-1)/(length(pbmc@assays$RNA@counts@p)-1)
(length(pbmc.nr3c1pos@assays$RNA@counts@p)-1)/(length(pbmc@assays$RNA@counts@p)-1)
(length(pbmc.nmur1pos@assays$RNA@counts@p)-1)/(length(pbmc@assays$RNA@counts@p)-1)
(length(pbmc.chrnb1pos@assays$RNA@counts@p)-1)/(length(pbmc@assays$RNA@counts@p)-1)



