library(Seurat)
library(tximport)
library(dplyr)
library(devtools)

wd <- setwd('~/trans_aligned_alevin_15000_output_transindex')

library(celda)

spleen.decon <- readRDS('spleen.soup.decon.rds')
#continuing with subset singlets
subset(spleen.decon, subset = status == 'singlet') -> singles


#find mito genes     Could be automised?
grep(rownames(singles), pattern = '-COX1-', value = T) -> Cox1
grep(rownames(singles), pattern = '-COX2-', value = T) -> Cox2
grep(rownames(singles), pattern = '-ATP8-', value = T) -> Atp8
grep(rownames(singles), pattern = '-ATP6-', value = T) -> Atp6
grep(rownames(singles), pattern = '-COX3-', value = T) -> Cox3
grep(rownames(singles), pattern = '-NU1M-', value = T) -> nu1m
grep(rownames(singles), pattern = '-NU2M-', value = T) -> nu2m
grep(rownames(singles), pattern = '-NU3M-', value = T) -> nu3m
grep(rownames(singles), pattern = '-NU4M-', value = T) -> nu4m


#not found but in mito genome
grep(rownames(singles), pattern = '-NU4LM-', value = T) -> nu4lm
grep(rownames(singles), pattern = '-NU5M-', value = T) -> nu5m
grep(rownames(singles), pattern = '-NU6M-', value = T) -> nu6m
grep(rownames(singles), pattern = '-CYB-', value = T) -> cyb
print(cyb)
#add these all together - put into a list
c(Cox1,Cox2,Atp8,Atp6,Cox3,nu1m,nu2m,nu3m,nu4m,nu5m,nu6m,cyb) -> mito.genes

#make percent.mito
percent.mito <- Matrix::colSums(GetAssayData(singles, slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(singles, slot = "counts"))

#colsums computes the sums of the matrix, takes singles
#getassadata pulls info for specified stored dimensional reduction analysis, slot specifies
#the kind of data to retrieve, singles is the(seurat) object

#add to seurat object
singles$percent.mito <- percent.mito


VlnPlot(singles, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

#more viz
plot1 <- FeatureScatter(singles, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(singles, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#remove outliers based on the earlier created scatter plots
singles <- subset(singles, subset = nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mito < 0.06)
#8002 with these stringent setting, this is up from 6k before
singles <- NormalizeData(singles, normalization.method = "LogNormalize", scale.factor = 10000)


##############
singles <- FindVariableFeatures(singles, selection.method = "vst", nfeatures = 8000)
#singles_mvp <- FindVariableFeatures(singles, selection.method = "mvp", nfeatures = 8000)

#singles_disp <- FindVariableFeatures(singles, selection.method = "disp", nfeatures = 8000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(singles), 10)

write.csv(top10, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/top10variablegenes.csv')

# plot variable features with and without labels
pdf(file="plots2.pdf")
par(mfrow=(c(1,3)))


plot1 <- VariableFeaturePlot(singles)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #label the top10 genes
plot1 + plot2

#scale the data of singles
singles <- ScaleData(singles)


#perform pca analysis based on variablefeatures 
singles <- RunPCA(singles, features = VariableFeatures(object = singles))

print(singles[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(singles, dims = 1:2, reduction = "pca")

DimPlot(singles, reduction = "pca")

DimHeatmap(singles, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(singles, dims = 1:15, cells = 500, balanced = TRUE)

dev.off()
##########################
args = commandArgs(trailingOnly=TRUE)

###deciding on the numbner of dimensions for clustering
#Elbowplot


ElbowPlot(singles, ndims = 50)

pbmc <- JackStraw(singles, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc
JackStrawPlot(pbmc, dims = 1:15)


# dims <- seq[from = args[1], to = args[2], by=args[3]]
# res <- seq[from = args[4], to = args[5], by=args[6]]

dims <- seq[from = 10, to = 40, by=5]
res <- seq[from = 0.5, to = 3, by=0.5]

source('do_scatter.R')
source('RunHarmony.R')


harmony_and_umap <- function(dims, res){
###Harmony step
pdf(file='harmony_', d, '_', r, '.pd')
par(mfrow=(c(1,3)))

library(cowplot)
library(harmony)
library(Rcpp)
source('do_scatter.R')
source('RunHarmony.R')

singles <- RunHarmony(singles, 'assignment')
singles <- FindNeighbors(singles, dims = 1:dim_val)
singles <- FindClusters(singles, resolution = res)
singles <- RunUMAP(singles, reduction = "harmony", dims = 1:dim_val) 

DimPlot(singles, reduction = "harmony", label = T)
FeaturePlot(singles, features = 'GFP') -> plot1
FeaturePlot(singles, features = 'EBFP') -> plot2
FeaturePlot(singles, features = 'CHERRY') -> plot3
DimPlot(singles, group.by = 'assignment') -> plot4 

plot_grid(plot1, plot2, plot3, plot4, ncol = 2)

dev.off()

### Finding differentially expressed features
#find markers
singles.unstranded.markers <- FindAllMarkers(singles, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.4)
t <- table(Idents(singles), singles$assignment)

#t <- table(Idents(singles), singles$assignment)
write.csv(t, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/cluster_x_assignment_', d, '_', r, '.csv')


#write into an output file
#into csv
write.csv(singles.unstranded.markers, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/markers.harmony_', d, '_', r, '.csv')

###################

singles <- FindNeighbors(singles, dims = 1:d)
singles <- FindClusters(singles, resolution = r)

head(Idents(singles), 5)

singles <- RunUMAP(singles, dims = 1:d)

pdf(file="plots3_allclusters_', d, '_', r, '.pdf")
par(mfrow=(c(1,3)))


DimPlot(singles, reduction = "umap", label = T)
FeaturePlot(singles, features = 'GFP') -> plot1
FeaturePlot(singles, features = 'EBFP') -> plot2
FeaturePlot(singles, features = 'CHERRY') -> plot3
DimPlot(singles, group.by = 'assignment') -> plot4 

library(cowplot)

#plot all 4 side by side for viewing
plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
dev.off()
}


for (d in dims){
  for(r in res){
    harmonyfunction(d,r)
  }
} 




