library(Seurat)
library(ggplot2)
setwd("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/matrix")

cg_rat<-read.table("HBCA-16R.CG.ratio.mtx")
cg_cov<-read.table("HBCA-16R.CG.cov.mtx")

dat <- CreateSeuratObject(counts = cg_cov, data=cg_rat)
dat<- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 20000)

# Identify the 10 most highly variable genes

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
dat <- ScaleData(dat)
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
plot2<-ElbowPlot(dat)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.5)
dat <- RunUMAP(dat, dims = 1:10)
plot3<-DimPlot(dat, reduction = "umap")
plot4<-FeaturePlot(dat, reduction = "umap",feature="nCount_RNA")



out_plot<-plot1/plot2/plot3/plot4
ggsave(out_plot,file="test.pdf")
