library(Seurat)
library(data.table)

library(ggplot2)
library(patchwork)

meta_files<-list.files(recursive=TRUE,path=".",pattern="*.passingCellsMapMethylStats.csv",full.names=TRUE)


met<-do.call("rbind",lapply(meta_files,function(x) read.csv(sep=",",file=x,header=T)))
row.names(met)<-paste0(met$sampleName,met$BC)



read_in_data<-function(x){
	dat<-as.data.frame(fread(x,showProgress=TRUE,nThread=5))
	row.names(dat)<-dat[,1]
	dat<-dat[2:ncol(dat)]
	return(dat)
}

hundokb_counts<-read_in_data("total_count.100kb.merged.tsv.gz")
hundokb_dat<-read_in_data("mc_posteriorest.100kb.merged.tsv.gz")
gene_counts<-read_in_data("total_count.genebody.merged.tsv.gz")
gene_dat<-read_in_data("mc_posteriorest.genebody.merged.tsv.gz")

seurat_object <- CreateSeuratObject(counts = t(hundokb_counts),meta.data=met,assay="met_100kb")
seurat_object[["met_100kb"]]
seurat_object[["met_gene"]] <- CreateAssayObject(counts=t(gene_counts),data=t(gene_dat))

saveRDS(seurat_object,file="met.SeuratObject.rds")

seurat_object <- FindVariableFeatures(seurat_object, selection.method = "dispersion", nfeatures = 10000,assay="met_100kb")
seurat_object <- RunPCA(seurat_object, assay="met_100kb",features = VariableFeatures(object = seurat_object))

pca_plot<-DimHeatmap(seurat_object, dims = 1:10, cells = 500, balanced = TRUE)
elbow_plot<-ElbowPlot(seurat_object)
ggsave(pca_plot/elbow_plot,file="met_pca__100kb_10dims.pdf")
