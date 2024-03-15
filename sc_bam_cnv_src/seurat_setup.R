#load modules
#module load nextflow/23.04.3
#module load singularity

cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cg_dataframes
singularity shell \
--bind /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2 \
~/multiome_bc.sif

library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)
library(dplyr)
library(cisTopic)

#LOAD IN REFERENCE HBCA 
#ref<-readRDS("/volumes/seq/projects/metACT/ref/hbca.rds")

#LOAD IN META FILES
meta_files<-list.files(recursive=TRUE,path="/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/report",pattern="*.passingCellsMapMethylStats.csv",full.names=TRUE)
meta<-do.call("rbind",lapply(meta_files,function(x) read.csv(sep=",",file=x,header=T)))
row.names(meta)<-paste0(meta$sampleName,meta$BC)


#LOAD IN ALL DATA
read_in_data<-function(x){
	dat<-as.data.frame(fread(x,showProgress=TRUE,nThread=20))
	sample<-unlist(lapply(strsplit(dat[,1],"_"),"[",1))
	bc<-unlist(lapply(dat[,1],function(x) substr(x,nchar(x)-29,nchar(x))))
	row.names(dat)<-paste0(sample,bc)
	dat<-dat[2:ncol(dat)]
	dat<-as.data.frame(t(dat))
	return(dat)
}

hundokb_counts<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/total_count.100kb.merged.tsv.gz")
hundokb_counts<-hundokb_counts[row.names(meta)]
hundokb_dat<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_posteriorest.100kb.merged.tsv.gz")
hundokb_dat<-hundokb_dat[row.names(meta)]
hundokb_dat[which(is.na(hundokb_dat),arr.ind = TRUE)] <- 1.0
row.names(hundokb_dat)<-gsub("_", "-", row.names(hundokb_dat))

#hundokb_rate<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_simplerate.100kb.merged.tsv.gz")
#hundokb_rate<-hundokb_rate[row.names(meta)]

gene_counts<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/total_count.genebody.merged.tsv.gz")
gene_dat<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_posteriorest.genebody.merged.tsv.gz")
gene_dat<-gene_dat[row.names(meta)]

#CONVERTING ENSEMBL TO GENE SYMBOLS
gene_names<-AnnotationDbi::select(org.Hs.eg.db, keys = row.names(gene_dat), keytype = 'ENSEMBL', columns = 'SYMBOL',multiVals="first")
gene_names<-gene_names[which(!duplicated(gene_names$ENSEMBL)),]
gene_names<-gene_names[which(!duplicated(gene_names$SYMBOL)),]
genes.filter <- unique(gene_names[!is.na(gene_names$SYMBOL),]$ENSEMBL )
gene_dat <- gene_dat[genes.filter,]
row.names(gene_dat)<-gene_names[!is.na(gene_names$SYMBOL),]$SYMBOL
gene_dat[which(is.na(gene_dat),arr.ind = TRUE)] <- 1.0

#CREATE SEURAT OBJECT
obj <- CreateSeuratObject(counts = as.matrix(hundokb_counts),meta.data=meta,assay="met_100kb")
obj <- SetAssayData(object = obj, slot = "scale.data", new.data = as.matrix(hundokb_dat), assay="met_100kb")
#change slot to layers for newer seurat
obj[["met_gene"]] <- CreateAssayObject(data=as.matrix(gene_dat))
obj <- SetAssayData(object = obj, slot = "scale.data", new.data = as.matrix(gene_dat),assay="met_gene")

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000,assay="met_100kb")
pdf("met_variable_100kb.pdf",width=10)
VariableFeaturePlot(obj,selection.method="vst")
dev.off()

obj <- RunPCA(obj, assay="met_100kb",features = VariableFeatures(object = obj))

pdf("met_pca_heatmap.pdf")
pca_plot<-DimHeatmap(obj, dims = 1:5, cells = 500, balanced = TRUE,assay="met_100kb")
dev.off()

elbow_plot<-ElbowPlot(obj)
ggsave(elbow_plot,file="met_pca_gene_elbow.pdf")

obj <- FindNeighbors(obj, dims = 1:5,assay="met_100kb")
obj <- FindClusters(obj, resolution = 0.2,assay="met_100kb")
obj <- RunUMAP(obj, dims = 1:5,assay="met_100kb")
obj <- RunTSNE(obj, dims = 1:5,assay="met_100kb")

pdf("met_umap.pdf",width=10)
DimPlot(obj, reduction = "umap",group.by=c("sampleName","seurat_clusters"))
dev.off()

pdf("met_tsne.pdf",width=10)
DimPlot(obj, reduction = "tsne",group.by=c("sampleName","seurat_clusters"))
dev.off()

Idents(obj)<-obj$sampleName
obj.markers <- FindAllMarkers(
	obj, 
	assay="met_gene",
	slot="scale.data",
	only.pos = TRUE,
	test.use = "wilcox", 
	min.pct=0,
	logfc.threshold = 0)
top_5<-as.data.frame(obj.markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup)
 
pdf("met_umap_feat.pdf",height=20,width=20)
FeaturePlot(obj, reduction = "umap",feature=top_5$gene,order=TRUE,ncol=5,cols=c("white","black"))
dev.off()
saveRDS(obj,file="hbca_met.Rds")

###just HBCA
dat<-subset(obj, cells=Cells(obj)[startsWith(prefix="HBCA",Cells(obj))])

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 5000,assay="met_100kb")
pdf("hbca_variable_100kb.pdf",width=10)
VariableFeaturePlot(dat,selection.method="vst")
dev.off()

dat <- RunPCA(dat, assay="met_100kb",features = VariableFeatures(object = dat))

pdf("hbca_pca_heatmap.pdf")
pca_plot<-DimHeatmap(dat, dims = 1:5, cells = 500, balanced = TRUE,assay="met_100kb")
dev.off()

elbow_plot<-ElbowPlot(dat)
ggsave(elbow_plot,file="met_pca_gene_elbow.pdf")

dat <- FindNeighbors(dat, dims = 1:5,assay="met_100kb")
dat <- FindClusters(dat, resolution = 0.2,assay="met_100kb")
dat <- RunUMAP(dat, dims = 1:5,assay="met_100kb")
dat <- RunTSNE(dat, dims = 1:5,assay="met_100kb")

pdf("hbca_umap.pdf",width=10)
DimPlot(dat, reduction = "umap",group.by=c("sampleName","seurat_clusters"))
dev.off()

pdf("hbca_tsne.pdf",width=10)
DimPlot(dat, reduction = "tsne",group.by=c("sampleName","seurat_clusters"))
dev.off()

Idents(dat)<-dat$seurat_clusters
dat.markers <- FindAllMarkers(
	dat, 
	assay="met_gene",
	slot="scale.data",
	only.pos = TRUE,
	test.use = "wilcox", 
	min.pct=0,
	logfc.threshold = 0)
top_5<-as.data.frame(dat.markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup)
 
pdf("hbca_umap_feat.pdf",height=20,width=20)
FeaturePlot(dat, reduction = "umap",feature=top_5$gene,order=TRUE,ncol=5,cols=c("white","black"))
dev.off()
saveRDS(dat,file="hbca_met.Rds")



cistopic_generation<-function(x){
	x<-dat
	out_cistopic_obj="hbca_cistopic.rds"
	out_seurat_object="hbca_met.rds"
	cistopic_counts<-GetAssayData(object = x, layer = "scale.data")
	row.names(cistopic_counts)<-sub("-", ":", row.names(cistopic_counts))
	#cistopic assumes beta values, inputting directly
	sub_cistopic <- cisTopic::createcisTopicObject(count.matrix = cistopic_counts, project.name = "met", min.cells = 0, min.regions = 0)
	print("Made cistopic object.")
	sub_cistopic_models<-cisTopic::runModels(sub_cistopic,
	topic=seq(from=5, to=50, by=10),
	nCores=5,
	addModels=FALSE) #using v2 of cistopic (we are only using single core anyway)
	saveRDS(sub_cistopic_models,file=out_cistopic_obj)

	sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =x@meta.data)
	sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')

	print("Finshed running cistopic.")

	#Add cell embeddings into seurat
	cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
	colnames(cell_embeddings)<-sub_cistopic_models@cell.names
	n_topics<-nrow(cell_embeddings)
	row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
	cell_embeddings<-as.data.frame(t(cell_embeddings))

	#Add feature loadings into seurat
	feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
	row.names(feature_loadings)<-paste0("topic_",1:n_topics)
	feature_loadings<-as.data.frame(t(feature_loadings))

	#combined cistopic results (cistopic loadings and umap with seurat object)
	cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="met_100kb",key="topic_")
	print("Cistopic Loading into Seurat")
	x[["cistopic"]]<-cistopic_obj
	n_topics<-ncol(Embeddings(x,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
	print("Running UMAP")
	x<-RunUMAP(x,reduction="cistopic",dims=1:n_topics,reduction.name="cistopic_umap")
	x <- FindNeighbors(object = x, reduction = 'cistopic', dims = 1:n_topics ) 
	x <- FindClusters(object = x, verbose = TRUE, graph.name="met_100kb_snn", resolution=0.2 ) 
	print("Plotting UMAPs")
	plt1<-DimPlot(x,reduction="umap",group.by=c("sampleName","seurat_clusters"))
	pdf("hbca.cistopic.umap.pdf",width=10)
	print(plt1)
	dev.off()
	saveRDS(sub_cistopic_models,file=out_cistopic_obj)
	saveRDS(x,file=out_seurat_object)

	Idents(x)<-x$seurat_clusters
	x.markers <- FindAllMarkers(
		x, 
		assay="met_gene",
		slot="scale.data",
		only.pos = TRUE,
		test.use = "wilcox", 
		min.pct=0,
		logfc.threshold = 0)
	top_5<-as.data.frame(x.markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup)
	 
	pdf("hbca_umap_feat.pdf",height=20,width=20)
	FeaturePlot(x, reduction = "umap",feature=top_5$gene,order=TRUE,ncol=length(unique(Idents(x))),cols=c("white","black"))
	dev.off()

ref<-readRDS("/volumes/seq/projects/metACT/ref/hbca.rds")
Idents(ref)<-ref$celltype
hbca_markers<-FindAllMarkers(ref)
top_5<-as.data.frame(hbca_markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup)





library(plyr)
#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-unlist(llply(hbca_snmarkers, unlist))

pdf("hbca_umap_markers.pdf",height=40,width=40)
FeaturePlot(x, reduction = "umap",pt.size=0.5,feature=features,order=TRUE,ncol=length(unique(Idents(x))),cols=c("#5e4fa2","#ffffbf","#9e0142"))
dev.off()

#color gene plots by hbca snrna top marker genes
#try to coembedd???
#try with cistopic rather than just pca