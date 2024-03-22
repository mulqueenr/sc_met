```bash
#load modules
#module load nextflow/23.04.3
#module load singularity

cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cg_dataframes
singularity shell \
--bind /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2 \
~/multiome_bc.sif
```


Load data
```R
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
hundokb_rate<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_simplerate.100kb.merged.tsv.gz")
hundokb_rate<-hundokb_rate[row.names(meta)]
cg_cov_per_cell<-data.frame(cg_count=as.numeric(colSums(hundokb_counts)),sample=meta$sampleName)
plt<-ggplot(cg_cov_per_cell)+geom_density(aes(x=cg_count,color=sample))
ggsave(plt,file="cg_cov_per_cell.pdf")

cells_cov_per_window<-data.frame(cell_cov=as.numeric(rowSums(is.na(hundokb_dat))))
cells_cov_per_window$cell_cov_perc<-cells_cov_per_window$cell_cov/ncol(hundokb_dat)
cells_cov_per_window$est_rate<-as.numeric(rowMeans(hundokb_dat,na.rm=T))
cells_cov_per_window$simple_rate<-as.numeric(rowMeans(hundokb_rate,na.rm=T))
plt1<-ggplot(cells_cov_per_window)+geom_density(aes(x=cell_cov))
plt2<-ggplot(cells_cov_per_window)+geom_density(aes(x=cell_cov_perc))
plt3<-ggplot(cells_cov_per_window)+geom_density(aes(x=est_rate))
plt4<-ggplot(cells_cov_per_window)+geom_density(aes(x=simple_rate))
ggsave(plt1/plt2/plt3/plt4,file="cells_cov_per_window.pdf")
row.names(hundokb_dat)<-gsub("_", "-", row.names(hundokb_dat))
row.names(hundokb_rate)<-gsub("_", "-", row.names(hundokb_rate))

#hundokb_dat[which(is.na(hundokb_dat),arr.ind = TRUE)] <- 1.0



gene_counts<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/total_count.genebody.merged.tsv.gz")
gene_counts<-gene_counts[row.names(meta)]
gene_dat<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_posteriorest.genebody.merged.tsv.gz")
gene_dat<-gene_dat[row.names(meta)]
#gene_rate<-read_in_data("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_dataframes/mc_simplerate.genebody.merged.tsv.gz")
#gene_rate<-gene_rate[row.names(meta)]

cells_cov_per_gene<-data.frame(cell_cov=as.numeric(rowSums(is.na(gene_dat))))
cells_cov_per_gene$cell_cov_perc<-cells_cov_per_gene$cell_cov/ncol(gene_dat)
cells_cov_per_gene$est_rate<-as.numeric(rowMeans(gene_dat,na.rm=T))
#cells_cov_per_gene$simple_rate<-as.numeric(rowMeans(gene_rate,na.rm=T))
plt1<-ggplot(cells_cov_per_gene)+geom_density(aes(x=cell_cov))
plt2<-ggplot(cells_cov_per_gene)+geom_density(aes(x=cell_cov_perc))
plt3<-ggplot(cells_cov_per_gene)+geom_density(aes(x=est_rate))
#plt4<-ggplot(cells_cov_per_gene)+geom_density(aes(x=simple_rate))
ggsave(plt1/plt2/plt3,file="cells_cov_per_gene.pdf")


#CONVERTING ENSEMBL TO GENE SYMBOLS
gene_names<-AnnotationDbi::select(org.Hs.eg.db, keys = row.names(gene_dat), keytype = 'ENSEMBL', columns = 'SYMBOL',multiVals="first")
gene_names<-gene_names[which(!duplicated(gene_names$ENSEMBL)),]
gene_names<-gene_names[which(!duplicated(gene_names$SYMBOL)),]
genes.filter <- unique(gene_names[!is.na(gene_names$SYMBOL),]$ENSEMBL )
gene_dat <- gene_dat[genes.filter,]
row.names(gene_dat)<-gene_names[!is.na(gene_names$SYMBOL),]$SYMBOL
#gene_dat[which(is.na(gene_dat),arr.ind = TRUE)] <- 1.0

#CREATE SEURAT OBJECT
obj <- CreateSeuratObject(counts = as.matrix(hundokb_counts),meta.data=meta,assay="met_100kb")
obj <- SetAssayData(object = obj, layer = "scale.data", new.data = as.matrix(hundokb_dat), assay="met_100kb")
#change slot to layers for newer seurat
obj[["met_gene"]] <- CreateAssayObject(data=as.matrix(gene_dat))
obj <- SetAssayData(object = obj, slot = "scale.data", new.data = as.matrix(gene_dat),assay="met_gene")

saveRDS(obj,file="met.Rds")

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000,assay="met_100kb",layer="data")
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
```

```R
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
 
pdf("hbca_pca_feat.pdf",height=20,width=20)
FeaturePlot(dat, reduction = "umap",feature=top_5$gene,order=TRUE,ncol=5,cols=c("white","black"))
dev.off()
saveRDS(dat,file="hbca_met.Rds")

```

Trying LDA via TITAN code

```R
library(TITAN)
library(parallel)
library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

dat<-readRDS(file="met.Rds")


model_maker <- function(topics) {
    selected.Model <- lda::lda.collapsed.gibbs.sampler(
    	cellList, topics, Genes, 
    	num.iterations = iterations, 
    	alpha = alpha, eta = beta, 
    	compute.log.likelihood = TRUE, 
    	burnin = burnin)[-1]
      saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
    }

RPC_calculation<-function(model_file) {
    topic_num <- as.numeric(gsub("[^0-9]+([0-9]+).*", "\\1", 
        model_file))
    topic_numbers <- c(topic_numbers, topic_num)
    model <- readRDS(paste0(outDir, "/", model_file))
    docterMat <- t(as.matrix(data.use))
    docterMat <- as(docterMat, "sparseMatrix")
    topworddist <- normalize(model$topics, byrow = T)
    doctopdist <- normalize(t(model$document_sums), byrow = T)
    perp <- text2vec::perplexity(docterMat, topworddist, doctopdist)
    return(c(topic_num,perp))
}

runTITAN<- function(Object,
	assay="met_100kb",
	seed.number=123,
	iterations=500,
	burnin=250,
	alpha=50,
	beta=0.1,
	nfeat=10000,
	outDir="./TITAN_LDA",
	topic_counts=seq(from=5, to=80, by=5)){
	
	set.seed(seed.number)
	print(paste0("Setting ",assay," as assay..."))
	Object[[assay]]<- as(object = Object[[assay]], Class = "Assay")

	print("Filling NA values with 1.0 for Neg Binom Estimate...")
	Object[[assay]]@data[which(is.na(Object[[assay]]@data),arr.ind = TRUE)] <- 1.0
	Object <- SetAssayData(
		object = Object, 
		slot = "scale.data", 
		new.data =GetAssayData(Object, slot = "data", assay = assay),
		assay=assay)

	print(paste0("Finding ",as.character(nfeat)," variable features..."))
	var <- FindVariableFeatures(
		GetAssayData(Object, slot = "data", assay = assay), 
		selection.method = "dispersion",
		nfeatures = nfeat, 
		slot="data",
		assay=assay)

	print(paste0("Setting up data for TITAN LDA..."))
	Object[[assay]]@var.features<-row.names(var[var$variable,])
	Object.sparse <- GetAssayData(Object, slot = "data", assay = assay)
	Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assay), ]
	data.use <- Matrix::Matrix(Object.sparse, sparse = T)
	data.use <- data.use * 10 #not sure if needed
	data.use <- round(data.use) #not sure if needed
	data.use <- Matrix::Matrix(data.use, sparse = T)
	sumMat <- Matrix::summary(data.use)
	cellList <- split(as.integer(data.use@i), sumMat$j)
	ValueList <- split(as.integer(sumMat$x), sumMat$j)
	cellList <- mapply(rbind, cellList, ValueList, SIMPLIFY = F)
	Genes <- rownames(data.use)
	cellList <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1, ] + 1]; x})

	print(paste0("Running ",as.character(length(topic_counts)), " topic models..."))
	lda_out<-mclapply(topic_counts, model_maker, mc.cores = length(topic_counts))

	files <- list.files(path = outDir, pattern = "Model_")
	perp_list <- NULL
	topic_numbers <- NULL
	RPC <- NULL
	files <- files[order(nchar(files), files)]

	print(paste0("Generating perplexity estimate for model selection..."))
	perp_out<-as.data.frame(do.call("rbind",lapply(files,RPC_calculation)))
	colnames(perp_out) <- c("Topics", "RPC")
	perp_out$Topics<-as.numeric(as.character(perp_out$Topics))
	perp_out$RPC<-as.numeric(as.character(perp_out$RPC))

	perp_out$perp<-as.numeric(unlist(c("NA",lapply(2:nrow(perp_out),function(i){
		rpc_dif<-perp_out[i-1,"RPC"]-perp_out[i,"RPC"]
		topic_dif<-perp_out[i-1,"Topics"]-perp_out[i,"Topics"]
		return(abs(rpc_dif)/topic_dif)
		}))))

	#select topics from model based on elbow
	print(paste0("Found ",as.character(elbow_topic), " topics as best topic model based on elbow plot..."))
	elbow_topic<-perp_out$Topics[which(min(diff(perp_out$perp),na.rm=T)==diff(perp_out$perp))+1]
	plt1 <- ggplot(data = perp_out, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point() +geom_vline(xintercept=elbow_topic,color="red")
	plt2 <- ggplot(data = perp_out, aes(x = Topics, y = perp, group = 1)) + geom_line() + geom_point()+geom_vline(xintercept=elbow_topic,color="red")
	ggsave(plt1/plt2,file=paste0("TITAN_",assay,"_elbow.pdf"))

	top_topics<-elbow_topic #set this up as autoselect based on elbow of plt2 in future
	top_model<-readRDS(paste0(outDir, "/", "Model_",as.character(top_topics),"topics.rds"))

	print(paste0("Adding topic model to Seurat Object as lda reduction..."))
	Object <- addTopicsToSeuratObject(model = top_model, Object = Object)
	#GeneDistrubition <- GeneScores(top_model)

	print(paste0("Running UMAP and clustering..."))
	Object<-RunUMAP(Object,reduction="lda",dims=1:top_topics,reduction.name="lda_umap")
	Object <- FindNeighbors(object = Object, reduction = 'lda', dims = 1:top_topics ,graph.name="lda_snn")
	Object <- FindClusters(object = Object, verbose = TRUE, graph.name="lda_snn", resolution=0.2 ) 
	print("Plotting UMAPs...")
	plt1<-DimPlot(Object,reduction="lda_umap",group.by=c("sampleName","seurat_clusters"))
	ggsave(plt1,file=paste0("TITAN_",assay,"_umap.pdf"))
	print("Done!")
	return(Object)
}

out<-runTITAN(Object=dat,assay="met_100kb")

hbca<-subset(dat, cells=Cells(dat)[startsWith(prefix="HBCA",Cells(dat))])
out<-runTITAN(Object=hbca,assay="met_100kb")

```

CisTOPIC LDA based analysis

```R
#cistopic_generation<-function(x){
	x<-dat
	out_cistopic_obj="hbca_cistopic.rds"
	out_seurat_object="hbca_met.rds"
	cistopic_counts<-GetAssayData(object = x, layer = "scale.data",assay="met_gene")

	row.names(cistopic_counts)<-sub("-", ":", row.names(cistopic_counts))
	cistopic_counts[which(is.na(cistopic_counts),arr.ind = TRUE)] <- 1.0

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
	cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay=assay,key="topic_")
	print("Cistopic Loading into Seurat")
	x[["cistopic"]]<-cistopic_obj
	n_topics<-ncol(Embeddings(x,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
	print("Running UMAP")
	x<-RunUMAP(x,reduction="cistopic",dims=1:n_topics,reduction.name="cistopic_umap")
	x <- FindNeighbors(object = x, reduction = 'cistopic', dims = 1:n_topics ) 
	x <- FindClusters(object = x, verbose = TRUE, graph.name="met_100kb_snn", resolution=0.2 ) 
	print("Plotting UMAPs")
	plt1<-DimPlot(x,reduction="cistopic_umap",group.by=c("sampleName","seurat_clusters"))
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
	top_5<-as.data.frame(x.markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 10) %>% ungroup)
	x <- SetAssayData(object = x, layer = "scale.data", new.data = GetAssayData(object = x, layer = "data",assay="met_gene"), assay="met_gene")
	x <- ScaleData(x,assay="met_gene",do.center=T,do.scale=T)
	pdf("hbca_umap_feat.pdf",height=20,width=20)
	FeaturePlot(x, reduction = "cistopic_umap",feature=top_5$gene,order=TRUE,cols=c("blue","white","red"))*DarkTheme()
	dev.off()






specific samples in /Volumes/seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/cellranger_v3_10/BCMHBCA
ref<-readRDS("/volumes/seq/projects/metACT/ref/hbca.rds")
Idents(ref)<-ref$celltype
ref<-subset(x = ref, downsample = 2000)#just subsettting because full data set is giant


hbca_markers<-FindAllMarkers(ref)
top_5<-as.data.frame(hbca_markers %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 20) %>% ungroup)
pdf("hbca_umap_feat_posteriori.pdf",height=20,width=20)
FeaturePlot(x, reduction = "cistopic_umap",feature=top_5$gene,order=TRUE,ncol=length(unique(Idents(x))),cols=c("white","black"))
dev.off()
saveRDS(hbca_markers,file="hbca_markers.rds")



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
features<-llply(hbca_snmarkers, unlist)

x<-ScaleData(x,assay="met_gene",do.scale=TRUE,do.center=TRUE)

pdf("hbca_umap_markers.pdf",height=60,width=40)
FeaturePlot(x, reduction = "cistopic_umap",pt.size=0.75,ncol=7,feature=unlist(features),order=TRUE,cols=c("black","red","yellow"))*DarkTheme() #the *DarkThemes is a weird patchwork syntax to apply to all plots
dev.off()

bulk <- AggregateExpression(x, group.by = "seurat_clusters", return.seurat = TRUE)
bulk

#summarize gene methylation rates
cluster_rates<-gene_dat[colnames(gene_dat) %in% names(x$seurat_clusters)] %>% 
	t %>% 
	as.data.frame %>%
	group_by(x$seurat_clusters) %>%
	dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
	t %>%
	as.data.frame
colnames(cluster_rates)<-cluster_rates[1,]
cluster_rates<-cluster_rates[2:nrow(cluster_rates),]

library(ComplexHeatmap)
mat<-cluster_rates[which(row.names(cluster_rates) %in% unlist(features)),]
mat<-mat[complete.cases(mat),] %>% mutate_all(function(x) as.numeric(as.character(x))) %>% t
row.names(mat)<-colnames(cluster_rates)
pdf("hbca_umap_markers_summarized.pdf")
Heatmap(scale(mat,scale=F),
column_order=colnames(mat),
column_names_gp = grid::gpar(fontsize = 8),
row_names_gp = grid::gpar(fontsize = 8))

dev.off()

# perform standard preprocessing on each object
ref<- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)

#ALREADY RUN
#x <- NormalizeData(x)
#x<- FindVariableFeatures(x)
#x <- ScaleData(x)

# find anchors
anchors <- FindTransferAnchors(reference = ref, query = x, reference.assay ="RNA", query.assay = "met_gene", reduction="cca")

predictions <- TransferData(
  anchorset = anchors,
  refdata = ref$celltype,
  weight.reduction="cca",
)

# transfer labels
x <- AddMetaData(object = x, metadata = predictions)

pdf("hbca_umap_predictedcells.pdf",height=20,width=20)
DimPlot(x, reduction = "cistopic_umap",group.by=c("sampleName","seurat_clusters","predicted.id"))
dev.off()

pdf("hbca_umap_predictedcells_dotplot.pdf",height=10,width=50)
Idents(x)<-x$seurat_clusters
DotPlot(x,assay="met_gene",features=unlist(features),scale=F,cluster.idents=T)
dev.off()

```
#color gene plots by hbca snrna top marker genes
#try to coembedd???
#try with cistopic rather than just pca