library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-c", "--task_cpus"), type="integer", default=NULL, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd(opt$input_dir)
cpu_count=opt$task_cpus

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

dat <- runVarbin(".",
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)


dat$cell_line<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",1))
dat$well<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",2))
dat$idx<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",3))

# Mark euploid cells if they exist
dat <- findAneuploidCells(dat)

# Mark low-quality cells for filtering
dat <- findOutliers(dat)

pdf("outlier_qc.heatmap.pdf")
plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones
pdf("all_cells.subclone.umap.pdf")
plotUmap(dat, label = 'subclones')
dev.off()

pdf("all_cells.superclone.umap.pdf")
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones','cell_line','reads_total'),order='hclust')
dev.off()

pdf("all_cells.subclone.phylo.pdf")
plotPhylo(dat, label = 'subclones')
dev.off()

saveRDS(dat,file="all_cells.scCNA.rds")

```



#additional processing after moving to geo

cd /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing
cp /volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/all_cells.scCNA.rds /volumes/seq/projects/metACT/ref
mv /volumes/seq/projects/metACT/ref/all_cells.scCNA.rds /volumes/seq/projects/metACT/ref/230612_MDAMB231_SKBR3_Wafergentest.gccact.scCNA.rds

cp /volumes/seq/projects/wgd/231228_RM_WGD_ACT/scCNA.rds /volumes/seq/projects/metACT/ref
mv /volumes/seq/projects/metACT/ref/scCNA.rds /volumes/seq/projects/metACT/ref/231228_RM_WGD_ACT.mcf10a.act.scCNA.rds

singularity shell \
--bind /volumes/seq/projects/metACT \
~/singularity/copykit.sif

cd /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing

R

library(copykit)
library(BiocParallel)
library(patchwork)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing")
dat<-readRDS("all_cells.scCNA.rds")
colData(dat)$assay<-"scimet"
colData(dat)$platform <- 'combinatorial_indexing'
mdamb231_skbr3_control<-readRDS("/volumes/seq/projects/metACT/ref/230612_MDAMB231_SKBR3_Wafergentest.gccact.scCNA.rds")
colData(mdamb231_skbr3_control)$cell_line<-"MDA-MB-231"
colData(mdamb231_skbr3_control[,SummarizedExperiment::colData(mdamb231_skbr3_control)$superclones == "s1"])$cell_line<- "SKBR3"
colData(mdamb231_skbr3_control)$assay<-"gccact"
colData(mdamb231_skbr3_control)$platform<-"icell8"

mcf10a_control<-readRDS("/volumes/seq/projects/metACT/ref/231228_RM_WGD_ACT.mcf10a.act.scCNA.rds")
colData(mcf10a_control)$cell_line<-"MCF10A"
colData(mcf10a_control)$assay<-"act"
colData(mcf10a_control)$platform<-"echo"

common_coldata<-Reduce(intersect,list(colnames(colData(dat)),colnames(colData(mdamb231_skbr3_control)),colnames(colData(mcf10a_control))))
colData(dat)<-colData(dat)[common_coldata]
colData(mdamb231_skbr3_control)<-colData(mdamb231_skbr3_control)[common_coldata]
colData(mcf10a_control)<-colData(mcf10a_control)[common_coldata]
mcf10a_control@assays@data<-mcf10a_control@assays@data[1:6]

merged_copykit <- cbind(dat, mdamb231_skbr3_control,mcf10a_control)

# Create a umap embedding 
merged_copykit <- runUmap(merged_copykit)
k_clones<-findSuggestedK(merged_copykit) #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
merged_copykit  <- findClusters(merged_copykit, k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones

# Plot a copy number heatmap with clustering annotation
pdf("scalemet_and_ref.subclone.heatmap.pdf")
plotHeatmap(merged_copykit, label = c('subclones','cell_line','assay','platform','reads_total'),order='hclust',n_threads=50)
dev.off()


mdamb231 <- merged_copykit[,SummarizedExperiment::colData(merged_copykit)$cell_line == "MDA-MB-231"]

# Create a umap embedding 
mdamb231 <- runUmap(mdamb231)
k_clones<-findSuggestedK(mdamb231) #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
mdamb231  <- findClusters(mdamb231, k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones

pdf("scalemet_and_ref.mdamb231.heatmap.pdf")
plotHeatmap(mdamb231, label = c('subclones','cell_line','assay','platform','reads_total'),order='hclust',n_threads=50)
dev.off()

