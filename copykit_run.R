library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 60), default = T)
BiocParallel::bpparam()
library(optparse)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="List of single-cell bam files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd(opt$input_dir)

dat <- runVarbin(".",
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)

dat$cell_line<-unlist(lapply(strsplit(dat$samples,"[.]"),"[[",1))

# Mark euploid cells if they exist
dat <- findAneuploidCells(dat)

# Mark low-quality cells for filtering
dat <- findOutliers(dat)

# Visualize cells labeled by filter and aneuploid status
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
dat  <- findClusters(dat,k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones
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

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones','cell_line'),order='hclust')
dev.off()

saveRDS(dat,file="all_cells.scCNA.rds")