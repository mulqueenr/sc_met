library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="allcells", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=1, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd(opt$input_dir)
cpu_count=opt$task_cpus
prefix=opt$output_prefix

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

pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat,k=10)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".superclone.umap.pdf"))
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')

# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, label = c('superclones','subclones','cell_line','reads_total'),order='hclust')
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)