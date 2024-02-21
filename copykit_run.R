library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
library(scquantum)
register(MulticoreParam(progressbar = T, workers = 5), default = T)
BiocParallel::bpparam()

setwd("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDedup")
dat <- runVarbin(".",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

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

dat$tn5<-unlist(lapply(strsplit(dat$sample,"[.]"),"[",2))
dat$tn5_plate<-substr(dat$tn5,1,1)
dat$tn5_row<-substr(dat$tn5,2,2)
dat$tn5_col<-substr(dat$tn5,3,4)

dat$cell_line<-"NA"
dat@colData[which(dat@colData$tn5_col %in% c("05","06","07","08","09","10")),]$cell_line<-"mda_mb_231"
dat@colData[which(dat@colData$tn5_col %in% c("01","02","03","04")),]$cell_line<-"mcf10a"
dat@colData[which(dat@colData$tn5_col %in% c("11","12")),]$cell_line<-"gm12878"

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones','tn5_col','cell_line','reads_total'),order='hclust')
dev.off()

pdf("all_cells.subclone.phylo.pdf")
plotPhylo(dat, label = 'subclones')
dev.off()

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones',),order='hclust')
dev.off()

saveRDS(dat,file="all_cells.scCNA.rds")

colData(dat)$info<- 'scalemet'

################################################
###Merge with standard act seq for comparison###
################################################
#scalebio dat

dat <- runVarbin(".",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

dat$tn5<-unlist(lapply(strsplit(dat$sample,"[.]"),"[",2))
dat$tn5_plate<-substr(dat$tn5,1,1)
dat$tn5_row<-substr(dat$tn5,2,2)
dat$tn5_col<-substr(dat$tn5,3,4)
dat$info <- 'scalebio'

dat$cell_line<-"NA"
dat@colData[which(dat@colData$tn5_col %in% c("05","06","07","08","09","10")),]$cell_line<-"mda_mb_231"
dat@colData[which(dat@colData$tn5_col %in% c("01","02","03","04")),]$cell_line<-"mcf10a"
dat@colData[which(dat@colData$tn5_col %in% c("11","12")),]$cell_line<-"gm12878"

#act set
act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
colData(act)$info <- 'act'
colData(act)$cell_line<-"mda_mb_231"

colData(dat)<-colData(dat)[colnames(colData(dat)) %in% colnames(colData(act))]
colData(act)<-colData(act)[colnames(colData(act)) %in% colnames(colData(dat))]

merged_copykit <- cbind(act, dat)

# Mark euploid cells if they exist
merged_copykit <- findAneuploidCells(merged_copykit) 
merged_copykit <- findOutliers(merged_copykit)
merged_copykit <- knnSmooth(merged_copykit)
merged_copykit <- runUmap(merged_copykit)
k_clones<-findSuggestedK(merged_copykit)
merged_copykit <- findClusters(merged_copykit ,k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones

pdf("act_merged.subclone.umap.pdf")
plotUmap(merged_copykit, label = 'subclones')
dev.off()

pdf("act_merged.cellline.umap.pdf")
plotUmap(merged_copykit, label = 'cell_line')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
merged_copykit <- calcConsensus(merged_copykit)
merged_copykit <- runConsensusPhylo(merged_copykit)
merged_copykit <- runPhylo(merged_copykit, metric = 'manhattan')

# Plot a copy number heatmap with clustering annotation
pdf("act_merged.subclone.heatmap.pdf")
plotHeatmap(merged_copykit, label = c('superclones','subclones','cell_line','info','reads_total'),order='hclust')
dev.off()

saveRDS(dat,file="act_merged.scCNA.rds")
