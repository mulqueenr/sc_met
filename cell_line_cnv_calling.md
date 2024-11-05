```bash
singularity shell \
--bind ~/projects \
--bind /volumes/seq/projects/metACT \
--bind /volumes/seq/projects/gccACT \
~/singularity/copykit.sif
```

```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")
scalemet<-readRDS("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/all_cells.scCNA.rds")
act_cells_dir<-"/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_nextflow/data/cells"

cpu_count=100
prefix="act_data"
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

act_dat <- runVarbin(act_cells_dir,
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)

# Mark euploid cells if they exist
act_dat <- findAneuploidCells(act_dat )

# Mark low-quality cells for filtering
act_dat <- findOutliers(act_dat)

pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(act_dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
act_dat <- knnSmooth(act_dat)

# Create a umap embedding 
act_dat <- runUmap(act_dat)
k_clones<-findSuggestedK(act_dat) #10

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
act_dat  <- findClusters(act_dat, 
                    k_superclones=k_clones@metadata$suggestedK, 
                    k_subclones=k_clones@metadata$suggestedK+10)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(act_dat, label = 'subclones')
plotUmap(act_dat, label = 'superclones')
dev.off()

# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(act_dat, label = c('superclones','subclones','reads_total'),order='hclust',n_threads=100)
dev.off()

act_dat$cell_line<-ifelse(act_dat$superclones %in% c("s4","s3"),"MDA-MB-231","SKBR-3")
saveRDS(act_dat,file=paste0(prefix,".scCNA.rds"))

#PROCESS MERGED FILES ALL CELL LINES
prefix="all_cellline_withACT"
scalemet<-readRDS("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/all_cells.scCNA.rds")
act_dat<-readRDS(file="act_data.scCNA.rds")
scalemet@colData<-scalemet@colData[colnames(scalemet@colData) %in% colnames(act_dat@colData)]
scalemet$method<-"scale_met"
act_dat$method<-"ACT"
dat_merged<-cbind(scalemet,act_dat)
dat_merged<-dat_merged[,colData(dat_merged)$cell_line %in% c("MCF10A","MCF7","MDA-MB-231","SKRB-3")]
dat_merged <- knnSmooth(dat_merged)
dat_merged <- runUmap(dat_merged)
k_clones<-findSuggestedK(dat_merged) #10
dat_merged  <- findClusters(dat_merged, 
                    k_superclones=k_clones@metadata$suggestedK+20, 
                    k_subclones=k_clones@metadata$suggestedK+30)#output from k_clones
saveRDS(dat_merged,file=paste0(prefix,".scCNA.rds"))

pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat_merged, label = 'subclones')
plotUmap(dat_merged, label = 'superclones')
dev.off()

pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat_merged, label = c('superclones','subclones','reads_total','cell_line','method'),order='hclust',n_threads=100)
dev.off()

dat_merged <- runMetrics(dat_merged)
plt<-ggplot(dat_merged@colData,aes(x=paste(cell_line,method),y=dat_merged@colData$overdispersion,fill=paste(cell_line,method),color=paste(cell_line,method)))+geom_violin()+geom_jitter()
ggsave(plt,file="overdispersion.pdf")

#PROCESS MERGED FILES JUST 231
prefix="mda_mb_231_withACT"
scalemet<-readRDS("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/all_cells.scCNA.rds")
act_dat<-readRDS(file="act_data.scCNA.rds")
scalemet@colData<-scalemet@colData[colnames(scalemet@colData) %in% colnames(act_dat@colData)]
scalemet$method<-"scale_met"
act_dat$method<-"ACT"
dat_merged<-cbind(scalemet,act_dat)
dat_merged<-dat_merged[,colData(dat_merged)$cell_line %in% c("MDA-MB-231")]
dat_merged <- runMetrics(dat_merged)
dat_merged<-dat_merged[,dat_merged@colData$overdispersion<0.03]

dat_merged <- knnSmooth(dat_merged)
dat_merged <- runUmap(dat_merged)
k_clones<-findSuggestedK(dat_merged) #10
dat_merged  <- findClusters(dat_merged, 
                    k_superclones=k_clones@metadata$suggestedK+20, 
                    k_subclones=k_clones@metadata$suggestedK+30)
saveRDS(dat_merged,file=paste0(prefix,".scCNA.rds"))


pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat_merged, label = 'subclones')
plotUmap(dat_merged, label = 'superclones')
dev.off()

pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat_merged, label = c('superclones','subclones','reads_total','cell_line','method'),order='hclust',n_threads=100)
dev.off()


#PROCESS MERGED FILES JUST HBCA
prefix="hbca"
scalemet<-readRDS("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/all_cells.scCNA.rds")
scalemet$method<-"scale_met"
dat_merged<-scalemet
dat_merged<-dat_merged[,startsWith(colData(dat_merged)$cell_line,prefix="HBCA")]
dat_merged <- runMetrics(dat_merged)
dat_merged<-dat_merged[,dat_merged@colData$overdispersion<0.03]
dat_merged <- knnSmooth(dat_merged)
dat_merged <- runUmap(dat_merged)
k_clones<-findSuggestedK(dat_merged) #10
dat_merged  <- findClusters(dat_merged, 
                    k_superclones=k_clones@metadata$suggestedK+20, 
                    k_subclones=k_clones@metadata$suggestedK+30)
saveRDS(dat_merged,file=paste0(prefix,".scCNA.rds"))


pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat_merged, label = 'subclones')
plotUmap(dat_merged, label = 'superclones')
dev.off()

pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat_merged, label = c('superclones','subclones','reads_total','cell_line','method'),order='hclust',n_threads=100)
dev.off()

