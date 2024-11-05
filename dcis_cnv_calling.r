```bash
singularity shell \
--bind ~/projects \
--bind /volumes/seq/projects/metACT \
--bind /volumes/seq/projects/gccACT \
~/singularity/copykit.sif
```

## ALL HBCA
```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)
library(colorRamp2)
cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")
prefix<-"HBCA_All"
hbca_16_83 <- runVarbin("/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/sc_bam",
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
hbca_16_83<-hbca_16_83[,endsWith(hbca_16_83@colData$sample,suffix=".sorted")]
colData(hbca_16_83)$cell_line<-unlist(lapply(strsplit(hbca_16_83@colData$sample,"[.]"),"[",1))
hbca_16_83<-hbca_16_83[,hbca_16_83@colData$cell_line %in% c("HBCA-16R","HBCA-83L")]
cells_dir<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup"
list.files(path=cells_dir,pattern="_split_bam")
hbca_17 <- runVarbin(paste(cells_dir,"HBCA-17T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
hbca_19 <- runVarbin(paste(cells_dir,"HBCA-19T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
colData(hbca_17)$cell_line<-"HBCA-17L"
colData(hbca_19)$cell_line<-"HBCA-19L"
colData(hbca_16_83)<-colData(hbca_16_83)[,colnames(colData(hbca_16_83 )) %in% colnames(colData(hbca_17))]
colData(hbca_17)<-colData(hbca_17)[,colnames(colData(hbca_17)) %in% colnames(colData(hbca_16_83 ))]
colData(hbca_19)<-colData(hbca_19)[,colnames(colData(hbca_19)) %in% colnames(colData(hbca_16_83 ))]
dat<-cbind(hbca_17,hbca_19)
dat<-cbind(dat,hbca_16_83)
dat  <- findAneuploidCells(dat  )
dat <- findOutliers(dat)
dat<- runMetrics(dat)
dat<-dat[,dat@colData$overdispersion<0.03]
pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('outlier', 'is_aneuploid'), 
        row_split = 'outlier')
dev.off()
dat<- knnSmooth(dat)
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10
dat <- findClusters(dat, 
                    k_superclones=k_clones@metadata$suggestedK+20, 
                    k_subclones=k_clones@metadata$suggestedK+30)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
plotUmap(dat, label = 'superclones')
dev.off()
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('superclones','subclones','reads_total','cell_line'),
        order='hclust',
        col=colorRamp2(breaks=c(-2,0,2),colors=c("blue","white","red")),
        n_threads=100)
dev.off()
dat <- runMetrics(dat)
plt<-ggplot(dat@colData,
        aes(x=paste(cell_line),
            y=dat@colData$overdispersion,
            fill=paste(cell_line),
            color=paste(cell_line)))+
        geom_violin()+geom_jitter()
ggsave(plt,file=paste0(prefix,"_overdispersion.pdf"))
saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(colData(dat),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T)

```

#DCIS_IDC_79

```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)
library(colorRamp2)

cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")

prefix<-"DCIS_IDC_79"
cells_dir<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup"
dcis <- runVarbin(paste(cells_dir,"DCIS-79T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
idc <- runVarbin(paste(cells_dir,"IDC-79T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
colData(dcis)$cell_line<-"DCIS-79T"
colData(idc)$cell_line<-"IDC-79T"
dat<-cbind(dcis,idc)
dat  <- findAneuploidCells(dat  )
dat <- findOutliers(dat)
dat<- runMetrics(dat)
dat<-dat[,dat@colData$overdispersion<0.03]
pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('outlier', 'is_aneuploid'), 
        row_split = 'outlier')
dev.off()
dat<- knnSmooth(dat)
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10
dat <- findClusters(dat, 
                    k_superclones=k_clones@metadata$suggestedK, 
                    k_subclones=k_clones@metadata$suggestedK)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
plotUmap(dat, label = 'superclones')
dev.off()
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('superclones','subclones','reads_total','cell_line'),
        order='hclust',
        col=colorRamp2(breaks=c(-2,0,2),colors=c("blue","white","red")),
        n_threads=100)
dev.off()
dat <- runMetrics(dat)
plt<-ggplot(dat@colData,
        aes(x=paste(cell_line),
            y=dat@colData$overdispersion,
            fill=paste(cell_line),
            color=paste(cell_line)))+
        geom_violin()+geom_jitter()
ggsave(plt,file=paste0(prefix,"_overdispersion.pdf"))

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(colData(dat),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T)


````

## DCIS 92T

```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)
library(colorRamp2)

cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")

prefix<-"DCIS-92T"
cells_dir<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup"
list.files(path=cells_dir,pattern="_split_bam")
dcis <- runVarbin(paste(cells_dir,"DCIS-92T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
colData(dcis)$cell_line<-"DCIS-92T"
dat<-dcis
dat  <- findAneuploidCells(dat  )
dat <- findOutliers(dat)
dat<- runMetrics(dat)
dat<-dat[,dat@colData$overdispersion<0.03]
pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('outlier', 'is_aneuploid'), 
        row_split = 'outlier')
dev.off()
dat<- knnSmooth(dat)
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10
dat <- findClusters(dat, 
                    k_superclones=k_clones@metadata$suggestedK, 
                    k_subclones=k_clones@metadata$suggestedK)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
plotUmap(dat, label = 'superclones')
dev.off()
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('superclones','subclones','reads_total','cell_line'),
        order='hclust',
        col=colorRamp2(breaks=c(-2,0,2),colors=c("blue","white","red")),
        n_threads=100)
dev.off()
dat <- runMetrics(dat)
plt<-ggplot(dat@colData,
        aes(x=paste(cell_line),
            y=dat@colData$overdispersion,
            fill=paste(cell_line),
            color=paste(cell_line)))+
        geom_violin()+geom_jitter()
ggsave(plt,file=paste0(prefix,"_overdispersion.pdf"))
saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(colData(dat),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T)

````

#DCIS-66T

```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)
library(colorRamp2)

cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")

prefix<-"DCIS-66T"
cells_dir<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup"
list.files(path=cells_dir,pattern="_split_bam")
dcis <- runVarbin(paste(cells_dir,"DCIS-66T_split_bam",sep="/"),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)
colData(dcis)$cell_line<-"DCIS-6T"
dat<-dcis
dat  <- findAneuploidCells(dat  )
dat <- findOutliers(dat)
dat<- runMetrics(dat)
dat<-dat[,dat@colData$overdispersion<0.03]
pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('outlier', 'is_aneuploid'), 
        row_split = 'outlier')
dev.off()
dat<- knnSmooth(dat)
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10
dat <- findClusters(dat, 
                    k_superclones=k_clones@metadata$suggestedK, 
                    k_subclones=k_clones@metadata$suggestedK)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
plotUmap(dat, label = 'superclones')
dev.off()
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('superclones','subclones','reads_total','cell_line'),
        order='hclust',
        col=colorRamp2(breaks=c(-2,0,2),colors=c("blue","white","red")),
        n_threads=100)
dev.off()
dat <- runMetrics(dat)
plt<-ggplot(dat@colData,
        aes(x=paste(cell_line),
            y=dat@colData$overdispersion,
            fill=paste(cell_line),
            color=paste(cell_line)))+
        geom_violin()+geom_jitter()
ggsave(plt,file=paste0(prefix,"_overdispersion.pdf"))
saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(colData(dat),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T)

````

#DCIS_41T

```R
library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)
library(ggplot2)
library(colorRamp2)

cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing")

prefix<-"DCIS_41T"
cells_dir<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams"
dcis <- runVarbin(paste(cells_dir),
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)

colData(dcis)$cell_line<-unlist(lapply(strsplit(colData(dcis)$sample,"[.]"),"[",1))
dat<-dcis[,dcis@colData$cell_line %in% c("DCIS-41T")]
dat  <- findAneuploidCells(dat  )
dat <- findOutliers(dat)
dat<- runMetrics(dat)
dat<-dat[,dat@colData$overdispersion<0.03]
pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('outlier', 'is_aneuploid'), 
        row_split = 'outlier')
dev.off()
dat<- knnSmooth(dat)
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10
dat <- findClusters(dat, 
                    k_superclones=k_clones@metadata$suggestedK, 
                    k_subclones=k_clones@metadata$suggestedK)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
plotUmap(dat, label = 'superclones')
dev.off()
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, 
        label = c('superclones','subclones','reads_total','cell_line'),
        order='hclust',
        col=colorRamp2(breaks=c(-2,0,2),colors=c("blue","white","red")),
        n_threads=100)
dev.off()
dat <- runMetrics(dat)
plt<-ggplot(dat@colData,
        aes(x=paste(cell_line),
            y=dat@colData$overdispersion,
            fill=paste(cell_line),
            color=paste(cell_line)))+
        geom_violin()+geom_jitter()
ggsave(plt,file=paste0(prefix,"_overdispersion.pdf"))

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(colData(dat),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T)
````
