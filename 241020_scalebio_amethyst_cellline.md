
# Generate CNV Clones

```bash
singularity shell \
--bind ~/projects \
--bind /volumes/seq/projects/metACT \
--bind /volumes/seq/projects/gccACT \
~/singularity/copykit.sif
source ~/.bashrc
```

```R
library(copykit)
library(BiocParallel)
library(optparse)
library(circlize)
BiocParallel::bpparam()

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="allcells", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=50, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=100
prefix="cellline"

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

#act mdamb231
act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
colData(act)$info <- 'act'
colData(act)$cell_line<-"mda_mb_231"

#scalebio cellline
cellline <- runVarbin("/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/sc_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE,
                  method="multipcf",
                  gamma=20)

#scalebio DCIS round 1
dcis1 <- runVarbin("/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE,
                  method="multipcf",
                  gamma=20)

#scalebio DCIS round 1
dcis2 <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE,
                  method="multipcf",
                  gamma=20)
dat  <- runMetrics(dat)

pdf(paste0(prefix,".qc_metrics.pdf"))
plotMetrics(dat, metric = c("overdispersion", 
                              "breakpoint_count",
                              "reads_total",
                              "reads_duplicates",
                              "reads_assigned_bins",
                              "percentage_duplicates"),
            label = "reads_total")
dev.off()

dat <- dat[,colData(dat)$reads_assigned_bins > 200000 ]
# Mark euploid cells if they exist
dat <- findAneuploidCells(dat,resolution=0.8)

# Mark low-quality cells for filtering
#dat <- findOutliers(dat)
#pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
#plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
#dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat,k=20)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output 

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')
dat <- calcInteger(dat, method = 'scquantum', assay = 'smoothed_bincounts')
pdf("cell_ploidy.pdf")
plotMetrics(dat, metric = 'ploidy', label = 'ploidy_score')
dev.off()

dat <- dat[,colData(dat)$ploidy_score < 0.25]
dat <- calcInteger(dat, method = 'fixed', ploidy_value = 2)
dat <- calcConsensus(dat, consensus_by = 'subclones', assay = 'integer')

# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat,label="subclones")
dev.off()

pdf(paste0(prefix,".subclone.heatmap_int.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones'),
    assay = 'integer',
    order='consensus_tree',
    col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    n_threads=50)
dev.off()

pdf(paste0(prefix,".subclone.heatmap_seg.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones'),
    order='consensus_tree',
    n_threads=50,
    col=colorRamp2(c(-2,-1, 0, 1,2), c("darkblue","blue", "white", "red","darkred")))
dev.off()


pdf(paste0(prefix,".subclone.heatmap_seg.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones'),
    order='consensus_tree',
    n_threads=50,
    col=colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")))
dev.off()


pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".subclone.consensus.pdf"))
plotHeatmap(dat,
            consensus = TRUE,
            assay = 'integer',
            label = 'subclones')
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)
```


# Amethyst run with CNV clones

```bash
singularity shell \
--bind ~/projects/ \
--bind /volumes/seq/projects/metACT \
~/singularity/amethyst.sif
```

```R
library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(plyr)
library(parallel)

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"

setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing")

#read in metadata
samp_dir2=list.files(path=paste0(in_dir2,"/report/"))

met2<-lapply(samp_dir2,function(i){
  samp_met<-read.csv(paste0(in_dir2,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat2<-do.call("rbind",met2)
metadat2$run<-2
metadat2<-metadat2[metadat2$CG_Cov>10000,]
metadat2<-metadat2[metadat2$sampleName %in% c("mcf10a","mcf7","mdamb231"),]
colnames(metadat2)<-tolower(colnames(metadat2))

#rename metadat3 colnames to be consistent with legacy processing (in future just reprocess old samples with new pipeline)
colnames(metadat2)<-c("cell_id","cov","cg_cov","mcg_pct","ch_cov","mch_pct","sample","total_reads","passing_reads","unique_reads","mito_reads","percent","pct_uniq_pass","pct_pass_total","pct_pass_cell","pass","threshold","tgmt","tgmt_well","i5","i5_well","i7","i7_well","cg_total_ratio","run")

#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
obj <- createObject()

metadat2$cov<-metadat2$cg_cov

#filter cells
h5paths2<-paste0(in_dir2,"/h5_files/",metadat2$sample,".",metadat2$tgmt_well,".h5")

obj@h5paths <- data.frame(row.names = rownames(metadat2), paths = c(h5paths2))
obj@metadata<-metadat2 #run from same kit

# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 100) 


#obj<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$sample %in% c("mcf10a"),]))


#split out just MDA cells for clustering
window_name="cg_100k_score";stepsize.=100000;threads.=100;prefix="mcf10a"
obj;bed.=NULL;metric.="score";threads.=100;neighbors.=50
  print(paste("Making window summaries for ",window_name))
  obj@genomeMatrices[[window_name]] <- makeWindows(obj,
                                                      stepsize = stepsize., 
                                                      type = "CG", 
                                                      metric = metric., 
                                                      bed = bed.,
                                                      threads = threads., 
                                                      index = "chr_cg", 
                                                      nmin = 2) 
  print(paste("Estimating dimensions..."))                                           
  #filter windows by cell coverage
  obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
  est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(10), threshold = 0.95)
  print(est_dim)
  set.seed(111)
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim+5, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 80, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    obj <- runTsne(obj, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(prefix,window_name,"_umap.pdf"))     

celltype1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "cluster_id",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_celltype_tracks"]] <- celltype1kbwindows[["pct_matrix"]]

library(rtracklayer) #local
for(i in 4:ncol(obj@genomeMatrices$cg_celltype_tracks)){
cluster=names(obj@genomeMatrices$cg_celltype_tracks)[i]
out_bw<-as.data.frame(obj@genomeMatrices$cg_celltype_tracks)
out_bw<-out_bw[c("chr","start","end",cluster)]
out_bw<-GRanges(out_bw[complete.cases(out_bw),])
names(out_bw@elementMetadata)<-"score"
out_bw<-out_bw[unique(findOverlaps(out_bw, type = "any", select = "first"))]
out_bw <- resize(out_bw, width=1000, fix='start') #resize to avoid 1base overlap
genome(out_bw)<-"hg38"
hg38_seq_info<-Seqinfo(genome="hg38")
seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
print(paste("Saving bigwig for...",cluster))
export(out_bw,con=paste0(cluster,"_mcf10a.bw"),format='bigWig')
}
