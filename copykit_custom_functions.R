#singularity shell \
#--bind ~/projects \
#--bind /volumes/seq/projects/metACT \
#--bind /volumes/seq/projects/gccACT \
#~/singularity/copykit.sif
#source ~/.bashrc

library(copykit)
library(BiocParallel)
library(optparse)
library(circlize)
BiocParallel::bpparam()

#source("~/projects/metact/src/copykit_custom_functions.R") to load in


generate_copykit_object<-function(indir="/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                                    paired_end=TRUE, info="scalebio",
                                    experiment_name="gccact", sample="patient"){
    dat <- runVarbin(indir,remove_Y = TRUE,genome="hg38",is_paired_end=paired_end)
    colData(dat)$info <- info
    colData(dat)$experiment<-experiment_name
    colData(dat)$sample_name<-sample #copykit uses metadata column named sample
    return(dat)
}


process_copykit_object<-function(dat=cellline,prefix="scalebio_cellline",min_readcount=250000){
    dat <- dat[,colData(dat)$reads_assigned_bins > min_readcount ]
    dat  <- runMetrics(dat)
    dat <- findAneuploidCells(dat) #mark aneuploid cells, i dont take much stock in this
    dat <- knnSmooth(dat,k=10) # kNN smooth profiles
    dat <- runUmap(dat) # Create a umap embedding 
    k_clones<-findSuggestedK(dat) 
    dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3) #use suggested k for clustering
    dat <- calcInteger(dat, method = 'scquantum',assay = 'segment_ratios')

    pdf(paste0(prefix,".qc_metrics.pdf"))
    plt<-plotMetrics(dat, metric = c("overdispersion", 
                                    "reads_assigned_bins",
                                    "ploidy"),
                label = "reads_total")
    print(plt)
    dev.off()
    dat <- calcConsensus(dat) # Calculate consensus profiles for each subclone, 
    dat <- runConsensusPhylo(dat) # and order cells by cluster for visualization with plotHeatmap
    dat <- runPhylo(dat, metric = 'manhattan')

    dat <- calcConsensus(dat, consensus_by = 'subclones', assay = 'integer')

    # Plot a copy number heatmap with clustering annotation
    pdf(paste0(prefix,".subclone.umap.pdf"))
    plt1<-plotUmap(dat,label="subclones")
    plt2<-plotUmap(dat,label='sample_name')
    print(plt1)
    print(plt2)
    dev.off()

    #plot segmentation ratio
    pdf(paste0(prefix,".subclone.heatmap_segratio.pdf"))
    plt<-plotHeatmap(dat, label = c('reads_total','subclones','sample_name','experiment'),
                    order_cells = 'consensus_tree',
                    col=colorRamp2(c(-2,-1, 0, 1,2), c("darkblue","blue", "white", "red","darkred")),
                    n_threads=50)
    print(plt)
    dev.off()

    #plot smoothed bins
    pdf(paste0(prefix,".subclone.heatmap_smoothedbin.pdf"))
    plt<-plotHeatmap(dat, label = c('reads_total','subclones','sample_name','experiment'),
                    assay = 'smoothed_bincounts',
                    order_cells = 'consensus_tree',
                    n_threads=50)
    print(plt)
    dev.off()

    #plot integer
    pdf(paste0(prefix,".subclone.integer.pdf"))
    plt<-plotHeatmap(dat, label = c('reads_total','subclones','sample_name','experiment'),
                consensus = FALSE,
                order_cells = 'consensus_tree',
                assay = 'integer')
    print(plt)
    dev.off()

    

    saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
    write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)
    return(dat)
}
