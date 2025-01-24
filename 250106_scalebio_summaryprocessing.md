# Reprocessing of all scalebio data generated thus far. 

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
prefix="allcells"

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

#act mdamb231
act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
colData(act)$info <- 'act'
colData(act)$experiment<-"gccact"
colData(act)$sample<-"mda-mb-231"

#scalebio cellline
cellline <- runVarbin("/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/sc_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(cellline)$info <- 'scalebio'
colData(cellline)$experiment<-"cellline_scalebio"
colData(cellline)$sample<-tolower(unlist(lapply(strsplit(colnames(cellline),"[.]"),"[",1)))

#scalebio DCIS round 1
dcis1 <- runVarbin("/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(dcis1)$info <- 'scalebio'
colData(dcis1)$experiment<-"dcis1_scalebio"
colData(dcis1)$sample<-tolower(unlist(lapply(strsplit(colnames(dcis1),"[.]"),"[",1)))

#scalebio DCIS round 2
dcis_66t <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-66T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(dcis_66t)$info <- 'scalebio'
colData(dcis_66t)$experiment<-"dcis2_scalebio"
colData(dcis_66t)$sample<-"dcis_66t"

dcis_79t <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-79T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(dcis_79t)$info <- 'scalebio'
colData(dcis_79t)$experiment<-"dcis2_scalebio"
colData(dcis_79t)$sample<-"dcis_79t"

dcis_92t <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-92T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(dcis_92t)$info <- 'scalebio'
colData(dcis_92t)$experiment<-"dcis2_scalebio"
colData(dcis_92t)$sample<-"dcis_92t"

hbca_17 <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/HBCA-17T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(hbca_17)$info <- 'scalebio'
colData(hbca_17)$experiment<-"dcis2_scalebio"
colData(hbca_17)$sample<-"hbca_17"

hbca_19 <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/HBCA-19T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(hbca_19)$info <- 'scalebio'
colData(hbca_19)$experiment<-"dcis2_scalebio"
colData(hbca_19)$sample<-"hbca_19"

idc_79t <- runVarbin("/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/IDC-79T_split_bam",
                  remove_Y = TRUE, 
                  genome="hg38",
                  is_paired_end=TRUE)
colData(idc_79t)$info <- 'scalebio'
colData(idc_79t)$experiment<-"dcis2_scalebio"
colData(idc_79t)$sample<-"idc_79t"

dat<-cbind(act, cellline, dcis1, 
            dcis_66t, dcis_79t, dcis_92t,
            hbca_17, hbca_19, idc_79t)

dat  <- runMetrics(dat)

saveRDS(dat,"/volumes/USR2/Ryan/projects/metact/cnv_processing/scalebio_allmerged.rds")
setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing/")
prefix="all_samples"

pdf(paste0(prefix,".qc_metrics.pdf"))
plt<-plotMetrics(dat, metric = c("overdispersion", 
                              "reads_assigned_bins"),
            label = "reads_total")

dev.off()

dat <- dat[,colData(dat)$reads_assigned_bins > 250000 ]
colData(dat)<-colData(dat)[row.names(colData(dat)) %in% colnames(dat),]

# Mark euploid cells if they exist
#dat <- findAneuploidCells(dat,resolution=0.8)

# Mark low-quality cells for filtering
#dat <- findOutliers(dat)
#pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
#plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
#dev.off()

#WRITE FUNCTION PER SAMPLE


# kNN smooth profiles
dat <- knnSmooth(dat,k=10)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output 

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')

colData(dat)<-colData(dat)[match(colnames(dat@assays@data$smoothed_bincounts),row.names(colData(dat))),] #reorder metadata
dat <- calcInteger(dat)

#pdf("cell_ploidy.pdf")
#plotMetrics(dat, metric = 'ploidy', label = 'ploidy_score')
#dev.off()

#dat <- dat[,colData(dat)$ploidy_score < 0.25]
#dat <- calcInteger(dat, method = 'fixed', ploidy_value = 2)
dat <- calcConsensus(dat, consensus_by = 'subclones', assay = 'integer')

# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat,label="subclones")
plotUmap(dat,label="sample")
dev.off()

pdf(paste0(prefix,".subclone.heatmap_segratio.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones','sample','experiment'),
    order='hclust',
    col=colorRamp2(c(-2,-1, 0, 1,2), c("darkblue","blue", "white", "red","darkred")),
    n_threads=50)
dev.off()


pdf(paste0(prefix,".subclone.heatmap_smoothedbin.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones','sample','experiment'),
    order='hclust',assay = 'smoothed_bincounts',
    col=colorRamp2(c(-2,-1, 0, 1,2), c("darkblue","blue", "white", "red","darkred")),
    n_threads=50)
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".subclone.integer.pdf"))
plotHeatmap(dat, label = c('reads_total','subclones','sample','experiment'),
            consensus = FALSE,
            assay = 'integer')
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
#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
dcis_cnv<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams/all_cells.scCNA.tsv"
in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
in_dir3="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst"

setwd(in_dir3)
#obj<-readRDS("hbca_dcis.met.Rds")

#read in metadata
samp_dir=list.files(path=paste0(in_dir,"/report/"))
samp_dir2=list.files(path=paste0(in_dir2,"/report/"))
samp_dir3=list.files(path="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples",pattern="*csv",full.names=T)

met1<-lapply(samp_dir,function(i){
  samp_met<-read.csv(paste0(in_dir,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat1<-do.call("rbind",met1)
metadat1$run<-1
metadat1<-metadat1[metadat1$CG_Cov>10000,]
metadat1<-metadat1[metadat1$sampleName %in% c("hbca-83l","hbca-16r","DCIS-41T","DCIS-66T"),]
colnames(metadat1)<-tolower(colnames(metadat1))

met2<-lapply(samp_dir2,function(i){
  samp_met<-read.csv(paste0(in_dir2,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat2<-do.call("rbind",met2)
metadat2$run<-2
metadat2<-metadat2[metadat2$CG_Cov>10000,]
metadat2<-metadat2[metadat2$sampleName %in% c("hbca-83l" ,"hbca-16r" ,"DCIS-41T", "DCIS-66T"),]
colnames(metadat2)<-tolower(colnames(metadat2))

met3<-lapply(samp_dir3,function(i){
  samp_met<-read.csv(i)
  row.names(samp_met)<-samp_met$cell_id
  return(samp_met)
})
metadat3<-do.call("rbind",met3)
metadat3$run<-3
metadat3<-metadat3[metadat3$cg_cov>10000,]
metadat3<-metadat3[complete.cases(metadat3),]

#rename metadat3 colnames to be consistent with legacy processing (in future just reprocess old samples with new pipeline)
colnames(metadat1)<-c("cell_id","cov","cg_cov","mcg_pct","ch_cov","mch_pct","sample","total_reads","passing_reads","unique_reads","mito_reads","percent","pct_uniq_pass","pct_pass_total","pct_pass_cell","pass","threshold","tgmt","tgmt_well","i5","i5_well","i7","i7_well","cg_total_ratio","run")
colnames(metadat2)<-c("cell_id","cov","cg_cov","mcg_pct","ch_cov","mch_pct","sample","total_reads","passing_reads","unique_reads","mito_reads","percent","pct_uniq_pass","pct_pass_total","pct_pass_cell","pass","threshold","tgmt","tgmt_well","i5","i5_well","i7","i7_well","cg_total_ratio","run")
metadat1<-metadat1[colnames(metadat1) %in% colnames(metadat3)]
metadat2<-metadat2[colnames(metadat2) %in% colnames(metadat3)]
metadat3<-metadat3[colnames(metadat3) %in% colnames(metadat2)]

metadat3<-metadat3[complete.cases(metadat3),]

#####NOT RUN#####
#Add CNV metadata for DCIS
#cnv_dat<-read.table(dcis_cnv,sep="\t",header=T)
#cnv_dat$idx<-gsub("_","+",cnv_dat$idx)
#metadat$cnv_clone<-"diploid"
#metadat[match(cnv_dat$idx,row.names(metadat)),]$cnv_clone<-cnv_dat$superclones

#head(obj@metadata)
#obj@metadata<-obj@metadata[obj@metadata$CG_Cov>10000,]
#plt<-ggplot(obj@metadata, aes(x=sampleName, y = CG_Cov)) +geom_violin() + geom_jitter()
#ggsave(file="cov_plot.pdf",plt)
#################

#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
obj <- createObject()

metadat1$cov<-metadat1$cg_cov
metadat2$cov<-metadat2$cg_cov
metadat3$cov<-metadat3$cg_cov
#row.names(metadat1)<-paste(row.names(metadat1),metadat1$run,sep="_")
#row.names(metadat2)<-paste(row.names(metadat2),metadat2$run,sep="_")
row.names(metadat3)<-paste(row.names(metadat3),metadat3$run,sep="_")

#filter cells
#obj@metadata <- obj@metadata |> dplyr::filter(cov > 100000 & cov < 40000000)
h5paths<-paste0(in_dir,"/cg_sort_cov/h5_files/",metadat1$sample,".",metadat1$tgmt_well,".h5")
h5paths2<-paste0(in_dir2,"/h5_files/",metadat2$sample,".",metadat2$tgmt_well,".h5")
h5paths3<-paste0(in_dir3,"/",metadat3$sample,"/",metadat3$sample,".",metadat3$tgmt_well,"_cov.h5")

obj@h5paths <- data.frame(row.names = rownames(rbind(metadat1,metadat2,metadat3)), paths = c(h5paths,h5paths2,h5paths3))
obj@metadata<-rbind(metadat1,metadat2,metadat3) #run from same kit

#RUN ONCE
#do a for loop to write the corrected cell name into the h5 files
#excluding CH methylation for now
#correct_h5_cellnames<-function(h5){
#    print(paste("Correcting names for...", basename(h5)))
#    run_number=as.character(3)
#    h5list = h5ls(h5)
#    for (i in 1:nrow(h5list)){
#        tryCatch({if(endsWith(h5list[i,"group"],"CG")){
#            if( !(paste0(h5list[i,"name"],"_",run_number) %in% h5list$name) &&
#              !endsWith(h5list[i,"name"],run_number)){
#                celldat<-h5read(h5,paste0("CG/",h5list[i,"name"]))
#                h5write(celldat, file=h5, name=paste0("CG/",h5list[i,"name"],"_",run_number))
#            }
#        }
#        },error =function(e) { cat("Proceeding past line",i,"for",basename(h5),"\n")} )
#    }
#    }

#mclapply(unique(h5paths3),correct_h5_cellnames,mc.cores=50)
#h5closeAll()


# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 200) 


cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=200,neighbors.=50){
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
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 175, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = 100, dist = 0.2, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    

for( j in c(0.05,0.1,0.2,0.3)){
plt_list<-lapply(c(10,25,50,75,100,150), function(i){
        obj <- runUmap(obj, neighbors = as.numeric(i), dist = as.numeric(j), 
          method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
        p <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(i,j))
        return(p)})
plt<-patchwork::wrap_plots(plt_list)
ggsave(plt,file=paste0(window_name,"_umap.",j,".pdf")) 
}
 

    obj <- runTsne(obj, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj<-cluster_by_windows(obj,window_name="cg_100k_score",stepsize.=100000,threads.=200)

#rerun clustering (overcluster and then merge them visually)
window_name="cg_100k_score"
obj <- runCluster(obj, k_phenograph = 175, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
obj@metadata$cluster_broad<-obj@metadata$cluster_id
obj <- runCluster(obj, k_phenograph = 25, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
obj@metadata$cluster_fine<-obj@metadata$cluster_id
obj@metadata$cluster_id<-obj@metadata$cluster_broad

obj <- runTsne(obj, perplexity=50, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
obj <- runUmap(obj, neighbors = 50, dist = 0.0001, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 

print("Plotting...")
p1 <- dimFeature(obj, colorBy = cluster_broad, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
p5 <- dimFeature(obj, colorBy = cluster_fine, reduction = "tsne") + ggtitle(paste(window_name,"TSNE Clusters"))
p6 <- dimFeature(obj, colorBy = cluster_broad, reduction = "tsne") + ggtitle(paste(window_name,"TSNE Clusters Broad"))
plt<-plot_grid(p1, p2,p3, p4, p5,p6,ncol=2)
ggsave(plt,file=paste0(window_name,"_umap.pdf"),height=10)    

table(obj@metadata$cluster_id)
table(obj@metadata$sample)

#protein coding gene promoters
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
cgi<-rtracklayer::import(BEDFile("/volumes/USR2/Ryan/projects/metact/ref/cgi_hg38.bed.gz"))

for (i in c("gene_name", "exon_number")) {gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
obj@ref<-gtf
protein_coding <- unique(obj@ref |> 
                dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> 
                dplyr::pull(gene_name))


#promoter extension and cgi filter
extend <- function(x, upstream=0, downstream=0)     
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

protein_coding_gtf<-gtf[gtf$type=="gene" & gtf$gene_type=="protein_coding" & gtf$seqid!="chrM",]
protein_coding_gtf<-protein_coding_gtf[!duplicated(protein_coding_gtf$gene_name),]
protein_coding_gtf<-makeGRangesFromDataFrame(protein_coding_gtf,keep.extra.columns=TRUE)
protein_coding_gtf<-extend(protein_coding_gtf,downstream=20000) #extend into gene body
protein_coding_gtf$cgi_overlap<-!is.na(findOverlaps(protein_coding_gtf,cgi,ignore.strand=TRUE,type="any",select="first")) #filter out promoters with CGI
promoter_bed<-data.frame(seqnames=seqnames(protein_coding_gtf),start=start(protein_coding_gtf)-1,end=end(protein_coding_gtf))
obj@genomeMatrices[["cg_promoter"]] <- makeWindows(obj, 
                                                     bed = promoter_bed, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 



obj@genomeMatrices[["cg_genebody"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

cluster_by_windows<-function(obj,window_name=cg_genebody,est_dim=5,threads.=100,neighbors.=50){
  print(est_dim)
  set.seed(111)
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 50, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(window_name,"_umap.pdf"))     
  return(obj)                               
}

#############Gene body methylation#############
dimEstimate(obj, genomeMatrices = c("cg_promoter"), dims = c(10), threshold = 0.95)

cluster_by_windows(obj,est_dim=5,window_name="cg_promoter")

#cluster_by_windows(obj,est_dim=3, window_name="cg_genebody")
saveRDS(obj,file="hbca_dcis.met.Rds")


###### DMR Calculation ####

cluster1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "cluster_id",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_cluster_tracks"]] <- cluster1kbwindows[["pct_matrix"]]

clusterfine1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "cluster_fine",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_clusterfine_tracks"]] <- clusterfine1kbwindows[["pct_matrix"]]


saveRDS(obj,file="hbca_dcis.met.Rds")

pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
dmrs<-testDMR(clusterfine1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = FALSE) #add additional columns direction column
celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c","",colnames(clusterfine1kbwindows[["sum_matrix"]])[grepl(pattern="_c",colnames(clusterfine1kbwindows[["sum_matrix"]]))]))
dmrs2$celltype<-celltype_test[dmrs2$test]
save(dmrs2,file="hbca_dcis.celltype_fine.dmr.Rds")

collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 1000, reduce = T, annotate = T) 
collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = makePalette(option = 7, n = 13) ) + theme_classic()
ggsave(plt,file="met_per_dmr_celltype.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(celltype, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(celltype, direction) |> slice_min(n = 20, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)
write.table(top_dmrs,file="cluster_dmr_fine.tsv")

plt<-heatMap(obj, matrix = "cg_celltype_tracks", regions = top_dmrs$location, nrow = 10, legend = FALSE, width = 1000, arrowOverhang = 2000)
ggsave(plt,file="dmrs_heatmap_celltype.pdf",width=50,height=50,limitsize=FALSE)


#subset 9 4 8 for immune
#arbitrary set from fine clusters just for plotting poster

celltype=c(
"1"="lumHR1","3"="lumHR1","12"="lumHR1","15"="lumHR1",
"6"="lumHR2",
"14"="lumHR3",
"17"="lumHR4","19"="lumHR4","23"="lumHR4",
"21"="lumHR5",
"22"="lumHR6",
"24"="lumHR7",
"9"="basal1",
"16"="lumsec",

"2"="perivasc",
"5"="fibro1",
"7"="endo",
"10"="fibro2",
"11"="stroma4",

"4"="tnkcell",
"8"="myeloid2",
"13"="myeloid1",
"18"="mast",
"20"="bcell")


obj@metadata$celltype<-celltype[obj@metadata$cluster_fine]

celltype1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "celltype",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_celltype_tracks"]] <- celltype1kbwindows[["pct_matrix"]]


saveRDS(obj,file="hbca_dcis.met.Rds")


###############################################################################################################
#################################################BIGWIG OUTPUT#################################################
###############################################################################################################
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
library(rtracklayer) #local

in_dir3="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst"

setwd(in_dir3)
#obj<-readRDS(file="hbca_dcis.met.Rds")

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
export(out_bw,con=paste0(cluster,"_celltype.bw"),format='bigWig')
}

###############################################################################################################
###############################################################################################################
###############################################################################################################

###############################################################################################################
########################################HYPOMETHYLATION PLOT###################################################
###############################################################################################################




colors=c(
"lumHR1"="#e33c96",
"lumHR2"="#ad2286",
"lumHR3"="#c03a95",
"lumHR4"="#ed3d94",
"lumHR5"="#bc70ae",
"lumHR6"="#7b2881",
"lumHR7"="#ae62a7",
"basal1"="#7f479c",
"lumsec"="#5cbb5b",

"perivasc"="#fed700",
"fibro1"="#ed1d7e",
"fibro2"="#b51f61",
"endo"="#a1d39d",
"stroma4"="#35b551",

"tnkcell"="#489bd5",
"myeloid2"="#70c6aa",
"myeloid1"="#35b551",
"mast"="#c4aed4",
"bcell"="#3953a4")

genes=c("KIT","PTPRC","EPCAM","KRT14","COL1A1","ELF5")

library(tidyverse)
library(colorspace)

cgisland<-"/volumes/USR2/Ryan/projects/metact/ref/cpgIslandExt.bed"

#set colors for each facet and use greyscaling for higher met levels
#https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
#histograM <- function(obj,

              genes=c("COL1A1")       
                      matrix = "cg_celltype_tracks";
                      #colors = NULL;
                      trackOverhang = 5000;
                      arrowOverhang = 3000;
                      ncol = length(genes);
                      invert_met = TRUE;
                      legend = TRUE;
                      removeNA = TRUE;
                      width = 1000;
                      trackScale = 1.5;
                      colorMax = 100;#) {

  if (!is.null(colors)) {
    pal <- colors
  } else {
    pal <- rev(c("#FF0082", "#dbdbdb", "#cccccc", "#999999"))
  }

  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }

  p <- vector("list", length(genes)) # empty plot list
  for (i in 1:length(genes)) {
    ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
    aggregated <- obj@genomeMatrices[[matrix]]

    toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                              aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                              aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
    ngroups <- ncol(toplot) - 3
    trackHeight <- ngroups * trackScale
    toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
    if (removeNA) {
      toplot <- toplot |> dplyr::filter(!is.na(pct_m))
    }
    if (invert_met){
      toplot$pct_m<-toplot$pct_m-100
      colorMax=-100
    }
    #initialize first row
    p[[i]] <- vector("list", length(colors)) # empty plot list
    
    for (j in names(colors)){
    toplot_sub<-toplot[toplot$group==j,]
    p[[i]][[j]] <- ggplot2::ggplot() + 
      ggplot2::geom_col(data = toplot_sub, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = width) + 
      scale_fill_gradient2(high="black",mid="grey",low=colors[j],midpoint=-60) + 
      theme_void() +
      theme(legend.position="none") + ylab(j)
    if(j==names(colors)[1]){
    p[[i]][[j]] <- p[[i]][[j]] + 
      ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |>
                          dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)),
                          promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                          ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = trackHeight*3, ymax = trackHeight)) +
      ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"),
                         ggplot2::aes(xmin = start, xmax = end, ymin = trackHeight*3, ymax = trackHeight)) +
      ggplot2::geom_segment(data = ref, aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang),
                          y = (trackHeight*1.5),
                          xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang),
                          yend = (trackHeight*1.5)), arrow = arrow(length = unit(trackHeight/40, "cm"))) + 
                          xlab(paste(genes[i],toplot$chr[1],toplot$start[1],"-",toplot$end[nrow(toplot)])) +
                          theme_void() + theme(legend.position="none")
      p[[i]][[j]] <- p[[i]][[j]]+ 
      ggplot2::scale_fill_gradientn(colors = rev(pal), limits = c(colorMax,0)) + 
      theme(axis.title.y = element_blank(),axis.text.y=element_blank()) +
      ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank()) + ylab(j) + ggtitle(genes[i]) + theme(axis.title.y = element_blank(),axis.text.y=element_blank(),legend.position="none")
  }
  if(j==last(names(colors))){
    p[[i]][[j]] <- p[[i]][[j]] + theme_void() + theme(axis.title.y = element_blank(),axis.text.y=element_blank(),legend.position="none")
 #to put genome location
  }
}
  ggsave(gridExtra::grid.arrange(grobs = p[[i]], ncol = ncol),file="test.pdf",limitsize=F,height=10)

}


###############################################################################################################
###############################################################################################################
###############################################################################################################

```bash
singularity shell \
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
library(rtracklayer) #local

setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing")
obj<-readRDS(file="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst/hbca_dcis.met.Rds")

#CHANGE THIS TO ALL CNVS
cnv<-read.table("/volumes/USR2/Ryan/projects/metact/cnv_processing/DCIS_41T.scCNA.tsv",sep="\t",header=T)
cnv$cellid<-gsub("_","+",unlist(lapply(strsplit(row.names(cnv),"[.]"),"[",3)))

obj@metadata$cnv<-"NA"
row.names(cnv)<-cnv$cellid
obj@metadata[row.names(obj@metadata) %in% row.names(cnv),]$cnv<-cnv$subclones
obj_cnv<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$sample %in% c("DCIS-41T"),]))

cluster_by_windows<-function(obj,window_name,stepsize.=NULL,
bed.=NULL,metric.="score",threads.=200,neighbors.=50,prefix="DCIS_41T"){
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
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 175, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    obj <- runTsne(obj, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(prefix,".",window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj_cnv<-cluster_by_windows(obj_cnv,window_name="cg_100k_score",stepsize.=100000,threads.=200,prefix="DCIS_41T")

p1 <- dimFeature(obj_cnv, colorBy = cluster_id, reduction = "umap") + ggtitle(paste("Clusters"))
p2 <- dimFeature(obj_cnv, colorBy = cnv, reduction = "umap") + ggtitle(paste("Samples"))
plt<-plot_grid(p1, p2,ncol=2)

ggsave(plt,file="DCIS_41T_umap.pdf")    






setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing")
obj<-readRDS(file="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst/hbca_dcis.met.Rds")

cnv<-read.table("/volumes/USR2/Ryan/projects/metact/cnv_processing/DCIS_IDC_79.scCNA.tsv",sep="\t",header=T)
cnv$cellid<-gsub("_","+",unlist(lapply(strsplit(row.names(cnv),"[.]"),"[",3)))

obj@metadata$cnv<-"NA"
row.names(cnv)<-paste0(cnv$cellid,"_3")
obj@metadata[row.names(obj@metadata) %in% row.names(cnv),]$cnv<-cnv$subclones
obj_cnv<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$sample %in% c("DCIS-79T","IDC-79T"),]))


cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=200,neighbors.=50,prefix="DCIS_41T"){
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
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 175, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    obj <- runTsne(obj, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(prefix,".",window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj_cnv<-cluster_by_windows(obj_cnv,window_name="cg_100k_score",stepsize.=100000,threads.=200,prefix="DCIS_IDC_79T")

p1 <- dimFeature(obj_cnv, colorBy = cluster_id, reduction = "umap") + ggtitle(paste("Clusters"))
p2 <- dimFeature(obj_cnv, colorBy = cnv, reduction = "umap") + ggtitle(paste("Samples"))
plt<-plot_grid(p1, p2,ncol=2)

ggsave(plt,file="DCIS_IDC_79T_umap.pdf")    

library(reshape2)
obj@metadata[startsWith(obj@metadata$celltype,prefix="lumHR"),]$celltype<-"lumHR"
out<-table(obj@metadata$celltype,obj@metadata$sample)
out<-melt(out)
colnames(out)<-c("celltype","sample","value")
plt<-ggplot(out, aes(fill=celltype, y=value,x=sample)) + 
    geom_bar(position="fill", stat="identity")
ggsave(plt,file="celltype_stackedbar.pdf")    
