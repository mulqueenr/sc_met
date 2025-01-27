
# Amethyst run with CNV clones

```bash
singularity shell \
--bind ~/projects/ \
--bind /volumes/seq/projects/metACT \
~/singularity/amethyst.sif
```

```R
source("/volumes/USR2/Ryan/projects/metact/src/amethyst_custom_functions.R") #this loads all necessary R packages and custom functions
setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
indir1<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
indir2<-"/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
indir3<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples"

in_dir3="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst"
```
### Read in metadata for all samples processed so far (excluding cell lines here)

```R
met1<-setup_amethyst_metadata(in_dir=indir1,sample_list=c("hbca-83l","hbca-16r","DCIS-41T","DCIS-66T"),cg_cov_filter=50000)
met2<-setup_amethyst_metadata(in_dir=indir2, sample_list=c("HBCA-83L" ,"HBCA-16R" ,"DCIS-41T", "DCIS-66T"),runid="2",cg_cov_filter=50000)
met3<-setup_amethyst_metadata(in_dir=indir3, sample_list=c("DCIS-66T","DCIS-79T","DCIS-92T", "HBCA-17T", "HBCA-19T" , "IDC-79T") ,runid="3",cg_cov_filter=50000)

met1<-met1[colnames(met1) %in% colnames(met3)]
met2<-met2[colnames(met2) %in% colnames(met3)]
met3<-met3[colnames(met3) %in% colnames(met2)]
met3<-met3[complete.cases(met3),]

#################### Read in h5paths based on metadata names #################### 
#filter cells
h5paths1<-paste0(indir1,"/cg_sort_cov/h5_files/",met1$sample,".",met1$tgmt_well,".h5")
h5paths2<-paste0(indir2,"/h5_files/",met2$sample,".",met2$tgmt_well,".h5")
h5paths3<-paste0(indir3,"/",met3$sample,"/",met3$sample,".",met3$tgmt_well,"_cov.h5")
h5paths3<-paste0(indir3,"/methylation_coverage/amethyst/",met3$sample,"/",met3$sample,".",met3$tgmt_well,"_cov.h5")

#################### #################### Set up object #################### #################### 
#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
obj <- createObject()
row.names(met3)<-paste(row.names(met3),met3$run,sep="_")#this is so the same index can be used when multiple kits are run, no need here because we only have two kits

#do a for loop to write the corrected cell name into the h5 files
#excluding CH methylation for now
#mclapply(unique(h5paths3),function(x) {correct_h5_cellnames(h5=x,runid=3)},mc.cores=50)
#h5closeAll()

obj@h5paths <- data.frame(row.names = rownames(rbind(met1,met2,met3)), paths = c(h5paths1,h5paths2,h5paths3)) 
obj@metadata<-rbind(met1,met2,met3) 
obj@metadata <- obj@metadata |> dplyr::filter(cov > 50000 & cov < 40000000)
obj@metadata$cov<-obj@metadata$cg_cov

# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 200) 
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
obj@ref<-gtf

```

### HBCA ONLY
```R

#################### Initial Clustering 100kb and 50kb #################### 
hbca <- subsetObject(obj, cells = row.names(obj@metadata[startsWith(toupper(obj@metadata$sample),prefix="HBCA"),]))

#decided on 50kb given the data quality
hbca<-cluster_by_windows(obj=hbca,
                        prefix="hbca",
                        window_name="cg_100k_score",
                        stepsize=50000,
                        threads=200,
                        neighbors=20,
                        distance=0.05,
                        overwrite_windows=FALSE,
                        k_phenograph=25)

saveRDS(hbca,"hbca.amethyst.rds")


#################### DMR between clusters and 1kb for cell typing 
 
# Run 100kb windows, cluster, run 1kb windows on clusters, find DMRs, recluster on DMRs, output bigwig tracks for ID

hbca<-dmr_and_1kb_window_gen(obj=hbca,prefix="hbca")

#rerun umap with dmrs
dmr<-readRDS(file="hbca.dmr.cluster_id.collapsed.rds")
dmr_bed<-data.frame(seqnames=dmr$chr,start=dmr$dmr_start-1,end=dmr$dmr_end)

hbca<-cluster_by_windows(obj=hbca,
                        window_name="cg_cluster_dmrs",
                        prefix="hbca",
                        bed=dmr_bed,
                        threads=200,
                        neighbors=20,#20
                        k_phenograph=60,#70
                        distance=0.05 #0.05
                        )
#summarize to 1kb windows per cell type
cluster1kbwindows <- calcSmoothedWindows(hbca, 
                                        type = "CG", 
                                        threads = 300,
                                        step = 1000,
                                        smooth = 3,
                                        index = "chr_cg",
                                        groupBy = "cluster_id", 
                                        returnPctMatrix = TRUE)
hbca@genomeMatrices[["cg_cluster_dmrs_tracks"]] <- cluster1kbwindows[["pct_matrix"]]
bigwig_output(obj=hbca,tracks="cg_cluster_dmrs_tracks")

saveRDS(hbca,"hbca.amethyst.rds")

```
### Additional Clustering With DCIS

```R

table(toupper(obj@metadata$sample))
#decided on 50kb given the data quality
obj<-cluster_by_windows(obj=obj,
                        prefix="dcis_hbca",
                        window_name="cg_100k_score",
                        stepsize=100000,
                        threads=300,
                        neighbors=20,
                        distance=0.05,
                        overwrite_windows=FALSE,
                        k_phenograph=50)

#read and filter dmrs
obj<-dmr_and_1kb_window_gen(obj=obj,prefix="dcis_hbca")

bigwig_output(obj=obj,tracks="cg_cluster_id_tracks")


#dmrs2<-readRDS(paste0(prefix,".dmr.",groupBy,".rds"))
#dmrs2<-filterDMR(dmrs2, method = "bonferroni", filter = FALSE) #add additional columns direction column
#dmrs2<-dmrs2[dmrs2$pval<0.01 & dmrs2$direction=="hypo",]
#dmrs2<-makeGRangesFromDataFrame(dmrs2,keep.extra.columns=TRUE)
#dmrs2<-reduce(GRanges(dmrs2))
#dmrs2<-dmrs2[width(dmrs2)<50000,]

#164731 dmrs
#median size 5kb, mak is 49kb
#dmr_bed<-data.frame(seqnames=seqnames(dmrs2),start(dmrs2),end=end(dmrs2))


#obj<-cluster_by_windows(obj=obj,
#                        window_name="cg_cluster_dmrs",
#                        prefix="dcis_hbca",
#                        bed=dmr_bed,
#                        threads=200,
#                        neighbors=20,#20
#                        k_phenograph=60,#70
#                        distance=0.05 #0.05
#                        )

saveRDS(obj,"dcis_hbca.amethyst.rds")

```

```R
source("~/projects/metact/src/amethyst_custom_functions.R")
setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
obj<-readRDS("dcis_hbca.amethyst.rds")

#Read in CNV clones for clonal analysis
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
obj@ref<-gtf
cnv<-list.files(path="/volumes/USR2/Ryan/projects/metact/cnv_processing",pattern="*.scCNA.tsv") #set wd
obj@metadata$subclones="NULL"

#41
dcis_41<-read.table("/volumes/USR2/Ryan/projects/metact/cnv_processing/dcis-41t.scCNA.tsv")
row.names(dcis_41)<-gsub(pattern="_",replacement="+",x=unlist(lapply(strsplit(row.names(dcis_41),"[.]"),"[",3)))
clone<-setNames(nm=row.names(dcis_41),paste0("dcis_41t_",dcis_41$subclone))
obj@metadata[names(clone),]$subclones<-clone

#66
dcis_66<-read.table("/volumes/USR2/Ryan/projects/metact/cnv_processing/scalebio_dcis66t.scCNA.tsv")
row.names(dcis_66)<-gsub(pattern="_",replacement="+",x=unlist(lapply(strsplit(row.names(dcis_66),"[.]"),"[",3)))
row.names(dcis_66)<-paste0(row.names(dcis_66),"_3")
clone<-setNames(nm=row.names(dcis_66),paste0("dcis_66t_",dcis_66$subclone))
obj@metadata[names(clone),]$subclones<-clone

#79
dcis_79<-read.table("/volumes/USR2/Ryan/projects/metact/cnv_processing/scalebio_idcdcis79t.scCNA.tsv")
row.names(dcis_79)<-gsub(pattern="_",replacement="+",x=unlist(lapply(strsplit(row.names(dcis_79),"[.]"),"[",3)))
row.names(dcis_79)<-paste0(row.names(dcis_79),"_3")
clone<-setNames(nm=row.names(dcis_79),paste0("dcis_79t_",dcis_79$subclone))
obj@metadata[names(clone),]$subclones<-clone

table(obj@metadata[,c("subclones","sample")])
saveRDS(obj,"dcis_hbca.amethyst.rds")

```


## Clone analysis per sample

```R
#################### DCIS CLONE ANALYSIS 
source("~/projects/metact/src/amethyst_custom_functions.R")
setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
obj<-readRDS("dcis_hbca.amethyst.rds")

clone_dmr(obj=obj,prefix="DCIS-41T",sample="DCIS-41T",k_phenograph=50)
clone_dmr(obj=obj,prefix="DCIS-66T",sample="DCIS-66T",k_phenograph=50)
clone_dmr(obj=obj,prefix="IDC_DCIS-79T",sample=c("IDC-79T","DCIS-79T"),k_phenograph=50)

```

## Methyltree per clone, made a custom function to mirror inputs to methyltree and do a per-sample window selection based on the paper's method

```R
source("/volumes/USR2/Ryan/projects/metact/src/amethyst_custom_functions.R") #this loads all necessary R packages and custom functions
setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
obj<-readRDS("dcis_hbca.amethyst.rds")

methyltree_output(obj=obj,prefix="DCIS-41T",sample="DCIS-41T",filt_min_pct=10,filt_max_pct=80,threads=100)
methyltree_output(obj=obj,prefix="DCIS-66T",sample="DCIS-66T",filt_min_pct=10,filt_max_pct=80,threads=100)
methyltree_output(obj=obj,prefix="IDC_DCIS-79T",sample=c("IDC-79T","DCIS-79T"),filt_min_pct=10,filt_max_pct=80,threads=100)

```


# Meme suite
# Want FIMO or AME i think
# DMRs (small peaks) per clsuter
# get_sequence()
# AME or STREME/XSTREME
# Hypergeo test?

```R
source("/volumes/USR2/Ryan/projects/metact/src/amethyst_custom_functions.R") #this loads all necessary R packages and custom functions

options(future.globals.maxSize = 20000 * 1024^2) #set parallelization max size to 20GB

setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
obj<-readRDS(file="dcis_hbca.celltypepeaks.amethyst.rds")
threads=50
prefix="allcells"

#chromvar preparation for all cells
chromvar_met_per_cell<-function(obj,stepsize=500,threads=50,percent_cell_coverage=2.5){
  #run window scoring for all cells
  chromvar_windows <- makeWindows(obj,
                                  stepsize = stepsize, 
                                  type = "CG", 
                                  metric = "score", 
                                  threads = threads, 
                                  index = "chr_cg", 
                                  nmin = 2) 
  #binarized based on hypomethylation score (1 for hypo, 0 for hyper)
  counts<-ifelse(chromvar_windows<0.25,1,0) #score ranges from -1 to 1, using 0.25 for cutoff
  #require 2.5% cell coverage for windows, this is kinda on par with ATAC data, kinda an arbitrary cutoff. but mostly to decrease computational time on which windows we scan for motifs (5% is 10k window, 2% is 909k windows)
  counts<-counts[rowSums(!is.na(counts))>=(ncol(counts)/100)*percent_cell_coverage,] 
  counts[is.na(counts)]<-0
  #remove purely hypermethylated windows
  counts<-counts[rowSums(counts)>1,] 
  #dim(counts)
  return(counts)}


#run on the cluster level for less coverage noise
chromvar_met_per_cluster<-function(obj,stepsize=500,threads=50,mincov=2,percent_cell_coverage=5,groupBy="cluster_id"){
  #run window scoring for all cells
  clusterwindows <- calcSmoothedWindows(obj, 
              type = "CG", 
              threads = threads,
              step = 500,
              smooth = 1,
              index = "chr_cg",
              groupBy = groupBy, 
              returnSumMatrix = TRUE, 
              returnPctMatrix = TRUE)
  cov_mat<-clusterwindows[["sum_matrix"]]
  cov_mat<-cov_mat[,grepl(colnames(cov_mat,pattern="_t$"))]
  pct_mat<-clusterwindows[["pct_matrix"]]
  pct_mat[cov_mat<=min_cov]<-NA
  #binarized based on hypomethylation score (1 for hypo, 0 for hyper)
  counts<-ifelse(chromvar_windows<0.25,1,0) #score ranges from -1 to 1, using 0.25 for cutoff
  #require 2.5% cell coverage for windows, this is kinda on par with ATAC data, kinda an arbitrary cutoff. but mostly to decrease computational time on which windows we scan for motifs (5% is 10k window, 2% is 909k windows)
  counts<-counts[rowSums(!is.na(counts))>=(ncol(counts)/100)*percent_cell_coverage,] 
  counts[is.na(counts)]<-0
  #remove purely hypermethylated windows
  counts<-counts[rowSums(counts)>1,] 
  #dim(counts)
  return(counts)}

  chromvar_methylation<-function(obj,counts,prefix="allcells",threads){
    #prepare summarized experiment for chromvar
    peaks<-GenomicRanges::makeGRangesFromDataFrame(data.frame(
      seqnames=unlist(lapply(strsplit(row.names(counts),"_"),"[",1)),
      start=unlist(lapply(strsplit(row.names(counts),"_"),"[",2)),
      end=unlist(lapply(strsplit(row.names(counts),"_"),"[",3))))

    #prepare motifs
    opts <- list()
    opts[["species"]] <- "Homo sapiens"
    opts[["collection"]] <- "CORE"
    opts[["all_versions"]] <- FALSE
    motifs <- getMatrixSet(JASPAR2020,opts)

  #split peaks evenly into chunks so we can multicore the motif scanning
  motif_matches<-mclapply(split(peaks,  cut(seq_along(peaks), threads, labels = FALSE)),
                          function(x){
                          matchMotifs(motifs, x, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff=0.01)},
                          mc.cores=threads)
  motif_ix<-do.call("rbind",motif_matches)

  #create summarized experiment
  rse <- SummarizedExperiment::SummarizedExperiment(
                                  assays=list(counts=as(counts, "sparseMatrix")),
                                  rowRanges=peaks)
  colData(rse)<-as(obj@metadata[colnames(counts),],"DataFrame")
  rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
  dev <- computeDeviations(object = rse, annotations = motif_ix)
  saveRDS(dev,file=paste0(prefix,".chromvar.rds"))
  dev<-readRDS(file=paste0(prefix,".chromvar.rds"))

  #calculate variability
  variability <- computeVariability(dev)
  ggsave(plotVariability(variability, use_plotly = FALSE),file="chromvar_variability.pdf")

  #Differential motif analysis
  diff_acc <- differentialDeviations(dev, "cluster_id")
  diff_var <- differentialVariability(dev, "cluster_id")

  #differential tfbs by highest variability
  diff_tfbs<-row.names(variability[variability$variability>0.3,])
  devs<-deviations(dev)
  devs[is.na(devs)]<-0 #fill in NA for dev scores
  #dim_out<-irlba::irlba(devs[diff_tfbs,], 30)
  dim_out<-t(devs[diff_tfbs,])
  dim<-uwot::umap(dim_out)
  dim<-as.data.frame(dim)
  colnames(dim)<-c("chromvar_umap_x","chromvar_umap_y")
  row.names(dim)<-colnames(devscores)
  dim$cluster_id<-obj@metadata[row.names(dim),]$cluster_id
  plt<-ggplot(dim,aes(x=chromvar_umap_x,y=chromvar_umap_y,color=cluster_id))+geom_point()+theme_minimal()
  ggsave(plt,file=paste0(prefix,".chromvar_umap.pdf"))

  sample_cor <- getSampleCorrelation(dev,threshold=0.3)
  sample_cor[is.na(sample_cor)]<-0 #fill in na as 0 for sites with no overlap
  plt<-pheatmap(as.dist(sample_cor), 
          annotation_row = as.data.frame(colData(dev)[c("cluster_id","sample")]),
          clustering_distance_rows = as.dist(1-sample_cor), 
          clustering_distance_cols = as.dist(1-sample_cor))
  ggsave(plt,file=paste0(prefix,".chromvar_motifs.heatmap.pdf"))
  saveRDS(dev,file=paste0(prefix,".chromvar.rds"))
}




#compute synergy between motifs
#getAnnotationSynergy(rse, motif_ix[,c(83,24,20)])


```
# Published ATAC peaks
# running window clustering now
```R
source("~/projects/metact/src/amethyst_custom_functions.R")
setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
options(future.globals.maxSize = 20000 * 1024^2) #set parallelization max size to 20GB
#BiocManager::install("plyranges")

obj<-readRDS("dcis_hbca.amethyst.rds")

peaks_list<-list.files("/volumes/USR2/Ryan/projects/metact/ref/nakshatri",pattern="narrowPeak$",full.names=TRUE)
bed_list<-mclapply(peaks_list,rtracklayer::import,mc.cores=10)
names(bed_list)<-basename(unlist(lapply(strsplit(peaks_list,"_peaks"),"[",1)))

#if average peaks are less than 500, i would expand
bed_list<-mclapply(names(bed_list),function(x){
        mcols(bed_list[[x]])$celltype <- x
        return(bed_list[[x]])},
        mc.cores=10)

celltype_bed<-do.call("c",bed_list)
celltype_bed_reduced<-GenomicRanges::reduce(celltype_bed, with.revmap=TRUE)
mcols(celltype_bed_reduced) <- do.call(rbind,
                                mclapply(mcols(celltype_bed_reduced)$revmap, function(i) {
                                data.frame(celltype = paste(mcols(celltype_bed)$celltype[ i ], collapse = ","))},
                                mc.cores=50))

saveRDS(celltype_bed_reduced,file="/volumes/USR2/Ryan/projects/metact/ref/nakshatri/celltype_peaks.500bp.rds")

celltype_bed_reduced<-readRDS(file="/volumes/USR2/Ryan/projects/metact/ref/nakshatri/celltype_peaks.500bp.rds")
celltype_bed<-data.frame(seqnames=seqnames(celltype_bed_reduced),
  start=start(celltype_bed_reduced)-1,
  end=end(celltype_bed_reduced))

obj@genomeMatrices[["cg_100k_score"]]<-NULL
obj@genomeMatrices[["cg_celltypepeaks"]]<-NULL

obj<-cluster_by_windows(obj=obj,
                        window_name="cg_celltypepeaks",
                        prefix="dcis_hbca.celltypepeaks",
                        bed=celltype_bed,
                        threads=100,
                        neighbors=20,#20
                        k_phenograph=60,#70
                        distance=0.05 #0.05
                        )


saveRDS(obj,file="dcis_hbca.celltypepeaks.amethyst.rds")

obj<-dmr_and_1kb_window_gen(obj=obj,prefix="dcis_hbca.celltypepeaks") #running line by line
bigwig_output(obj=obj,tracks="cg_cluster_id_tracks")
saveRDS(obj,file="dcis_hbca.celltypepeaks.amethyst.rds")

obj<-readRDS(file="dcis_hbca.celltypepeaks.amethyst.rds")
celltypes<-unique(unlist(lapply(strsplit(celltype_bed_reduced$celltype,","),"[")))

avg_met_list<-mclapply(celltypes,function(i){
  print(i)
  celltype_bed<-celltype_bed_reduced[grepl(i,celltype_bed_reduced$celltype),]
  celltype_row_names<-paste(seqnames(celltype_bed),start(celltype_bed)-1,end(celltype_bed),sep="_")
  celltype_met<-obj@genomeMatrices[["cg_celltypepeaks"]][celltype_row_names %in% row.names(obj@genomeMatrices[["cg_celltypepeaks"]]),]
  avgmet<-setNames(nm=colnames(obj@genomeMatrices[["cg_celltypepeaks"]]),colMeans(celltype_met,na.rm=T))
  return(avgmet)
  },mc.cores=10)

avg_met<-do.call("cbind",avg_met_list)
colnames(avg_met)<-paste0(celltypes,"_avgpeakmet")
avg_met<-as.data.frame(t(scale(t(avg_met))))#scaled per cell to account for average met
obj@metadata[paste0(celltypes,"_avgpeakmet")]<-avg_met[row.names(obj@metadata),]

plt_list<-lapply(paste0(celltypes,"_avgpeakmet"),function(i){
  plt<-ggplot() + 
  geom_point(aes(
    x=as.numeric(obj@metadata$umap_x),
    y=as.numeric(obj@metadata$umap_y),
    color=as.numeric(obj@metadata[[i]]))) + 
  theme_void() + 
  scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + 
  ggtitle(i)
  return(plt)})
plt<-patchwork::wrap_plots(plt_list)
ggsave(plt,file="avg_methylation.pdf",width=30,height=30)



#Methylation average over cell type peaks for plotting?

```
# TO DO: CELL TYPING BY PEAK CALLING (JUST MAKE A PEAK OVERLAP METRIC OR SOMETHING)  (AVG methylation over peakssets?)
# CHROMVAR RERUN ON ALL CELLS

# TO DO: ADD MEME SCANNING TO DMRS USING JASPAR 2024 MOTIFS



```bash
singularity shell \
--bind ~/projects/ \
--bind /volumes/seq/projects/metACT \
~/projects/metact/src/amethyst.sif
source activate base
conda activate MethylTree

cd /volumes/USR2/Ryan/projects/metact/amethyst_processing
python ~/projects/metact/src/methyltree_subclones.py -i DCIS-41T
python ~/projects/metact/src/methyltree_subclones.py -i DCIS-66T
python ~/projects/metact/src/methyltree_subclones.py -i IDC_DCIS-79T

```

```python

```


#################### #################### Additional Clustering Promoters, Genes, Promoters without CGI #################### #################### 

#protein coding gene promoters
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
hbca@ref<-gtf

cgi<-rtracklayer::import(BEDFile("/volumes/USR2/Ryan/projects/metact/ref/cgi_hg38.bed.gz"))

#Gene names for clustering
for (i in c("gene_name", "exon_number")) {gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
hbca@ref<-gtf
protein_coding <- unique(obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> dplyr::pull(gene_name))

#Promoter bed with cgi overlap info
protein_coding_gtf<-gtf[gtf$type=="gene" & gtf$gene_type=="protein_coding" & gtf$seqid!="chrM",]
protein_coding_gtf<-protein_coding_gtf[!duplicated(protein_coding_gtf$gene_name),]
protein_coding_gtf<-makeGRangesFromDataFrame(protein_coding_gtf,keep.extra.columns=TRUE)
protein_coding_gtf<-extend(protein_coding_gtf,downstream=20000) #extend into gene body
protein_coding_gtf$cgi_overlap<-!is.na(findOverlaps(protein_coding_gtf,cgi,ignore.strand=TRUE,type="any",select="first")) #filter out promoters with CGI
promoter_bed<-data.frame(seqnames=seqnames(protein_coding_gtf),start=start(protein_coding_gtf)-1,end=end(protein_coding_gtf))

hbca<-cluster_by_windows(obj=hbca,
                        window_name="cg_promoter_score",
                        bed=promoter_bed,
                        threads=200,
                        neighbors=20,
                        distance=0.05)


hbca<-cluster_by_windows(obj=hbca,
                        window_name="cg_genebody",
                        genes=protein_coding,
                        promoter=FALSE,
                        threads=200,
                        neighbors=20,
                        distance=0.05)
                        
