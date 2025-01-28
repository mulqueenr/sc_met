
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

# Methyltree run for clonal analysis

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

#cell type assignment from bigwigs
celltypes<-setNames(nm=as.character(1:15),c("LumHRCancer_4","Fibro","LumHRNormal","Tcell_3","Myeloid","Basal_1","Endo","LumSec_2","Tcell_1","LumHRCancer_2","LumHRCancer_1","LumSec_1","Basal_2","LumHRCancer_3","Tcell_2"))
obj@metadata$celltype<-celltypes[obj@metadata$cluster_id]
saveRDS(obj,file="dcis_hbca.celltypepeaks.amethyst.rds")

cluster1kbwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 300,
                                        step = 500,
                                        smooth = 3,
                                        index = "chr_cg",
                                        groupBy = "celltype", 
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)

obj@genomeMatrices[[paste0("cg_",groupBy,"_tracks")]] <- cluster1kbwindows[["pct_matrix"]]
saveRDS(obj,file="dcis_hbca.celltypepeaks.amethyst.rds")

plt<-histograM(obj,gene="COL1A1",matrix=paste0("cg_",groupBy,"_tracks"))
ggsave(plt,file="test_met_cov.pdf")
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


# Chromvar motif analysis

```R
source("/volumes/USR2/Ryan/projects/metact/src/amethyst_custom_functions.R") #this loads all necessary R packages and custom functions

options(future.globals.maxSize = 20000 * 1024^2) #set parallelization max size to 20GB

setwd("/volumes/USR2/Ryan/projects/metact/amethyst_processing") #set wd
obj<-readRDS(file="dcis_hbca.celltypepeaks.amethyst.rds")
threads=50
prefix="allcells"

#run per cell
counts<-chromvar_met_per_cell(obj,stepsize=500,threads=50,percent_cell_coverage=2.5)
chromvar_methylation(obj=obj,counts=counts,prefix="allcells_cells",threads=50)

#run on the cluster level for less coverage noise
counts<-chromvar_met_per_cluster(obj,stepsize=500,threads=50,mincov=30,percent_cell_coverage=50,groupBy="cluster_id")
chromvar_methylation(obj=obj,counts=counts,prefix="allcells_cluster",threads=50)

dcis<-subsetObject(obj, cells=row.names(obj@metadata[startsWith(toupper(obj@metadata$sample),prefix="DCIS-41T"),]))
counts<-chromvar_met_per_cluster(dcis,stepsize=500,threads=100,mincov=2,percent_cell_coverage=50,groupBy="subclones")
#run chromvar on just dcis-41t
chromvar_methylation(obj=obj,counts=counts,prefix="allcells_cluster",threads=50)

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
                        
