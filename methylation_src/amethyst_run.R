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

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
setwd(in_dir2)
obj <- createObject()

#read in metadata
samp_dir=list.files(path=paste0(in_dir,"/report/"))
samp_dir2=list.files(path=paste0(in_dir2,"/report/"))
met1<-lapply(samp_dir,function(i){
  samp_met<-read.csv(paste0(in_dir,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat1<-do.call("rbind",met1)
metadat1$run<-1
metadat1<-metadat1[metadat1$CG_Cov>10000,]
metadat1<-metadat1[metadat1$sampleName %in% c("hbca-83l","hbca-16r","DCIS-41T","DCIS-66T"),]

met2<-lapply(samp_dir2,function(i){
  samp_met<-read.csv(paste0(in_dir2,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat2<-do.call("rbind",met2)
metadat2$run<-2
metadat2<-metadat2[metadat2$CG_Cov>10000,]
metadat2<-metadat2[metadat2$sampleName %in% c("hbca-83l" ,"hbca-16r" ,"DCIS-41T", "DCIS-66T"),]

metadat<-rbind(metadat1,metadat2)
#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
metadat$mcg_pct<-metadat$CG_mC_Pct
metadat$cov<-metadat$CG_Cov
obj@metadata<-metadat

head(obj@metadata)
obj@metadata<-obj@metadata[obj@metadata$CG_Cov>10000,]
plt<-ggplot(obj@metadata, aes(x=sampleName, y = CG_Cov)) +geom_violin() + geom_jitter()
ggsave(file="cov_plot.pdf",plt)

#filter cells
#obj@metadata <- obj@metadata |> dplyr::filter(cov > 100000 & cov < 40000000)
h5paths<-paste0(in_dir,"/cg_sort_cov/h5_files/",metadat1$sampleName,".",metadat1$tgmt_well,".h5")
h5paths2<-paste0(in_dir2,"/h5_files/",metadat2$sampleName,".",metadat2$tgmt_well,".h5")
obj@h5paths <- data.frame(row.names = c(rownames(metadat)), paths = c(h5paths,h5paths2))
head(obj@h5paths)

# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 50) 

cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=100,neighbors.=50){
  print(paste("Making window summaries for ",window_name))
  obj@genomeMatrices[[window_name]] <- makeWindows(obj,
                                                      stepsize = stepsize., 
                                                      type = "CG", 
                                                      metric = metric., 
                                                      bed = bed.,
                                                      threads = threads., 
                                                      index = "chr_cg", 
                                                      nmin = 2, 
                                                      species = "human") 
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
  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sampleName, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj_out<-cluster_by_windows(obj,window_name="cg_50k_score",stepsize.=50000,threads.=150)


#run on windows (bed format)
#download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/pbmc_vignette/pbmc_highconfidence_dmrs.bed", "./Downloads/pbmc_dmr.bed") 
#obj@genomeMatrices[["pbmc_dmrs"]] <- makeWindows(obj, bed = "./Downloads/pbmc_dmr.bed", type = "CG", metric = "percent", threads = 1, index = "chr_cg", nmin = 2, species = "human")

#obj@genomeMatrices[["pbmc_dmrs"]] <- obj@genomeMatrices[["pbmc_dmrs"]][rowSums(!is.na(obj@genomeMatrices[["pbmc_dmrs"]])) >= 10, ] 
#dimEstimate(obj, genomeMatrices = c("pbmc_dmrs"), dims = c(10), threshold = 0.90)
#set.seed(111)
#obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = c("pbmc_dmrs"), dims = c(9), replaceNA = c(0))
#obj@reductions[["irlba_regressed"]] <- regressCovBias(obj, reduction = "irlba")
#obj <- runCluster(obj, k_phenograph = 10, reduction = "irlba_regressed") # consider increasing k_phenograph to 50 for larger datasets
#obj <- runUmap(obj, neighbors = 5, dist = 0.05, method = "euclidean", reduction = "irlba_regressed") 
#plt1<-dimFeature(obj, colorBy = cluster_id, reduction = "umap")
#plt2<-dimFeature(obj, colorBy = cluster_id) + facet_wrap(vars(batch)) # Batch is simulated to illustrate function utility. Any column in the metadata will work.
#plt<-plot_grid(plt1, plt2)
#ggsave(plt,file="./Downloads/umap_dmr.pdf")

#plot by feature
plt<-sampleComp(obj, groupBy = "sampleName", colorBy = "cluster_id") 
ggsave(plt,file="sample_comp.pdf")


#annotate data
#obj@ref <- makeRef(ref="HG38",gtf="/volumes/USR2/Ryan/ref/gencode.v43.annotation.gtf") #ref stored within amethyst.sif container
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
for (i in c("gene_name", "exon_number")) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, 
        i))
}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
obj@ref<-gtf

cluster1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 100,
                                         step = 1000,
                                         smooth = 3,
                                         species = "human",
                                         index = "chr_cg",
                                         groupBy = "cluster_id",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)

obj@genomeMatrices[["cg_cluster_tracks"]] <- cluster1kbwindows[["pct_matrix"]]


library(plyr)
#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("COBL","GABRP","ELF5","CCL28","KRT15","KIT") #"AC011247.1",
hbca_snmarkers[["basal"]]=c("CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2") #"AC044810.2",
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1") #"AL357507.1",
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1") #"AC012409.2"
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-llply(hbca_snmarkers, unlist)

protein_coding <- unique(obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> dplyr::pull(gene_name))

#plot histogram of markers 
plt_list<-lapply(names(features), function(i) {
  feat=unlist(features[[i]])
  feat<-feat[feat%in% protein_coding]
  plt<-histograM(obj, genes = feat, trackOverhang=20000,
        matrix = "cg_cluster_tracks", 
        width = 1000)
  ggsave(plt,file=paste0("marker_histograms.",i,".pdf"),width=30)
})

#plot heatmap of markers 
plt_list<-lapply(names(features), function(i) {
  feat=unlist(features[[i]])
  feat<-feat[feat%in% protein_coding]
  plt<-heatMap(obj, genes = feat,
        matrix = "cg_cluster_tracks", trackOverhang=20000,
        width = 1000)
  ggsave(plt,file=paste0("marker_heatmaps.",i,".pdf"),width=5,height=10)
})

obj@genomeMatrices[["cg_promoter"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = TRUE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 100, 
                                                     index = "chr_cg", 
                                                     nmin = 2, 
                                                     species = "human") 

obj@genomeMatrices[["cg_genebody"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 100, 
                                                     index = "chr_cg", 
                                                     nmin = 2, 
                                                     species = "human") 


plt_list<-lapply(names(features), function(i) {
  plt<-dotM(obj, genes = features[[i]], groupBy = "cluster_id", matrix = "cg_promoter") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 8)) +coord_flip() + ggtitle(i)
  return(plt)
})

plt_out<-plot_grid(plotlist=plt_list,ncol=1)
ggsave(plt_out,file="marker_promoter_dotplot.pdf",width=10,height=length(features)*5,limitsize=F)


plt_list<-lapply(names(features), function(i) {
  plt<-dotM(obj, genes = features[[i]], groupBy = "cluster_id", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 8)) +coord_flip() + ggtitle(i)
  return(plt)
})

plt_out<-plot_grid(plotlist=plt_list,ncol=1)
ggsave(plt_out,file="marker_genebody_dotplot.pdf",width=10,height=length(features)*5,limitsize=F)


#find cluster markers on promoters
cluster_promoter_markers <- findClusterMarkers(obj, 
                                               matrix = "cg_promoter", 
                                               genes = protein_coding, 
                                               threads = 100)
cluster_promoter_markers <- cluster_promoter_markers |> dplyr::filter(p.val < 0.1)
head(cluster_promoter_markers)
posteriori_markers<-cluster_promoter_markers |> dplyr::filter(p.val < 0.05)  |> dplyr::filter(is.finite(logFC)) |> dplyr::filter(!duplicated(gene)) |> dplyr::filter(direction=="hypomethylated") |> group_by(cluster_id) |> arrange(p.val) %>% slice(1:10)
plt_list<-lapply(unique(posteriori_markers$cluster_id), function(i) {
  plt<-dotM(obj, genes = posteriori_markers[posteriori_markers$cluster_id == i,]$gene, groupBy = "cluster_id", matrix = "cg_promoter") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "black")) + scale_size(range = c(1, 15)) + guides(fill="none")
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,nrow=1)
ggsave(plt_out,file="denovomarkers_promoters_dotplot.pdf",height=10,width=length(features)*5,limitsize=F)


#find cluster markers on genebodies
cluster_genebody_markers <- findClusterMarkers(obj, 
                                               matrix = "cg_genebody", 
                                               genes = protein_coding, 
                                               threads = 100)
cluster_genebody_markers <- cluster_genebody_markers |> dplyr::filter(p.val < 0.1)
head(cluster_genebody_markers)
posteriori_markers<-cluster_genebody_markers |> dplyr::filter(p.val < 0.05)  |> dplyr::filter(is.finite(logFC)) |> dplyr::filter(!duplicated(gene)) |> dplyr::filter(direction=="hypomethylated") |> group_by(cluster_id) |> arrange(p.val) %>% slice(1:10)
plt_list<-lapply(unique(posteriori_markers$cluster_id), function(i) {
  plt<-dotM(obj, genes = posteriori_markers[posteriori_markers$cluster_id == i,]$gene, groupBy = "cluster_id", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "black")) + scale_size(range = c(1, 15)) + guides(fill="none")
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,nrow=1)
ggsave(plt_out,file="denovomarkers_genebody_dotplot.pdf",height=10,width=length(features)*5,limitsize=F)



#DMR calculation
pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
dmrs <- testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 2, nminGroup = 2) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs2 <- filterDMR(dmrs, method = "bonferroni", filter = FALSE) #add additional columns direction column

head(dmrs2)
collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 2000, reduce = T, annotate = T) 

head(collapsed_dmrs)
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
       aes(y = test, x = n, fill = test)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = makePalette(option = 7, n = 8) ) + theme_classic()
ggsave(plt,file="met_per_dmr.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(test, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(test, direction) |> slice_min(n = 5, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)

plt<-heatMap(obj, matrix = "cg_cluster_tracks", regions = top_dmrs$location, nrow = 10, legend = FALSE, width = 1000, arrowOverhang = 2000)
ggsave(plt,file="dmrs_heatmap.pdf",width=20,height=30)


#get overlapping genes with DMRs
dmr_genes<-collapsed_dmrs |> dplyr::filter(dmr_padj < 0.05) |> dplyr::filter(is.finite(dmr_logFC)) |> dplyr::filter(gene_names!="NA") |> dplyr::filter(direction=="hypo") |> group_by(test) |> arrange(dmr_padj) %>% slice(1:10)
dmr_genes$gene<-unlist(lapply(dmr_genes$gene_names, function (x) {
  if(grepl(",",x)){
    out_gene<-strsplit(x,",")[[1]][1]
  }else {out_gene<-x}
  return(out_gene)}))

#plot heatmap of dmr marker genes
plt_list<-lapply(unique(dmr_genes$test), function(i) {
  feat=dmr_genes[dmr_genes$test==i,]$gene
  feat<-feat[feat%in% protein_coding]
  plt<-histograM(obj, genes = feat, trackOverhang=20000,
        matrix = "cg_cluster_tracks", 
        width = 1000,legend=FALSE) 
  ggsave(plt,file=paste0("dmr_histograms_wiggle.",i,".pdf"),width=30)
})



#Recluster on DMRs

#Recluster on DMRs given by Lauren

#Recluster on ATAC peaks



#DFFF00