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
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 200) 


cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=200,neighbors.=50,prefix){
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
  ggsave(plt,file=paste0(prefix,window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj<-cluster_by_windows(obj,window_name="cg_100k_score",stepsize.=100000,threads.=200,prefix="cell_line")

