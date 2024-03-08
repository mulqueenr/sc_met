#load modules
#module load nextflow/23.04.3
#module load singularity

#set up environment variables 
#export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"
#export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
#export srcDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src"
#export refDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref"
#export sif="${srcDir}/copykit.sif"

#cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cg_dataframes
#export sif="${srcDir}/multiome_bc.sif"

#singularity shell \
#--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2 \
#$sif


#zcat mc_posteriorest.100kb.merged.tsv.gz | head -n 1000 | gzip > TEST.mc_posteriorest.100kb.merged.tsv.gz
#zcat mc_posteriorest.genebody.merged.tsv.gz | head -n 1000 | gzip > TEST.mc_posteriorest.genebody.merged.tsv.gz

library(Seurat)
library(data.table)

meta_files<-list.files(recursive=TRUE,path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/report/",pattern="*.passingCellsMapMethylStats.csv",full.names=TRUE)


met<-do.call("rbind",lapply(meta_files,function(x) read.csv(sep=",",file=x,header=T)))
row.names(met)<-paste0(met$sampleName,met$BC)



read_in_data<-function(x){
	dat<-as.data.frame(fread(x,showProgress=TRUE,nThread=5))
	row.names(dat)<-dat[,1]
	dat<-dat[2:ncol(dat)]
	return(dat)
}

hundokb_counts<-read_in_data("total_counts.100kb.merged.tsv.gz")
hundokb_dat<-read_in_data("mc_posteriorest.100kb.merged.tsv.gz")
gene_counts<-read_in_data("total_counts.genebody.merged.tsv.gz")
gene_dat<-read_in_data("mc_posteriorest.genebody.merged.tsv.gz")

seurat_object <- CreateSeuratObject(counts = t(hundokb_counts),meta.data=met,assay="met_100kb")
seurat_object[["met_100kb"]]
seurat_object[["met_gene"]] <- CreateAssayObject(counts=t(gene_counts),data=t(gene_dat))

saveRDS(seurat_object,file="met.SeuratObject.rds")

