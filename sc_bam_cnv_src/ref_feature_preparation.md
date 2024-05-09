make transcript bed
```bash
gtf="/volumes/seq/projects/metACT/ref/grch38/filteredGTF/GRCh38_transcriptsOnly.gtf"

awk 'OFS="\t" {split($10,a,"\""); split($14,b,"\""); print "chr"$1,$4,$5,a[2],b[2]}' $gtf > GRCh38_transcripts.bed


select_longest() { 
	awk -v gene="$1" '{if($5==gene){print $0,$3-$2}}' GRCh38_transcripts.bed \
	| sort -k6,6n - \
	| tail -n 1
}

export -f select_longest

parallel -j 20 select_longest ::: $(awk '{print $5}' GRCh38_transcripts.bed | uniq) > GRCh38_transcripts.longest.bed

```

make 100kb bed
```bash
ref="/volumes/seq/projects/metACT/ref/grch38/genome.txt"
grep -v "^K" ${ref} \
| grep -v "^G" \
| awk 'OFS="\t" {print "chr"$1,$2}' > genome.filt.txt

bedtools makewindows -w 100000 -g genome.filt.txt | awk 'OFS="\t" {print $1,$2,$3,$1"_"$2"_"$3}' > genome_windows.100kb.bed
```

met atlas
```bash
bsub -Is -W 4:00 -q transfer -n 10 -M 10 -R rusage[mem=10] /bin/bash
dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects"
module load bedtools
module load parallel

cd ${dir}/metact/ref
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM5_ESM.zip
unzip 41586_2022_5580_MOESM5_ESM.zip -d met_atlas

merge_bed () {
	in=$1
	celltype=$(printf '%s\n' "$in" | cut -d. -f1)
	awk -v outcell="$celltype" 'OFS ="\t" {print $1,$2,$3,$6,$7,outcell}' $in | tail -n +2
}
# add header and sort bed file, merge overlapping sites
export -f merge_bed

#header=$(echo chr$'\t'start$'\t'end$'\t'feat_description$'\t'gene$'\t'celltype)
(parallel merge_bed ::: *bed | sort -k 1,1 -k2,2n -k3,3n - | grep "^chr[1-9|X]" | bedtools merge -c 4,5,6 -o first,first,collapse) > ../met_atlas.hg19.bed

#perform liftOver
#set up


cd ${dir}/metact/src
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+wrx liftOver
cd ${dir}/metact/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

#run
cd ${dir}/metact/ref
${dir}/metact/src/liftOver -tab -bedPlus=3 \
met_atlas.hg19.bed \
hg19ToHg38.over.chain.gz \
met_atlas.hg38.bed \
met_atlas.hg19.unmapped.bed

cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref
awk 'OFS="\t" {print $1,$2,$3,$1"_"$2"_"$3}' met_atlas.hg38.bed > met_atlas.hg38.feat.bed


```

atac-seq peaks from snubar HBCA

NOTE: I ended up gettting the data from our GEO server, since the release didn't have cell type information.
Original file was located here: /volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn/ready.coassay_hto.sr4.rds

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5484482
```bash
bsub -Is -W 4:00 -q transfer -n 10 -M 10 -R rusage[mem=10] /bin/bash
dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects"

mkdir ${dir}/metact/ref/snubar_atac
cd ${dir}/metact/ref/snubar_atac

sftp mulqueen@qcprpgeo.mdanderson.edu

get /volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn/ready.coassay_hto.sr4.rds
get /volumes/USR1/yyan/project/coda/rds/hbca_10x_hg19/tosignac/ready.signac_binary.rds


#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5484nnn/GSM5484482/suppl/GSM5484482%5FSNuBarARC%2DHBCA.scATAC.cell%5Fmetainfo.csv.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5484nnn/GSM5484482/suppl/GSM5484482%5FSNuBarARC%2DHBCA.scATAC.filtered%5Fpeak%5Fbc%5Fmatrix.barcodes.tsv.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5484nnn/GSM5484482/suppl/GSM5484482%5FSNuBarARC%2DHBCA.scATAC.filtered%5Fpeak%5Fbc%5Fmatrix.peaks.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5484nnn/GSM5484482/suppl/GSM5484482%5FSNuBarARC%2DHBCA.scATAC.fragments.tsv.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5484nnn/GSM5484482/suppl/GSM5484482%5FSNuBarARC%2DHBCA.scATAC.fragments.tsv.gz.tbi.gz

#gzip -d GSM5484482_SNuBarARC-HBCA.scATAC.fragments.tsv.gz.tbi.gz
```
-->
```bash
bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash
dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects" #
cd ${dir}/metact/ref/snubar_atac


#load modules
module load singularity

singularity shell \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact \
~/singularity/multiome_bc.sif

```

Seadragon Version
```R
library(Seurat)
library(Signac)

proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref/snubar_atac" #
setwd(proj_dir)

snubar<-readRDS("ready.coassay_hto.sr4.rds") #snubar<-readRDS("/volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn/ready.coassay_hto.sr4.rds")
hbca_atac<-readRDS("ready.signac_binary.rds") #geo version 
library(Signac); library(Seurat)
obja <- read_rds('/volumes/USR1/yyan/project/coda/rds/hbca_10x_hg19/tosignac/ready.signac_binary.rds')
UMAPPlot(obja, group.by='celltypes', label=T)
UMAPPlot(obja, group.by='cellstates', label=T)

Idents(hbca_atac)<-hbca_atac$celltypes

hbca_atac_markers<-FindAllMarkers(hbca_atac,test='LR',latent.vars="nCount_peaks")
```

GEO Version

```bash
singularity shell \
--bind /volumes/seq/projects/metACT/ \
--bind /volumes/USR1/yyan/project/coda/rds/hbca_10x_hg19/tosignac \
--bind /volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn \
~multiome_bc.sif
```
```R
library(Seurat)
library(Signac)

proj_dir="/volumes/seq/projects/metACT/ref/snubar_atac"
setwd(proj_dir)


snubar<-readRDS("/volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn/ready.coassay_hto.sr4.rds")
hbca_atac<- readRDS('/volumes/USR1/yyan/project/coda/rds/hbca_10x_hg19/tosignac/ready.signac_binary.rds')

Idents(hbca_atac)<-hbca_atac$celltypes; DefaultAssay(hbca_atac)<-"peaks"
Idents(snubar)<-snubar$celltype; DefaultAssay(snubar)<-"peaks"

snubar_atac_markers<-FindAllMarkers(snubar,test='LR',latent.vars="nCount_peaks")

hbca_atac_markers<-FindAllMarkers(hbca_atac,test='LR',latent.vars="nCount_peaks")

```

GEO Version with hg38 reference
Using ArchRtoSignac to convert Arrow files
https://github.com/swaruplabUCI/ArchRtoSignac

```bash
conda create -n scATAC -c conda-forge r-base r-essentials
conda activate scATAC
```
```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# install ArchRtoSignac
devtools::install_github("swaruplabUCI/ArchRtoSignac")
# load ArchRtoSignac
library(ArchRtoSignac)

packages <- c("ArchR","Seurat", "Signac","stringr") # required packages
loadinglibrary(packages)
```
singularity shell \
--bind /volumes/seq/projects/metACT/ \
--bind /volumes/USR1/yyan/project/coda/rds/hbca_10x_hg38/archr_project \
--bind /volumes/USR1/yyan/project/snubar_atac/rds_analysis_ready/coassay_4pat/coassay/postQC/wnn \
~multiome_bc.sif
```


/volumes/USR1/yyan/project/coda/rds/hbca_10x_hg38/archr_project/ArrowFiles



```


# ATAC Peaks from Multiome Paper

On exacloud
```R

library(Signac)
library(Seurat)
library(reshape2)
library(parallel)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.SeuratObject.rds"

dat=readRDS(opt$object_input)
dat<-subset(dat,Diagnosis %in% c("NAT","DCIS"))
#remove unnecessary assays, so I can better parallelize
dat[["RNA"]]<-NULL;dat[["SoupXRNA"]]<-NULL;dat[["chromvar"]]<-NULL;dat[["GeneActivity"]]<-NULL;dat[["SCT"]]<-NULL;
DefaultAssay(dat)<-"peaks"
Idents(dat)<-dat$HBCA_predicted.id

parallelized_markers<-function(obj,x){
  markers<-FindMarkers(
    obj, assay="peaks",
    slot = "data",
    ident.1 = x,
    test.use = "LR",
    only.pos=T
  )
  markers$ident<-x
return(markers)
}

out<-mclapply(unique(as.character(Idents(dat))),function(x) parallelized_markers(dat,x),mc.cores=5)


saveRDS(out,"multiome_peak_markers.rds")



```
Transferred to geo
```R
dat<-readRDS("multiome_peak_markers.rds")

idents<-c("lumepi","basalepi","myeloid", "endovascular", "adipo", "fibro", "endolymphatic", "pericyte", "Tcell")           # this was taken from unique(Idents(dat)) call order 
dat<-lapply(1:length(dat), function(x) {dat[[x]]$ident<-idents[x]; return(dat[[x]])})
dat<-do.call("rbind",dat)

dat<-dat[dat$p_val<0.05,]

#Taking even nominally significant output
#table(dat$ident)

#        Tcell         adipo      basalepi endolymphatic  endovascular 
#         4565         10243         35383          6732         18665 
#        fibro        lumepi       myeloid      pericyte 
#        32948         78947          9071          9454 

out<-data.frame(chr=unlist(lapply(strsplit(row.names(dat),"-"),"[",1)),
	start=unlist(lapply(strsplit(row.names(dat),"-"),"[",2)),
	end=unlist(lapply(strsplit(row.names(dat),"-"),"[",3)),
	feat_name=paste(dat$ident,1:nrow(dat),sep="_"))

write.table(out,file="multiome_bc.marker_atac_peaks.bed",sep="\t",col.names=F,row.names=F,quote=F)

#seadragon file located here:
#/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref/multiome_bc.marker_atac_peaks.bed
```