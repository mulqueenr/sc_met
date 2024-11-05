Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On seadragon.

```bash
#ssh RMulqueen@seadragon
proj_dir="/volumes/seq/projects/metACT"
mkdir ${proj_dir}

#download https://github.com/ScaleBio/ScaleMethyl/blob/main/docs/scaleBioTools.md
#download binaries from /volumes/seq/projects/metACT/tools/ScaleMethyl/envs/download-scale-tools.sh
/volumes/seq/projects/metACT/tools/ScaleMethyl/envs/download-scale-tools.sh
#conda update -n base -c defaults conda
#conda create -n scalemet
#conda activate scalemet
#conda install mamba
#mamba env create -f /volumes/seq/projects/metACT/tools/ScaleMethyl/envs/scaleMethyl.conda.yml
#mamba env create -f /volumes/seq/projects/metACT/tools/ScaleMethyl/envs/scaleMethylPyQc.yml


#install bcl-convert (need root access)
which bcl-convert
```

Running on geo

Generate job file scalemethyl_pipeline.sh
```bash

proj_dir="/volumes/seq/projects/metACT"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/240115_RMMM_scalebiotest2"
fastqDir="${proj_dir}/240111_VH00219_563_AAFFFCCM5"
samples="${proj_dir}/samples.csv"
genome="${proj_dir}/ref/genome.json"
fq_out="${runDir}/231204_scalebio_test/fastq"

#ln -s /volumes/seq/flowcells/MDA/nextseq2000/2024/240111_VH00219_563_AAFFFCCM5 ${proj_dir}/.

cd $runDir

#ran first and found all the reads assigned to undetermined, then looked up resulting indexes
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,1
TrimUMI,0
MinimumTrimmedReadLength,16
MaskShortReads,16
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,2
OverrideCycles,Y150;I10;I10;U8Y142
[BCLConvert_Data]
Sample_ID,index,index2
ScaleMethyl,GGAGGCCTCC,ATGGAGCTAC
ScaleMethyl,GGAGGCCTCC,TCTCATTGCC
ScaleMethyl,GGAGGCCTCC,GGCATAACCG
ScaleMethyl,GGAGGCCTCC,TAATTCCGGT
ScaleMethyl,GGAGGCCTCC,TGCGCCGAAG
ScaleMethyl,GGAGGCCTCC,CCTGCGTATT
ScaleMethyl,GGAGGCCTCC,CTAAGAGTTA
ScaleMethyl,GGAGGCCTCC,AAGCCTACGA""" > ${runDir}/samplesheet.csv

#index (this one is symmetric so hard to tell if it is actually revcomp)
#GGAGGCCTCC

#index2 (our index2 are revcomp from the listed ones on their github)
#https://github.com/ScaleBio/ScaleMethyl/blob/main/references/i5.txt
#GTAGCTCCAT
#GGCAATGAGA
#CGGTTATGCC
#ACCGGAATTA
#CTTCGGCGCA
#AATACGCAGG
#TAACTCTTAG
#TCGTAGGCTT
#run bcl2fastq myself 
bcl-convert \
--bcl-input-directory ${fastqDir} \
--bcl-num-conversion-threads 10 \
--bcl-num-compression-threads 10 \
--bcl-num-decompression-threads 10 \
--output-directory ${runDir} \
--sample-sheet ${runDir}/samplesheet.csv \
--force

rm -rf Undetermined*fastq.gz

echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-2H12,ScaleMethyl""" > ${runDir}/samples.csv

cd $runDir

#run til failure (combined report)
nextflow run ${scalebio_nf} \
--fastqDir ${runDir} \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--genome ${genome} \
--fastqOut true \
--trimOut true \
--bamOut true \
--bamDedupOut true \
--covOut true \
--max_memory 200.GB \
--max_cpus 40 \
-profile docker,singularity \
-w ~/tmp \
-resume

#get count of dedup reads per bam
proj_dir="/volumes/seq/projects/metACT"
runDir="${proj_dir}/240115_RMMM_scalebiotest2"
singularity shell /volumes/USR2/Ryan/tmp/singularity/public.ecr.aws-o5l3p3e4-scale_methyl-v1.5.img
export PATH="$PATH:/volumes/seq/projects/metACT/tools/ScaleMethyl/bin"

python library_complexity_and_stats.py --sampleName MCF10A --outDir MCF10A --tgmt_barcodes references/tgmt.txt --topCellPercentage 99 --threshold 0    --minCellRatio 20 --minReads 1000 --minUniqTotal 1 --maxUniqTotal 100 --i5_barcodes references/i5.txt --i7_barcodes references/i7.txt

```
```bash
count_scbam() {
	samtools view $1 | awk '{split($1,a,":");print a[1]}' | sort -T . | uniq -c | awk -v var="$1" 'OFS="\t" {print $1,$2,var}'
}
export -f count_scbam
bam_in=$(find -L -iname "*.bam")
parallel --jobs 50 count_scbam ::: $bam_in | sort -k1,1n > cell_readcounts.tsv
parallel --jobs 30 count_scbam ::: $bam_in | sort -k1,1n | awk '{if($1>=10000) print $0}' > cell_id.pf.tsv
```
```bash
#split out single cell bams
cd $runDir

mkdir -p $runDir/final_sc_bams
split_finalscbams() {
	proj_dir="/volumes/seq/projects/metACT"
	runDir="${proj_dir}/240115_RMMM_scalebiotest2"
	outname=$(basename $3)
	((samtools view -H $3) && (samtools view $3 | awk -v cell_id="$2" 'OFS="\t" {split($1,a,":"); if(a[1]==cell_id) print $0}')) | samtools view -bS | samtools sort -T . -o $runDir/final_sc_bams/${2}.${outname}

}

export -f split_finalscbams
parallel --colsep "\t" --jobs 30 split_finalscbams :::: cell_id.pf.tsv
```
```bash
#split out single cell bams (prededup)
cd $runDir

mkdir -p $runDir/prededup_sc_bams

split_prededupscbams() {
	proj_dir="/volumes/seq/projects/metACT"
	runDir="${proj_dir}/240115_RMMM_scalebiotest2"
	prededup_bam=$(echo $3 | awk '{gsub("bamDeDup","bsbolt",$0);gsub(".dedup.nsrt.bam",".bam",$0); print }')
	outname=$(basename $prededup_bam)
	((samtools view -H $prededup_bam) && (samtools view $prededup_bam | awk -v cell_id="$2" 'OFS="\t" {split($1,a,":"); if(a[1]==cell_id) print $0}')) | samtools view -bS | samtools sort -T . -o $runDir/prededup_sc_bams/${2}.${outname}
}

export -f split_prededupscbams
parallel --colsep "\t" --jobs 100 split_prededupscbams :::: cell_id.pf.tsv
```

```bash
#perform picard tools complexity test
dir=$runDir/prededup_sc_bams
cd $dir


project_count() {
outname=${1::-4}
samtools sort -T . -n -o - $1 | samtools fixmate -m - - | samtools sort -T . -o - - | samtools markdup -s - ${1::-4}.rmdup.bam 2> ${1::-4}.rmdup.stats.txt

java -jar /volumes/seq/code/3rd_party/picard/picard-2.20.4/picard.jar EstimateLibraryComplexity MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 I=${1::-4}.rmdup.bam O=${1::-4}.complex_metrics.txt

awk -v outname=$outname 'NR==8 {print outname,$3,$9,$10}' ${1::-4}.complex_metrics.txt
}
export -f project_count

bam_in=`ls *bam`
parallel --jobs 20 project_count ::: $bam_in > cell_library_sizes.txt
```

```R
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("/volumes/seq/projects/metACT/240115_RMMM_scalebiotest2/prededup_sc_bams")
dat<-read.table("cell_library_sizes.txt")
colnames(dat)<-c("cellid","uniq_reads","dup_perc","projected_uniq_reads")
dat$sample<-unlist(lapply(strsplit(dat$cellid,"[.]"),"[",2))

plt1<-ggplot(dat,aes(x=sample,y=uniq_reads))+geom_jitter(aes(color=sample))+geom_boxplot(aes(color=sample),fill=NA,outlier.shape=NA)+theme_minimal()
plt2<-ggplot(dat,aes(x=sample,y=dup_perc))+geom_jitter(aes(color=sample))+geom_boxplot(aes(color=sample),fill=NA,outlier.shape=NA)+theme_minimal()
plt3<-ggplot(dat,aes(x=sample,y=projected_uniq_reads))+geom_jitter(aes(color=sample))+geom_boxplot(aes(color=sample),fill=NA,outlier.shape=NA)+theme_minimal()

plt<-plt1/plt2/plt3
ggsave(plt,file="scalebio_alpha_library_complexity.pdf")

as.data.frame(dat  %>% summarize(
	mean_current_uniq_reads=mean(uniq_reads),
	median_current_uniq_reads=median(uniq_reads),
	mean_projected_uniq_reads=mean(projected_uniq_reads),
	median_projected_uniq_reads=median(projected_uniq_reads),
	cell_count=n()))

```

Currently not enough reads per cell to run.

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
library(scquantum)
register(MulticoreParam(progressbar = T, workers = 5), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/metACT/240115_RMMM_scalebiotest2/final_sc_bams")
#run in scalebio pipeline output: bamDeDup directory

dat <- runVarbin(".",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
dat <- findAneuploidCells(dat)

# Mark low-quality cells for filtering
dat <- findOutliers(dat)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat)


# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
dat  <- findClusters(dat,k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones

pdf("all_cells.subclone.umap.pdf")
plotUmap(dat, label = 'subclones')
dev.off()

pdf("all_cells.superclone.umap.pdf")
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')

dat$tn5<-unlist(lapply(strsplit(dat$sample,"[.]"),"[",2))
dat$tn5_plate<-substr(dat$tn5,1,1)
dat$tn5_row<-substr(dat$tn5,2,2)
dat$tn5_col<-substr(dat$tn5,3,4)

dat$cell_line<-"NA"
dat@colData[which(dat@colData$tn5_col %in% c("05","06","07","08","09","10")),]$cell_line<-"mda_mb_231"
dat@colData[which(dat@colData$tn5_col %in% c("01","02","03","04")),]$cell_line<-"mcf10a"
dat@colData[which(dat@colData$tn5_col %in% c("11","12")),]$cell_line<-"gm12878"

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones','tn5_col','cell_line','reads_total'),order='hclust')
dev.off()

pdf("all_cells.subclone.phylo.pdf")
plotPhylo(dat, label = 'subclones')
dev.off()

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(dat, label = c('superclones','subclones',),order='hclust')
dev.off()

saveRDS(dat,file="all_cells.scCNA.rds")

colData(dat)$info<- 'scalemet'

################################################
###Merge with standard act seq for comparison###
################################################
#scalebio dat

dat <- runVarbin(".",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

dat$tn5<-unlist(lapply(strsplit(dat$sample,"[.]"),"[",2))
dat$tn5_plate<-substr(dat$tn5,1,1)
dat$tn5_row<-substr(dat$tn5,2,2)
dat$tn5_col<-substr(dat$tn5,3,4)
dat$info <- 'scalebio'

dat$cell_line<-"NA"
dat@colData[which(dat@colData$tn5_col %in% c("05","06","07","08","09","10")),]$cell_line<-"mda_mb_231"
dat@colData[which(dat@colData$tn5_col %in% c("01","02","03","04")),]$cell_line<-"mcf10a"
dat@colData[which(dat@colData$tn5_col %in% c("11","12")),]$cell_line<-"gm12878"

#act set
act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
colData(act)$info <- 'act'
colData(act)$cell_line<-"mda_mb_231"

colData(dat)<-colData(dat)[colnames(colData(dat)) %in% colnames(colData(act))]
colData(act)<-colData(act)[colnames(colData(act)) %in% colnames(colData(dat))]

merged_copykit <- cbind(act, dat)

# Mark euploid cells if they exist
merged_copykit <- findAneuploidCells(merged_copykit) 
merged_copykit <- findOutliers(merged_copykit)
merged_copykit <- knnSmooth(merged_copykit)
merged_copykit <- runUmap(merged_copykit)
k_clones<-findSuggestedK(merged_copykit)
merged_copykit <- findClusters(merged_copykit ,k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK)#output from k_clones

pdf("act_merged.subclone.umap.pdf")
plotUmap(merged_copykit, label = 'subclones')
dev.off()

pdf("act_merged.cellline.umap.pdf")
plotUmap(merged_copykit, label = 'cell_line')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
merged_copykit <- calcConsensus(merged_copykit)
merged_copykit <- runConsensusPhylo(merged_copykit)
merged_copykit <- runPhylo(merged_copykit, metric = 'manhattan')

# Plot a copy number heatmap with clustering annotation
pdf("act_merged.subclone.heatmap.pdf")
plotHeatmap(merged_copykit, label = c('superclones','subclones','cell_line','info','reads_total'),order='hclust')
dev.off()

saveRDS(dat,file="act_merged.scCNA.rds")
```

Submit job
```bash
bsub -Is -W 36:00 -q long -n 10 -M 200 -R rusage[mem=200] scalemethyl_pipeline.sh
```

```bash
count_scbam() {
	echo "$1\t$(samtools view $1 | wc -l)"
}
export -f count_scbam
bam_in=$(find -iname "*.bam")
parallel --jobs 50 count_scbam ::: $bam_in | sort -k1,1n > cell_readcounts.tsv



gzip_small_bams() {
	gzip $1
}
export -f gzip_small_bams
ls -lS | awk '$5<30000 {print $9}' > small_file_list.txt
parallel --jobs 50 -a small_file_list.txt count_scbam

```
