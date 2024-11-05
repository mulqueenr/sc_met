## Setting up the ScaleMethyl Nextflow pipeline to work on a IBM LSF HPC

Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On mdandersons seadragon.

Transfer Our Data to our lab directory on seadragon
```bash
#get data from cell line run onto seadragon
proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
runDir="${proj_dir}/241007_RM_scalebio_dcis2"
mkdir -p ${runDir}
cd ${runDir}
#novaseq S2 flowcell run by sequencing core temporarily stored here: 
#/rsrch4/home/genetics/htep/PI/nnavin/241004_A01819_0637_BHY5MJDMXY/
mkdir ${runDir}/241004_A01819_0637_BHY5MJDMXY
cp -R /rsrch4/home/genetics/htep/PI/nnavin/241004_A01819_0637_BHY5MJDMXY/ ${runDir}/241004_A01819_0637_BHY5MJDMXY

###########transfer to geo##############


#used i5 index set 1 (plate 1) and set 2 (plate 2) and i7 index set 1
#Set up samples.csv for ScaleMethyl
echo """sample,barcodes,libName
DCIS-92T,1A01-1D12,ScaleMethyl
DCIS-66T,1E01-1H12,ScaleMethyl
DCIS-79T,2A01-2D12,ScaleMethyl
IDC-79T,2E01-2H12,ScaleMethyl
HBCA-19T,3A01-3D12,ScaleMethyl
HBCA-17T,3E01-3H12,ScaleMethyl""" > ${runDir}/samples.csv

#build proper formated singularity container
singularity build ~/singularity/scalemethyl_v1.6.sif ~/singularity/public.ecr.aws-o5l3p3e4-scale-methyl-tools@sha256-6fd63db48e8786ed1cfc17d7e3effd3fd696ccb8e5e54803959e2dcd2f794aec.img

#set up directories and variables
proj_dir="/volumes/USR2/Ryan/projects/metact"
scalebio_nf="${proj_dir}/tools2/ScaleMethyl" #tools is a symlink directory so it wasn't mounting properly for singularity
runDir="${proj_dir}/241007_RM_scalebio_dcis2"
genome="${proj_dir}/ref/reference/genome.json"
fastqDir="${runDir}/241004_A01819_0637_BHY5MJDMXY/241004_A01819_0637_BHY5MJDMXY"
samples="${runDir}/samples.csv"

export SCRATCH="/volumes/USR2/Ryan/scratch/scalemet_work"
export TMPDIR="/volumes/USR2/Ryan/scratch"
export NXF_SINGULARITY_CACHEDIR="/volumes/USR2/Ryan/singularity"
export SINGULARITY_BINDPATH="/volumes/seq/projects/metACT/tools/ScaleMethyl/bin" 

mkdir -p $SCRATCH

source activate conda #(to use more recent java version)
cd $runDir
nextflow run ${scalebio_nf} \
--runFolder ${fastqDir} \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--genome ${genome} \
--maxMemory 500.GB \
--maxCpus 300 \
-profile singularity \
-w ${SCRATCH}/scalemet_work \
-resume


```
```bash
mkdir -p transfer_dat
mkdir -p transfer_dat_bam
find ./cg_sort_cov -type l -exec bash -c 'cp -R "$(readlink -m "$0")" ./transfer_dat' {} \; #scale pipeline makes empty files, this throws errors for empty files (but can be ignored)
find ./bamDeDup -type l -exec bash -c 'cp -R "$(readlink -m "$0")" ./transfer_dat_bam' {} \; #scale pipeline makes empty files, this throws errors for empty files (but can be ignored)

bsub -Is -W 6:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
sftp mulqueen@qcprpgeo.mdanderson.edu

cp -R ./report ./transfer_dat
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240526_RMMM_scalebio_dcis/cg_sort_cov
```

#split bam (on geo also)

```bash
cd /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup

#get read counts of cellIDs, merge all tagmentation well read counts together, and then run scbam splitting
readcount_bam() {
	bam_name=$(basename $1)
    bam_path=$1
	export cellline=$(echo $bam_name | cut -d '.' -f 1 )
	export well=$(echo $bam_name | cut -d '.' -f 2 )
	samtools view $1 \
	| awk -v b=${bam_path} '{split($1,a,":"); print a[8],b}' \
	| sort \
	| uniq -c \
	| sort -k1,1n \
	| awk '$1>100000 {print $0}' > ${cellline}.${well}.tgmt_readcount.tsv
}

export -f readcount_bam
parallel -j 50 readcount_bam ::: $(ls -R ./*/*bam)

cat *tgmt_readcount.tsv > readcount.tsv

split_bam() {
test=$1
idx=$(echo $test | cut -d ' ' -f 2 )
outidx=$(echo $idx | sed -e 's/+/_/g' -)
bam=$(echo $test | cut -d ' ' -f 3)
cellline=$(basename $bam | cut -d '.' -f 1 )
well=$(basename $bam | cut -d '.' -f 2 )
mkdir -p ./${cellline}_split_bam
((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split($1,a,":"); if(a[8]==i) print $0}')) \
| samtools view -bS - \
| samtools sort -T . -O BAM -o ./${cellline}_split_bam/${cellline}.${well}.${outidx}.sorted.bam - 
}

export -f split_bam
parallel -j 100 -a readcount.tsv split_bam 

#index bams for igv
index_bam() {
outname=$1
samtools index -b -@ 1 $1 ${1}.bai
}

export -f index_bam
parallel -j 100 index_bam ::: $(ls -R ./*split_bam/*bam)

"""

# Run CNV calling
```bash
singularity shell \
--bind ~/projects \
~/singularity/copykit.sif

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./DCIS-66T_split_bam/ -p DCIS-66T -c 100

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./DCIS-92T_split_bam/ -p DCIS-92T -c 100

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./DCIS-79T_split_bam/ -p DCIS-79T -c 100

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./IDC-79T_split_bam/ -p IDC-79T -c 100

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./HBCA-17T_split_bam/ -p HBCA-17T -c 100

Rscript /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/copykit_run_scalemet.R \
-i ./HBCA-19T_split_bam/ -p HBCA-19T -c 100

```

```R

library(copykit)
library(BiocParallel)
BiocParallel::bpparam()
library(optparse)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="allcells", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=1, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=100
prefix="DCIS-92T"
setwd("./DCIS-92T_split_bam/")

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

dat <- runVarbin(".",
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE)


dat$cell_line<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",1))
dat$well<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",2))
dat$idx<-unlist(lapply(strsplit(dat$sample,"[.]"),"[[",3))

# Mark euploid cells if they exist
dat <- findAneuploidCells(dat)

# Mark low-quality cells for filtering
dat <- findOutliers(dat)
#dat<- dat[,colData(dat)$outlier == FALSE]

pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat,k=10)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".superclone.umap.pdf"))
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')

# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, label = c('superclones','subclones','cell_line','reads_total') ,order_cells='consensus_tree',n_threads=50)
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)
```