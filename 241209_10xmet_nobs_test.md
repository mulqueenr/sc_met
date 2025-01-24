Experiment: https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%2810xmet.one%7C6C1BDECA-4977-0A4C-AD58-1D18D92E3B99%2F241202%2010xMET%20ATAC%20Kit%205%27%20Splint%7CBCF23CDF-B6C8-B349-8BF2-4F98B9264C0E%2F%29

Sequencing folder: /volumes/seq/flowcells/MDA/nextseq2000/2024/20241207_Ryan_10xMET_WGD/241207_VH01788_81_AAG7V2TM5

Ran 10X ATAC kit with demetCTP addition. Using their supplied Tn5 and following the protocol as they describe. This was a control for the first 5' Splint Attempt. The 5' Split failed, but the control will give us information on complexity we can expect.

Use 10x cellranger for initial fastq generation

```bash
#set up environment variables 
export SCRATCH="/volumes/USR2/Ryan/work"
export projDir="/volumes/USR2/Ryan/projects/10x_MET"
export srcDir="${projDir}/src"
export refDir="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
export outDir="${projDir}/241105_10xmet_231_nobs"
export DNArunDir="/volumes/seq/flowcells/MDA/nextseq2000/2024/20241207_Ryan_10xMET_WGD/241207_VH01788_81_AAG7V2TM5"
export outname="10xmet_231"
export task_cpus=250
export index="/volumes/USR2/Ryan/projects/metact/ref/bsbolt/genome_index.pkl"

mkdir -p $outDir
cd $outDir

echo """Lane,Sample,Index
*,${outname}_noconv,SI-NA-C9""" > sample_met.csv
export sample_met="${outDir}/sample_met.csv"

mkdir -p ${SCRATCH}/10xmet_work

~/tools/cellranger-atac-2.1.0/cellranger-atac mkfastq \
--run=${DNArunDir} \
--csv=${sample_met} \
--jobmode="local" \
--localcores=${task_cpus} \
--delete-undetermined \
--output-dir="${outDir}" \
--id=${outname}_met

#can use regular cell ranger for this
~/tools/cellranger-atac-2.1.0/cellranger-atac count \
                        --id=${outname}_met \
                        --reference=${refDir} \
                        --fastqs="${outDir}/${outname}_met/outs/fastq_path/AAG7V2TM5/10xmet_231_noconv" \
                        --localcores=100 \
                        --localmem=200


#split out single cells to estimate complexity

singularity shell \
--bind /volumes/USR2/Ryan/projects/metact/ref/:/ref \
~/singularity/amethyst.sif
source ~/.bashrc

bam_in="/volumes/USR2/Ryan/projects/10x_MET/241105_10xmet_231_nobs/10xmet_231_met/outs/fastq_path/AAG7V2TM5/10xmet_231_noconv/10xmet_231_met/outs/possorted_bam.bam"

java -jar ~/tools/picard.jar EstimateLibraryComplexity \
--I ${bam_in} \
--O "allcells.metrics.txt" \
--BARCODE_TAG CB:Z

samtools split -d CB \
-u ${outDir}/unassigned_reads.bam \
-M 6000 --output-fmt BAM \
-f ${outDir}'/sc_bam/%*_%!.%.' ${bam_in}

#remove files with size less than 3M
find . -maxdepth 1 -type f -name "*.bam" -size -3000k -delete

#* Based on the Lander-Waterman equation that states:
#*     C/X = 1 - exp( -N/X )
#* where
#*     X = number of distinct molecules in library
#*     N = number of read pairs
#*     C = number of distinct fragments observed in read pairs

lib_complexity(){
    bam=$1
    samtools sort -n $bam | samtools fixmate -p -m - - | samtools sort | samtools markdup -r -S  --write-index - ${bam::-4}.bbrd.bam
    java -jar ~/tools/picard.jar EstimateLibraryComplexity \
        I="${bam::-4}.bbrd.bam" \
        O="${bam::-4}.metrics.txt"
    }

export -f lib_complexity
parallel -j 50 lib_complexity ::: $(ls *-1.bam)

for i in *metrics.txt; do echo -e "$i \t $(grep "^MissingLibrary.1" $i | awk 'OFS="\t" {print $3,$9,$10}')"; done > read_counts.tsv
```

```R
library(ggplot2)
dat<-read.table("read_counts.tsv",sep="\t")
setwd("/volumes/USR2/Ryan/projects/10x_MET/241105_10xmet_231_nobs")
colnames(dat)<-c("cellID","total_reads","duplication_rate","estimated_libsize")
dat$duplication_perc<-dat$duplication_rate*100
dat$uniq_perc<-100-dat$duplication_perc

plt<-ggplot(dat,aes(y=log10(total_reads),x=uniq_perc))+geom_point()+theme_minimal()+coord_cartesian(xlim=c(0,100),ylim=c(0,8))
ggsave(plt,file="complexity_plot.pdf")

plt<-ggplot(dat,aes(y=log10(estimated_libsize),x=1))+geom_violin()+geom_jitter()+theme_minimal()+coord_cartesian(ylim=c(0,8))
ggsave(plt,file="estimated_size.pdf")
```

```bash
singularity shell \
--bind ~/projects \
~/singularity/copykit.sif
source ~/.bashrc
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
setwd("/volumes/USR2/Ryan/projects/wgd/241210_wgd/data/cells")
cpu_count=100
prefix="241210"

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

dat <- runVarbin(".",
    remove_Y = TRUE,
    genome="hg38",
    is_paired_end=TRUE, min_bincount=0)

dat  <- runMetrics(dat)

#dat <- dat[,colData(dat)$reads_total > 100000]

pdf(paste0(prefix,".qc_metrics.pdf"))
plotMetrics(dat, metric = c("overdispersion", 
                              "breakpoint_count",
                              "reads_total",
                              "reads_duplicates",
                              "reads_assigned_bins",
                              "percentage_duplicates"),
            label = "reads_total")
dev.off()

# Mark euploid cells if they exist
dat <- findAneuploidCells(dat)

# Mark low-quality cells for filtering
dat <- findOutliers(dat)


pdf(paste0(prefix,".outlier_qc.heatmap.pdf"))
plotHeatmap(dat, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat,k=10)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones
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
plotHeatmap(dat, label = c('reads_total'),order='hclust')
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)
```