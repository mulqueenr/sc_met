Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On seadragon.

```bash
#ssh RMulqueen@seadragon
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
mkdir ${proj_dir}

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

bsub -Is -W 4:00 -q interactive -n 1 -M 16 -R rusage[mem=16] /bin/bash
```

Set up nextflow and dependencies
```bash
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
mkdir ${proj_dir}/tools
#clone git
git clone "https://github.com/ScaleBio/ScaleMethyl.git"

mkdir ${proj_dir}/ref
#download prebuilt human genome reference
wget http://scale.pub.s3.amazonaws.com/genomes/grch38.tgz
tar xvfz grch38.tgz

module load awscli
aws configure #set up access key and secrete id from AWS account

#install bsbolt
conda install -c cpfarrell bsbolt #installed to base conda env

mkdir tmp
TMPDIR="~/tmp"
EXPORT $TMPDIR
```

Download Test data and Reference Folder
```bash
module load awscli
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"

mkdir ${proj_dir}/test_run
aws s3 cp s3://scale.pub/testData/methylation . --recursive
cp -R ${proj_dir}/test_run/reference ${proj_dir}/ref
#set up their test run references with my own genome.json
```

Pull singularity files
```bash
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8
singularity pull docker://nfcore/bclconvert:3.9.3
#these are blocked on seadragon for some reason, so just transferring in with sftp
singularity pull docker://public.ecr.aws/o5l3p3e4/scale_methyl_python_dependencies:v1.2 #
singularity pull docker://public.ecr.aws/o5l3p3e4/scale_methyl:v1.5
```

Transfer Our Data
```bash
#get data from cell line run onto seadragon
mkdir ${proj_dir}/20231202_RMMM_scalebio_celllines
cd ${proj_dir}/20231202_RMMM_scalebio_celllines
sftp mulqueen@10.132.80.157
get /volumes/seq/flowcells/MDA/nextseq2000/2023/20231202_RMMM_scalebio_celllines/*
/volumes/seq/flowcells/MDA/nextseq2000/2023/20231202_RMMM_scalebio_celllines
```

Running the testrun on geo

```bash

#run on seadragon transfer and manually pull singularity files (just ran nextflow below and copied singularity pull requests from stderr)

cd ~/singularity
module load singularity/3.7.0
module load nextflow/22.10.1

proj_dir="/volumes/seq/projects/metACT"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
fastqDir="${proj_dir}/scalebio_test_run_data/pbmc_downsampled_outdir/downsampled_pbmcs"
samples="${proj_dir}/scalebio_test_run_data/pbmc_downsampled_outdir/samples.csv"
genome="${proj_dir}/ref/genome.json"
cd $fastqDir


#these make it run locally, and use the singularity containers we pulled manually above
export NXF_SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity/"
export TMPDIR="/rsrch4/home/genetics/rmulqueen/tmp"
export SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity/"
#mkdir -p $SINGULARITY_CACHEDIR #make sure these directories are made
#mkdir -p $SINGULARITY_CACHEDIR/tmp
#mkdir -p $SINGULARITY_CACHEDIR/pull
export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
export SINGULARITY_PULLDIR=$SINGULARITY_CACHEDIR/pull
export CWL_SINGULARITY_CACHE=$SINGULARITY_PULLDIR
export NXF_OFFLINE='TRUE' #https://nf-co.re/docs/usage/offline

cd $fastqDir
nextflow run ${scalebio_nf}/main.nf \
--fastqDir=$fastqDir \
--samples=$samples \
--genome=$genome \
-profile docker,singularity


#nano /volumes/seq/projects/metACT/tools/ScaleMethyl/main.nf 

```

Running our samples
```bash
#changed i5.txt to i5.cellline.txt (reverse compliments)
#changed i7.txt to i7.cellline.txt (reverse compliments)
#changed lib.json in /volumes/seq/projects/metACT/tools/ScaleMethyl/references/
#"index2Seqs": "i5.cellline.txt" etc

#current wd for figuring out whats going on with the indexes
proj_dir="/volumes/seq/projects/metACT"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5"
samples="${proj_dir}/231204_scalebio_test/samples.csv"
genome="${proj_dir}/ref/genome.json"
workingDir="${proj_dir}/231204_scalebio_test"

#run bcl2fastq myself 
flowcelldir="${proj_dir}/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5"
fq_out="${proj_dir}/231204_scalebio_test/fastq"

cd $workingDir

#these make it run locally, and use the singularity containers we pulled manually above
#export NXF_SINGULARITY_CACHEDIR="~/singularity/"
#export TMPDIR="/~/tmp"
#export SINGULARITY_CACHEDIR="~/singularity/"
#mkdir -p $SINGULARITY_CACHEDIR #make sure these directories are made
#mkdir -p $SINGULARITY_CACHEDIR/tmp
#mkdir -p $SINGULARITY_CACHEDIR/pull
#export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
#export SINGULARITY_PULLDIR=$SINGULsARITY_CACHEDIR/pull
#export CWL_SINGULARITY_CACHE=$SINGULARITY_PULLDIR
#export NXF_OFFLINE='TRUE' #https://nf-co.re/docs/usage/offline

#run bcl2fastq myself 
flowcelldir="/volumes/seq/projects/metACT/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5"
fq_out="/volumes/seq/projects/metACT/231204_scalebio_test/fastq"

#running without sample sheet first. Then running to assign undetermined reads.
bcl2fastq --sample-sheet samplesheet.csv \
-R 231202_VH00219_541_AAF3KVTM5 \
--no-lane-splitting \
-o fastq

samplesheet="/volumes/seq/projects/metACT/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5/samplesheet.csv"

#Get count of indexes
zcat Undetermined_S0_L001_I1_001.fastq.gz | sed -n '2~4p' | sort -T . | uniq -c | sort -k1,1n > i1_counts.txt
tail i1_counts.txt 
#  20503 CCTAGCGGTA
#  21072 CATAAGCGGT
#  30636 TGTTATGGAA
#  39070 NNNNNNNNNN
#  40746 TTGATATGGA
#  41823 TTTATATGAA
#  42583 TTTGATAGGA
#12007894 GGGGGGGGGG
#25933596 CCTAAGCGGT
#86347587 TTGATATGAA
#presumed 1 index for i1: TTGATATGAA

zcat Undetermined_S0_L001_I2_001.fastq.gz | sed -n '2~4p' | sort -T . | uniq -c | sort -k1,1n > i2_counts.txt
tail i2_counts.txt 
#2281570 GGGGGGGGGG
#3937031 CGTCGCAGCC
#9476537 CGGCGTAACT
#11112786 TAACGATCCA
#11674295 AGATCTCGGT
#13110316 TATCATGATC
#14440971 TTGCCTTGGC
#15614202 GTAGCTCCAT
#15890159 TGCTAATTCT
#25113911 GAGCATATGG
#Presumed 8 indexes for i2: CGGCGTAACT TAACGATCCA AGATCTCGGT TATCATGATC TTGCCTTGGC GTAGCTCCAT TGCTAATTCT GAGCATATGG

echo """[Settings]
CreateFastqForIndexReads,1
MinimumTrimmedReadLength,16
MaskShortReads,16
OverrideCycles,Y150;I10;I10;Y150
[Data]
Sample_ID,index,index2
cellline_I01,TTGATATGAA,CGGCGTAACT
cellline_I01,TTGATATGAA,TAACGATCCA
cellline_I01,TTGATATGAA,AGATCTCGGT
cellline_I01,TTGATATGAA,TATCATGATC
cellline_I01,TTGATATGAA,TTGCCTTGGC 
cellline_I01,TTGATATGAA,GTAGCTCCAT
cellline_I01,TTGATATGAA,TGCTAATTCT
cellline_I01,TTGATATGAA,GAGCATATGG""" > ${flowcelldir}/samplesheet.csv

bcl2fastq \
-R $flowcelldir \
--create-fastq-for-index-reads \
-o $fq_out \
--sample-sheet $samplesheet

mkdir Undetermined
mv Undetermined* Undetermined

#current wd for figuring out whats going on with the indexes
proj_dir="/volumes/seq/projects/metACT"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
fq_out="${proj_dir}/231204_scalebio_test/fastq"
samples="${proj_dir}/231204_scalebio_test/samples.csv"
runDir="${proj_dir}/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5"
genome="${proj_dir}/ref/genome.json"
workingDir="${proj_dir}/231204_scalebio_test"
flowcelldir="${proj_dir}/231204_scalebio_test/231202_VH00219_541_AAF3KVTM5"
export SINGULARITY_CACHEDIR="~/singularity/"
mkdir -p $SINGULARITY_CACHEDIR #make sure these directories are made

cd ${workingDir}



#conda install
bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
conda install mamba
mamba update -n base -c defaults conda

proj_dir="/volumes/seq/projects/metACT"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
bash ${scalebio_nf}/envs/download-scale-tools.sh 
mamba env create -f ${scalebio_nf}/envs/scaleMethyl.conda.yml 
mamba env create -f ${scalebio_nf}/envs/scaleMethylPyQc.yml 
#All tools need to be available on $PATH or in /PATH/TO/ScaleMethyl/bin/

/volumes/seq/projects/metACT/tools/ScaleMethyl/bin

#run til failure (combined report)
nextflow run ${scalebio_nf}/main.nf --fastqDir ${fq_out} --samples ${samples} --outDir "${fq_out}/231204_ScaleMethyl" --genome ${genome} --fastqOut true --trimOut true --bamOut true --bamDedupOut true --covOut true --max_memory 100.GB --max_cpus 50 -profile docker,singularity -resume

#it ends up stalling on the combine report due to some error. I think I'll prepull everything and run on exacloud via conda in the future instead.
```