## Settiing up the ScaleMethyl Nextflow pipeline to work on a IBM LSF HPC

Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On mdandersons seadragon2.

```bash
#ssh RMulqueen@seadragon2
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
mkdir -p ${proj_dir}

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
```

Transfer Our Data from our labs node 
```bash
#get data from cell line run onto seadragon
mkdir -p ${proj_dir}/240115_RMMM_scalebiotest2
cd ${proj_dir}/240115_RMMM_scalebiotest2
sftp mulqueen@qcprpgeo.mdanderson.edu
get -r /volumes/seq/flowcells/MDA/nextseq2000/2024/240111_VH00219_563_AAFFFCCM5 #grab full rundir
```

Pull all necessary files for scalemethyl nextflow pipeline
```bash
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
mkdir -p ~/singularity
cd ~/singularity
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8
singularity pull docker://nfcore/bclconvert:3.9.3

#these addresses blocked on seadragon for some reason, so just transferring in with sftp from geo (which has less restricted internet access)
sftp mulqueen@qcprpgeo.mdanderson.edu
cd ~/singularity
get scale_methyl_v1.5.sif
get scale_methyl_python_dependencies_v1.2.sif

#otherwise can use same pull commands
#singularity pull docker://public.ecr.aws/o5l3p3e4/scale_methyl_python_dependencies:v1.2
#singularity pull docker://public.ecr.aws/o5l3p3e4/scale_methyl:v1.5
```

## Running on seadragon2

Initial set up of files and running bcl-convert
```bash
bsub -Is -W 6:00 -q interactive -n 7 -M 100 -R rusage[mem=100] /bin/bash #small interactive node for the P1 300 cycle nextseq run of one column

module load nextflow/23.04.3

#set dir variables
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/240115_RMMM_scalebiotest2"
fastqDir="${proj_dir}/240115_RMMM_scalebiotest2/240111_VH00219_563_AAFFFCCM5"
samples="${proj_dir}/samples.csv"
genome="${proj_dir}/ref/genome.json"

cd $runDir

#ran first and found all the reads assigned to undetermined, then looked up resulting indexes to confirm indexes we expect
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

#spot check indexes
zcat Undetermined_S0_L001_I1_001.fastq.gz | sed -n '2~4p' - | head
zcat Undetermined_S0_L001_I2_001.fastq.gz | sed -n '2~4p' - | head 

#run bcl-convert
singularity exec --bind ${runDir} \
~/singularity/bclconvert_3.9.3.sif \
bcl-convert \
--bcl-input-directory ${fastqDir} \
--bcl-num-conversion-threads 1 \
--bcl-num-compression-threads 1 \
--bcl-num-decompression-threads 1 \
--output-directory ${runDir}/fq_out \
--sample-sheet ${runDir}/samplesheet.csv \
--force

#delete these from fq_out directory, otherwise they will throw errors in ScaleMethyl Pipeline (it tries to demux these)
rm -rf Undetermined*fastq.gz #write them out and then check size, then delete to save space

#Set up samples.csv for ScaleMethyl
echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-2H12,ScaleMethyl""" > ${runDir}/samples.csv

```

## Modify the ${scalebio_nf}/nextflow.config file to force offline singularity use.

```bash
#change the nextflow config file to use local SIF rather than external downloads
#format is file://<path from librarydir>

cp ${scalebio_nf}/nextflow.config ${scalebio_nf}/nextflow.config.backup #make backup

#export binding path for singularity, so home and project directory are always bound
echo "export SINGULARITY_BIND='$HOME,/rsrch4/home/genetics/rmulqueen/projects/metact'" >> ~/.bashrc
```

Remove these lines from ${scalebio_nf}/nextflow.config (i just did this with nano):
```bash
profiles {
  conda {
    conda.enabled = true
    process.conda = "$projectDir/envs/scaleMethyl.conda.yml"
    process {
      withLabel:py_process {
        conda = "$projectDir/envs/scaleMethylPyQc.yml"
      }
    }
  }
  docker {
     docker {
      enabled = true
    }
  }
  singularity {
    singularity.enabled = true 
    docker.enabled = false
  }
}
```

Add these lines to nextflow.config (i just did this with nano)

```bash
singularity.enabled = true
docker.enabled = false
singularity.autoMounts = true
singularity.noHttps = true
```

And replace the remote pulls of singularity/docker files with our local ones.
```bash
#"public.ecr.aws/o5l3p3e4/scale_methyl:v1.5"
sed -i 's/public.ecr.aws\/o5l3p3e4\/scale_methyl\:v1.5/file:\/\/\/rsrch4\/home\/genetics\/rmulqueen\/singularity\/scale_methyl_v1.5.sif/g' ${scalebio_nf}/nextflow.config

#"public.ecr.aws/o5l3p3e4/scale_methyl_python_dependencies:v1.2"
sed -i 's/public.ecr.aws\/o5l3p3e4\/scale_methyl_python_dependencies\:v1.2/file:\/\/\/rsrch4\/home\/genetics\/rmulqueen\/singularity\/scale_methyl_python_dependencies_v1.2.sif/g' ${scalebio_nf}/nextflow.config

#'nfcore/bclconvert:3.9.3'
sed -i 's/nfcore\/bclconvert\:3.9.3/file:\/\/\/rsrch4\/home\/genetics\/rmulqueen\/singularity\/bclconvert_3.9.3.sif/g' ${scalebio_nf}/nextflow.config

#'biocontainers/fastqc:v0.11.9_cv8'
sed -i 's/biocontainers\/fastqc\:v0.11.9_cv8/file:\/\/\/rsrch4\/home\/genetics\/rmulqueen\/singularity\/fastqc_v0.11.9_cv8.sif/g' ${scalebio_nf}/nextflow.config
```

## Now to submit to a seadragon2 node.

Consider max mem per queue
http://hpcweb/using_the_cluster.dir/lsf.dir/MaxMemoryChart.pdf
and max cores and wall time
http://hpcweb/using_the_cluster.dir/lsf.dir/QueueParameters_2023-12-08.pdf

Generate the 240115_scalemet.lsf file for bsub submission.

```bash
#BSUB -W 6:00
#BSUB -q e40medium
#BSUB -n 40
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -o /rsrch4/home/genetics/rmulqueen/projects/metact/240115_RMMM_scalebiotest2/bsub.log
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/projects/metact/240115_RMMM_scalebiotest2
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240115_scalemet

#load modules
module load nextflow/23.04.3
module load singularity/3.7.0
module load glib/2.51.5 #might not need this

#set up directories and variables
export proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
export scalebio_nf="${proj_dir}/tools/ScaleMethyl"
export runDir="${proj_dir}/240115_RMMM_scalebiotest2"
export fastqDir="${proj_dir}/240115_RMMM_scalebiotest2/240111_VH00219_563_AAFFFCCM5"
export samples="${proj_dir}/samples.csv"
export genome="${proj_dir}/ref/genome.json"

#set global cache and library directories
export NXF_SINGULARITY_LIBRARYDIR="/rsrch4/home/genetics/rmulqueen/singularity" #might not need these
export NXF_SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity" #might not need these

#--profile command can also probably be removed since we forced the nextflow.config file
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fq_out \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--genome ${genome} \
--fastqOut true \
--trimOut true \
--bamOut true \
--bamDedupOut true \
--covOut true \
--max_memory 100.GB \
--max_cpus 40 \
--profile docker,singularity \
-w ${runDir}/work \
-resume

```

## Submit job.

```bash
cd /rsrch4/home/genetics/rmulqueen/projects/metact/240115_RMMM_scalebiotest2/
rm -rf bsub.log
bsub < 240115_scalemet.lsf
```

## Move back to our local HPC

```bash
bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

cd /rsrch4/home/genetics/rmulqueen/projects/metact

#mostly just avoiding the big flowcell data
rsync -alPvz \
./240115_RMMM_scalebiotest2/bamDeDup \
./240115_RMMM_scalebiotest2/bsbolt \
./240115_RMMM_scalebiotest2/cg_sort_cov \
./240115_RMMM_scalebiotest2/ch_sort_cov \
./240115_RMMM_scalebiotest2/fastq \
./240115_RMMM_scalebiotest2/fastqc \
./240115_RMMM_scalebiotest2/fq_out \
./240115_RMMM_scalebiotest2/library_barcode_metrics \
./240115_RMMM_scalebiotest2/library_report \
./240115_RMMM_scalebiotest2/Logs \
./240115_RMMM_scalebiotest2/matrix \
./240115_RMMM_scalebiotest2/report \
./240115_RMMM_scalebiotest2/Reports \
./240115_RMMM_scalebiotest2/*csv \
./240115_RMMM_scalebiotest2/*lsf \
./240115_RMMM_scalebiotest2/*json \
mulqueen@qcprpgeo.mdanderson.edu:/volumes/seq/projects/metACT/240115_RMMM_scalebiotest2

cd /rsrch4/home/genetics/rmulqueen/projects/metact

rsync -alPvz \
./240115_RMMM_scalebiotest2/work \
mulqueen@qcprpgeo.mdanderson.edu:/volumes/seq/projects/metACT/240115_RMMM_scalebiotest2

```

# Continued work on qcprpgeo.mdanderson.edu

Link local run directory to working directory for posterity
```bash
proj_dir="/volumes/seq/projects/metACT"
runDir="${proj_dir}/240115_RMMM_scalebiotest2"
ln -s /volumes/seq/flowcells/MDA/nextseq2000/2024/240111_VH00219_563_AAFFFCCM5 ${runDir}/.
```

Fix all the symlink outputs with our local server address
Both seadragon and geo house 
```bash
proj_dir="/volumes/seq/projects/metACT"
runDir="${proj_dir}/240115_RMMM_scalebiotest2"

find $runDir -xtype l -exec bash -c 'target="$(readlink "{}")"; link="{}"; target="$(echo "$target" | sed "s/\/rsrch4\/home\/genetics\/rmulqueen\/projects\/metact/\/volumes\/seq\/projects\/metACT/g")"; ln -Tfs "$target" "$link"' \;
```