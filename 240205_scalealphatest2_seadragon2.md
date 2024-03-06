## Setting up the ScaleMethyl Nextflow pipeline to work on a IBM LSF HPC

Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On mdandersons seadragon2.

```bash
#ssh RMulqueen@seadragon2
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
mkdir -p ${proj_dir}

bsub -Is -W 6:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
```

Transfer Our Data from our labs node 
```bash
#get data from cell line run onto seadragon
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
mkdir -p ${proj_dir}/240205_RMMM_scalebiotest2
cd ${proj_dir}/240205_RMMM_scalebiotest2
sftp mulqueen@qcprpgeo.mdanderson.edu
get -r /volumes/seq/flowcells/MDA/NovaSeq/2024/PM2517 #grab full rundir
```

## Running on seadragon2
See 240115_scalealphatest2_seadragon2 for more details about singularity setup
Initial set up of files and running bcl-convert.
Note: I've been running this first in an interactive node to confirm the indexes off the first bit of reads. Then putting it onto a larger node with more cores.

```bash
bsub -Is -W 6:00 -q interactive -n 7 -M 100 -R rusage[mem=100] /bin/bash #small interactive node for bcl-convert 

module load nextflow/23.04.3

#set dir variables
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/240205_RMMM_scalebiotest2"
fastqDir="${proj_dir}/240205_RMMM_scalebiotest2/PM2517"
samples="${proj_dir}/samples.csv"
genome="${proj_dir}/ref/genome.json"

cd $runDir

#used i5 index set 1 and i7 index set 1
#ran first and found all the reads assigned to undetermined, then looked up resulting indexes to confirm indexes we expect
#use i7.txt and i5.txt for making sample sheet
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,1
TrimUMI,0
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,2
MinimumTrimmedReadLength,16
MaskShortReads,16
OverrideCycles,Y159;I10;I10;U8Y151
[BCLConvert_Data]
Sample_ID,index,index2
ScaleMethyl_GTAGCTCCAT,GGAGGCCTCC,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CAGCAGTATC,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CGAAGGCATG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TTCAATATAA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CGAATCTCCT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,GTCTCCAGAG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TCGAAGCGCG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,GCCGTCGCGT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TATTCCGTTA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,ACTGGTAGAT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,AATTGAGAGA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,ATGAGTTCTC,GTAGCTCCAT
ScaleMethyl_GGCAATGAGA,GGAGGCCTCC,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CAGCAGTATC,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CGAAGGCATG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TTCAATATAA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CGAATCTCCT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,GTCTCCAGAG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TCGAAGCGCG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,GCCGTCGCGT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TATTCCGTTA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,ACTGGTAGAT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,AATTGAGAGA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,ATGAGTTCTC,GGCAATGAGA
ScaleMethyl_CGGTTATGCC,GGAGGCCTCC,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CAGCAGTATC,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CGAAGGCATG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TTCAATATAA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CGAATCTCCT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,GTCTCCAGAG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TCGAAGCGCG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,GCCGTCGCGT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TATTCCGTTA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,ACTGGTAGAT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,AATTGAGAGA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,ATGAGTTCTC,CGGTTATGCC
ScaleMethyl_ACCGGAATTA,GGAGGCCTCC,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CAGCAGTATC,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CGAAGGCATG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TTCAATATAA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CGAATCTCCT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,GTCTCCAGAG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TCGAAGCGCG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,GCCGTCGCGT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TATTCCGTTA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,ACTGGTAGAT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,AATTGAGAGA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,ATGAGTTCTC,ACCGGAATTA
ScaleMethyl_CTTCGGCGCA,GGAGGCCTCC,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CAGCAGTATC,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CGAAGGCATG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TTCAATATAA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CGAATCTCCT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,GTCTCCAGAG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TCGAAGCGCG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,GCCGTCGCGT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TATTCCGTTA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,ACTGGTAGAT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,AATTGAGAGA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,ATGAGTTCTC,CTTCGGCGCA
ScaleMethyl_AATACGCAGG,GGAGGCCTCC,AATACGCAGG
ScaleMethyl_AATACGCAGG,CAGCAGTATC,AATACGCAGG
ScaleMethyl_AATACGCAGG,CGAAGGCATG,AATACGCAGG
ScaleMethyl_AATACGCAGG,TTCAATATAA,AATACGCAGG
ScaleMethyl_AATACGCAGG,CGAATCTCCT,AATACGCAGG
ScaleMethyl_AATACGCAGG,GTCTCCAGAG,AATACGCAGG
ScaleMethyl_AATACGCAGG,TCGAAGCGCG,AATACGCAGG
ScaleMethyl_AATACGCAGG,GCCGTCGCGT,AATACGCAGG
ScaleMethyl_AATACGCAGG,TATTCCGTTA,AATACGCAGG
ScaleMethyl_AATACGCAGG,ACTGGTAGAT,AATACGCAGG
ScaleMethyl_AATACGCAGG,AATTGAGAGA,AATACGCAGG
ScaleMethyl_AATACGCAGG,ATGAGTTCTC,AATACGCAGG
ScaleMethyl_TAACTCTTAG,GGAGGCCTCC,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CAGCAGTATC,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CGAAGGCATG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TTCAATATAA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CGAATCTCCT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,GTCTCCAGAG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TCGAAGCGCG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,GCCGTCGCGT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TATTCCGTTA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,ACTGGTAGAT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,AATTGAGAGA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,ATGAGTTCTC,TAACTCTTAG
ScaleMethyl_TCGTAGGCTT,GGAGGCCTCC,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CAGCAGTATC,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CGAAGGCATG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TTCAATATAA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CGAATCTCCT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,GTCTCCAGAG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TCGAAGCGCG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,GCCGTCGCGT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TATTCCGTTA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,ACTGGTAGAT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,AATTGAGAGA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,ATGAGTTCTC,TCGTAGGCTT""" > ${runDir}/samplesheet.csv


#spot check indexes
#zcat Undetermined_S0_L001_I1_001.fastq.gz | sed -n '2~4p' - | head #
#zcat Undetermined_S0_L001_I2_001.fastq.gz | sed -n '2~4p' - | head #i5 looks right (matches with https://github.com/ScaleBio/ScaleMethyl/blob/main/references/i5.txt)

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
rm -rf ${runDir}/fq_out/Undetermined*fastq.gz #write them out and then check size, then delete to save space

#Set up samples.csv for ScaleMethyl
echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-2H12,ScaleMethyl""" > ${runDir}/samples.csv

```




## Set up sample sheets
```bash
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/240205_RMMM_scalebiotest2"


#used i5 index set 1 and i7 index set 1
#ran first and found all the reads assigned to undetermined, then looked up resulting indexes to confirm indexes we expect
#use i7.txt and i5.txt for making sample sheet
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,1
TrimUMI,0
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,2
MinimumTrimmedReadLength,16
MaskShortReads,16
OverrideCycles,Y159;I10;I10;U8Y151
[BCLConvert_Data]
Sample_ID,index,index2
ScaleMethyl_GTAGCTCCAT,GGAGGCCTCC,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CAGCAGTATC,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CGAAGGCATG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TTCAATATAA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,CGAATCTCCT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,GTCTCCAGAG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TCGAAGCGCG,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,GCCGTCGCGT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,TATTCCGTTA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,ACTGGTAGAT,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,AATTGAGAGA,GTAGCTCCAT
ScaleMethyl_GTAGCTCCAT,ATGAGTTCTC,GTAGCTCCAT
ScaleMethyl_GGCAATGAGA,GGAGGCCTCC,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CAGCAGTATC,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CGAAGGCATG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TTCAATATAA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,CGAATCTCCT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,GTCTCCAGAG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TCGAAGCGCG,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,GCCGTCGCGT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,TATTCCGTTA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,ACTGGTAGAT,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,AATTGAGAGA,GGCAATGAGA
ScaleMethyl_GGCAATGAGA,ATGAGTTCTC,GGCAATGAGA
ScaleMethyl_CGGTTATGCC,GGAGGCCTCC,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CAGCAGTATC,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CGAAGGCATG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TTCAATATAA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,CGAATCTCCT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,GTCTCCAGAG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TCGAAGCGCG,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,GCCGTCGCGT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,TATTCCGTTA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,ACTGGTAGAT,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,AATTGAGAGA,CGGTTATGCC
ScaleMethyl_CGGTTATGCC,ATGAGTTCTC,CGGTTATGCC
ScaleMethyl_ACCGGAATTA,GGAGGCCTCC,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CAGCAGTATC,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CGAAGGCATG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TTCAATATAA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,CGAATCTCCT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,GTCTCCAGAG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TCGAAGCGCG,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,GCCGTCGCGT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,TATTCCGTTA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,ACTGGTAGAT,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,AATTGAGAGA,ACCGGAATTA
ScaleMethyl_ACCGGAATTA,ATGAGTTCTC,ACCGGAATTA
ScaleMethyl_CTTCGGCGCA,GGAGGCCTCC,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CAGCAGTATC,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CGAAGGCATG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TTCAATATAA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,CGAATCTCCT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,GTCTCCAGAG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TCGAAGCGCG,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,GCCGTCGCGT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,TATTCCGTTA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,ACTGGTAGAT,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,AATTGAGAGA,CTTCGGCGCA
ScaleMethyl_CTTCGGCGCA,ATGAGTTCTC,CTTCGGCGCA
ScaleMethyl_AATACGCAGG,GGAGGCCTCC,AATACGCAGG
ScaleMethyl_AATACGCAGG,CAGCAGTATC,AATACGCAGG
ScaleMethyl_AATACGCAGG,CGAAGGCATG,AATACGCAGG
ScaleMethyl_AATACGCAGG,TTCAATATAA,AATACGCAGG
ScaleMethyl_AATACGCAGG,CGAATCTCCT,AATACGCAGG
ScaleMethyl_AATACGCAGG,GTCTCCAGAG,AATACGCAGG
ScaleMethyl_AATACGCAGG,TCGAAGCGCG,AATACGCAGG
ScaleMethyl_AATACGCAGG,GCCGTCGCGT,AATACGCAGG
ScaleMethyl_AATACGCAGG,TATTCCGTTA,AATACGCAGG
ScaleMethyl_AATACGCAGG,ACTGGTAGAT,AATACGCAGG
ScaleMethyl_AATACGCAGG,AATTGAGAGA,AATACGCAGG
ScaleMethyl_AATACGCAGG,ATGAGTTCTC,AATACGCAGG
ScaleMethyl_TAACTCTTAG,GGAGGCCTCC,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CAGCAGTATC,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CGAAGGCATG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TTCAATATAA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,CGAATCTCCT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,GTCTCCAGAG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TCGAAGCGCG,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,GCCGTCGCGT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,TATTCCGTTA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,ACTGGTAGAT,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,AATTGAGAGA,TAACTCTTAG
ScaleMethyl_TAACTCTTAG,ATGAGTTCTC,TAACTCTTAG
ScaleMethyl_TCGTAGGCTT,GGAGGCCTCC,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CAGCAGTATC,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CGAAGGCATG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TTCAATATAA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,CGAATCTCCT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,GTCTCCAGAG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TCGAAGCGCG,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,GCCGTCGCGT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,TATTCCGTTA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,ACTGGTAGAT,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,AATTGAGAGA,TCGTAGGCTT
ScaleMethyl_TCGTAGGCTT,ATGAGTTCTC,TCGTAGGCTT""" > ${runDir}/samplesheet.csv

#Set up samples.csv for ScaleMethyl
echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-2H12,ScaleMethyl""" > ${runDir}/samples.csv

```

PM2517_bclconvert.lsf
```bash
#BSUB -W 6:00
#BSUB -q e40medium
#BSUB -n 40
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -o /rsrch4/home/genetics/rmulqueen/projects/metact/240205_RMMM_scalebiotest2/bsub.log
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/projects/metact/240205_RMMM_scalebiotest2
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240205_RMM_scalebiotest2_bclconvert

#load modules
module load nextflow/23.04.3

#set dir variables
proj_dir="/rsrch4/home/genetics/rmulqueen/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/240205_RMMM_scalebiotest2"
fastqDir="${proj_dir}/240205_RMMM_scalebiotest2/PM2517"
samples="${proj_dir}/samples.csv"
genome="${proj_dir}/ref/genome.json"

#run bcl-convert
singularity exec --bind ${runDir} \
~/singularity/bclconvert_3.9.3.sif \
bcl-convert \
--bcl-input-directory ${fastqDir} \
--bcl-num-conversion-threads 10 \
--bcl-num-compression-threads 10 \
--bcl-num-decompression-threads 10 \
--output-directory ${runDir}/fq_out \
--sample-sheet ${runDir}/samplesheet.csv \
--force

#delete these from fq_out directory, otherwise they will throw errors in ScaleMethyl Pipeline (it tries to demux these)
rm -rf ${runDir}/fq_out/Undetermined*fastq.gz #write them out and then check size, then delete to save space

```
Submit the job
```bash
cd /rsrch4/home/genetics/rmulqueen/projects/metact/240205_RMMM_scalebiotest2
rm -rf bsub.log
bsub < PM2517_bclconvert.lsf

"""
Successfully completed.

Resource usage summary:

    CPU time :                                   71442.00 sec.
    Max Memory :                                 72 GB
    Average Memory :                             66.24 GB
    Total Requested Memory :                     100.00 GB
    Delta Memory :                               28.00 GB
    Max Swap :                                   -
    Max Processes :                              6
    Max Threads :                                1109
    Run time :                                   5344 sec.
    Turnaround time :                            5345 sec.
"""
```

## Now to submit a ScaleMethyl job to a seadragon2 node.

Consider max mem per queue
http://hpcweb/using_the_cluster.dir/lsf.dir/MaxMemoryChart.pdf
and max cores and wall time
http://hpcweb/using_the_cluster.dir/lsf.dir/QueueParameters_2023-12-08.pdf

Generate the 240205_scalemet.lsf file for bsub submission.

```bash
#BSUB -W 240:00
#BSUB -q e40long
#BSUB -n 40
#BSUB -M 300
#BSUB -R rusage[mem=300]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bsub.log
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240205_RMM_scalebiotest2

#load modules
module load nextflow/23.04.3
module load singularity/3.7.0
module load glib/2.51.5 #might not need this

#set up directories and variables
export proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export scalebio_nf="${proj_dir}/tools/ScaleMethyl"
export runDir="${proj_dir}/240205_RMMM_scalebiotest2"
export fastqDir="${proj_dir}/240205_RMMM_scalebiotest2/PM2517"
export samples="${proj_dir}/samples.csv"
export genome="${proj_dir}/ref/genome.json"

#set global cache and library directories
export NXF_SINGULARITY_LIBRARYDIR="/rsrch4/home/genetics/rmulqueen/singularity" #might not need these
export NXF_SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity" #might not need these
export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"

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
--max_memory 300.GB \
--max_cpus 40 \
--profile docker,singularity \
-w ${SCRATCH}/work \
-resume

```

## Submit job.

```bash
cd /rsrch4/home/genetics/rmulqueen/projects/metact/240205_RMMM_scalebiotest2
rm -rf bsub.log
bsub < 240205_scalemet.lsf
```

<!-- 

Successfully completed.

Resource usage summary:

    CPU time :                                   6479521.00 sec.
    Max Memory :                                 133 GB
    Average Memory :                             77.06 GB
    Total Requested Memory :                     300.00 GB
    Delta Memory :                               167.00 GB
    Max Swap :                                   -
    Max Processes :                              616
    Max Threads :                                2344
    Run time :                                   121095 sec.
    Turnaround time :                            123590 sec.

The output (if any) is above this job summary.

 -->



## Custom Nextflow pipeline for downstream analysis
240205_scalemet_postprocessing.lsf
```bash
#BSUB -W 240:00
#BSUB -q e40long
#BSUB -n 40
#BSUB -M 300
#BSUB -R rusage[mem=300]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bsub_postprocessing.log
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240205_RMM_scalebiotest2

#load modules
module load nextflow/23.04.3

#set up environment variables 
export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export srcDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src"
export refDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref"
export sif="${srcDir}/copykit.sif"

#call nextflow
nextflow ${srcDir}/singlecell_met_nf.groovy \
--refdir $refDir \
-with-singularity $sif \
-w ${SCRATCH}/met_work \
--scalemethylout ${projDir}/240205_RMMM_scalebiotest2 \
-resume

```

```bash
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2
rm -rf bsub_postprocessing.log
bsub < /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/240205_scalemet_postprocessing.lsf
cat /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bsub_postprocessing.log

cat /rsrch4/scratch/genetics/rmulqueen/met_work/34/3d1b9f1c92717189742e435a09da6d/.command.log 

```









 ## Move to directory with larger storage

 ```bash
mkdir -p /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact
cp -r /rsrch4/home/genetics/rmulqueen/projects/metact /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact
```

## Get interactive node to explore data

```bash
#BSUB -W 24:00
#BSUB -q e80medium
#BSUB -n 80
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup/bsub.out
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240216_scalebio_scbam_split

#load modules
module load samtools
eval "$(/risapps/rhel8/miniconda3/py39_4.12.0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1

#set up functions
#count reads export
count_reads() { 
        samtools view $1 | awk -v b=$1 '{split($1,a,":"); print a[1],b}' | sort | uniq -c | sort -k1,1n
}

#split bams export
split_bams() { 
        test=$1
        idx=$(echo $test | cut -d ' ' -f 2 )
        outidx=$(echo $idx | sed -e 's/+/_/g' -)
        bam=$(echo $test | cut -d ' ' -f 3)
        outprefix=$(echo $bam | cut -d '/' -f 2)
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split($1,a,":"); if(a[1]==i); print $0}')) | samtools view -bS > ./sc_bams/${outprefix}.${idx}.bam
}

export -f count_reads
export -f split_bams

#parallelize it
parallel -j 60 count_reads ::: $(find ./ -maxdepth 2 -name '*bam' -type l ) | sort -k1,1n >unique_read_counts.tsv

#filter to bam files with >100000 unique reads
awk '$1>100000 {print $0}' unique_read_counts.tsv > cells_pf.txt
mkdir -p sc_bams
parallel -j 60 -a cells_pf.txt split_bams

#do CNV calling on split


```

```bash
#BSUB -W 24:00
#BSUB -q e80medium
#BSUB -n 80
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup/bsub.out
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 240216_scalebio_scbam_split

export SINGULARITY_TMPDIR=/singularity_tmpdir



```

```bash
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup
bsub < 240216_scalebio_scbam_split.lsf

#Resource usage summary:

#    CPU time :                                   971259.06 sec.
#    Max Memory :                                 1 GB
#    Average Memory :                             0.96 GB
#    Total Requested Memory :                     100.00 GB
#    Delta Memory :                               99.00 GB
#    Max Swap :                                   -
#    Max Processes :                              365
#    Max Threads :                                366
#    Run time :                                   13759 sec.
#    Turnaround time :                            18382 sec.

```



```bash
bsub -Is -W 6:00 -q transfer -n 4 -M 10 -R rusage[mem=10] /bin/bash

du -shL /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cg_sort_cov/
#13G cg_sort_cov/

rsync -aLPvz \
/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cg_sort_cov \
mulqueen@qcprpgeo.mdanderson.edu:/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2
```



```bash
bsub -Is -W 6:00 -q interactive -n 10 -M 200 -R rusage[mem=200] /bin/bash
singularity shell --bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact /rsrch4/home/genetics/rmulqueen/singularity/scmetR.sif

cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/bamDeDup/sc_bams
```
```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
library(scquantum)
register(MulticoreParam(progressbar = T, workers = 5), default = T)
BiocParallel::bpparam()

dat <- runVarbin(".",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
## Making a sif file (see sif_creation.md) to do further copykit processing.

```

```R


setwd("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/matrix")
Matrix::readMM('.')

```