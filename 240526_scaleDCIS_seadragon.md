## Setting up the ScaleMethyl Nextflow pipeline to work on a IBM LSF HPC

Setting up environment following https://github.com/ScaleBio/ScaleMethyl
On mdandersons seadragon.

Transfer Our Data from our labs node 
```bash
#get data from cell line run onto seadragon
proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/250526_RMMM_scalebio_dcis"
mkdir -p ${proj_dir}
cd ${proj_dir}
sftp mulqueen@qcprpgeo.mdanderson.edu
get -r /volumes/seq/flowcells/MDA/NovaSeq/2024/240523_VH00219_594_AAFLYGNM5 #grab full rundir
```

## Running on seadragon
See 240115_scalealphatest2_seadragon2 for more details about singularity setup
Initial set up of files and running bcl-convert.
Note: I've been running this first in an interactive node to confirm the indexes off the first bit of reads. Then putting it onto a larger node with more cores.

```bash
bsub -Is -W 24:00 -q medium -n 10 -M 100 -R rusage[mem=100] /bin/bash #small interactive node for bcl-convert 

module load nextflow/23.04.3

#set dir variables
proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
scalebio_nf="${proj_dir}/tools/ScaleMethyl"
runDir="${proj_dir}/250526_RMMM_scalebio_dcis"
fastqDir="${runDir}/240523_VH00219_594_AAFLYGNM5"
samples="${runDir}/samples.csv"
genome="${proj_dir}/ref/genome.json"

cd $runDir

#used i5 index set 2 and i7 index set 1
#ran first and found all the reads assigned to undetermined, then looked up resulting indexes to confirm indexes we expect
#use i7.txt and i5.txt for making sample sheet
#GTCCTACTTG  A02
#AGGAGAGAAC  B02
#CATAGTTCGG  C02
#CCATAAGTCC  D02
#TGGCCGGCCT  E02
#AGCTGCAATA  F02
#TCATCTCGAT  G02
#GATGATCGTA  H02

#NOTE: REVERSE COMPLEMENT OF I2 because its on the nextseq, for novaseq I2 is consistent with cell sample names
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,1
TrimUMI,0
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,2
MinimumTrimmedReadLength,16
MaskShortReads,16
OverrideCycles,Y151;I10;I10;U8Y143
[BCLConvert_Data]
Sample_ID,index,index2
ScaleMethyl_GTCCTACTTG,GGAGGCCTCC,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,CAGCAGTATC,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,CGAAGGCATG,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,TTCAATATAA,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,CGAATCTCCT,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,GTCTCCAGAG,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,TCGAAGCGCG,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,GCCGTCGCGT,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,TATTCCGTTA,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,ACTGGTAGAT,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,AATTGAGAGA,CAAGTAGGAC
ScaleMethyl_GTCCTACTTG,ATGAGTTCTC,CAAGTAGGAC
ScaleMethyl_AGGAGAGAAC,GGAGGCCTCC,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,CAGCAGTATC,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,CGAAGGCATG,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,TTCAATATAA,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,CGAATCTCCT,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,GTCTCCAGAG,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,TCGAAGCGCG,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,GCCGTCGCGT,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,TATTCCGTTA,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,ACTGGTAGAT,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,AATTGAGAGA,GTTCTCTCCT
ScaleMethyl_AGGAGAGAAC,ATGAGTTCTC,GTTCTCTCCT
ScaleMethyl_CATAGTTCGG,GGAGGCCTCC,CCGAACTATG
ScaleMethyl_CATAGTTCGG,CAGCAGTATC,CCGAACTATG
ScaleMethyl_CATAGTTCGG,CGAAGGCATG,CCGAACTATG
ScaleMethyl_CATAGTTCGG,TTCAATATAA,CCGAACTATG
ScaleMethyl_CATAGTTCGG,CGAATCTCCT,CCGAACTATG
ScaleMethyl_CATAGTTCGG,GTCTCCAGAG,CCGAACTATG
ScaleMethyl_CATAGTTCGG,TCGAAGCGCG,CCGAACTATG
ScaleMethyl_CATAGTTCGG,GCCGTCGCGT,CCGAACTATG
ScaleMethyl_CATAGTTCGG,TATTCCGTTA,CCGAACTATG
ScaleMethyl_CATAGTTCGG,ACTGGTAGAT,CCGAACTATG
ScaleMethyl_CATAGTTCGG,AATTGAGAGA,CCGAACTATG
ScaleMethyl_CATAGTTCGG,ATGAGTTCTC,CCGAACTATG
ScaleMethyl_CCATAAGTCC,GGAGGCCTCC,GGACTTATGG
ScaleMethyl_CCATAAGTCC,CAGCAGTATC,GGACTTATGG
ScaleMethyl_CCATAAGTCC,CGAAGGCATG,GGACTTATGG
ScaleMethyl_CCATAAGTCC,TTCAATATAA,GGACTTATGG
ScaleMethyl_CCATAAGTCC,CGAATCTCCT,GGACTTATGG
ScaleMethyl_CCATAAGTCC,GTCTCCAGAG,GGACTTATGG
ScaleMethyl_CCATAAGTCC,TCGAAGCGCG,GGACTTATGG
ScaleMethyl_CCATAAGTCC,GCCGTCGCGT,GGACTTATGG
ScaleMethyl_CCATAAGTCC,TATTCCGTTA,GGACTTATGG
ScaleMethyl_CCATAAGTCC,ACTGGTAGAT,GGACTTATGG
ScaleMethyl_CCATAAGTCC,AATTGAGAGA,GGACTTATGG
ScaleMethyl_CCATAAGTCC,ATGAGTTCTC,GGACTTATGG
ScaleMethyl_TGGCCGGCCT,GGAGGCCTCC,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,CAGCAGTATC,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,CGAAGGCATG,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,TTCAATATAA,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,CGAATCTCCT,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,GTCTCCAGAG,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,TCGAAGCGCG,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,GCCGTCGCGT,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,TATTCCGTTA,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,ACTGGTAGAT,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,AATTGAGAGA,AGGCCGGCCA
ScaleMethyl_TGGCCGGCCT,ATGAGTTCTC,AGGCCGGCCA
ScaleMethyl_AGCTGCAATA,GGAGGCCTCC,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,CAGCAGTATC,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,CGAAGGCATG,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,TTCAATATAA,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,CGAATCTCCT,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,GTCTCCAGAG,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,TCGAAGCGCG,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,GCCGTCGCGT,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,TATTCCGTTA,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,ACTGGTAGAT,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,AATTGAGAGA,TATTGCAGCT
ScaleMethyl_AGCTGCAATA,ATGAGTTCTC,TATTGCAGCT
ScaleMethyl_TCATCTCGAT,GGAGGCCTCC,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,CAGCAGTATC,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,CGAAGGCATG,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,TTCAATATAA,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,CGAATCTCCT,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,GTCTCCAGAG,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,TCGAAGCGCG,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,GCCGTCGCGT,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,TATTCCGTTA,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,ACTGGTAGAT,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,AATTGAGAGA,ATCGAGATGA
ScaleMethyl_TCATCTCGAT,ATGAGTTCTC,ATCGAGATGA
ScaleMethyl_GATGATCGTA,GGAGGCCTCC,TACGATCATC
ScaleMethyl_GATGATCGTA,CAGCAGTATC,TACGATCATC
ScaleMethyl_GATGATCGTA,CGAAGGCATG,TACGATCATC
ScaleMethyl_GATGATCGTA,TTCAATATAA,TACGATCATC
ScaleMethyl_GATGATCGTA,CGAATCTCCT,TACGATCATC
ScaleMethyl_GATGATCGTA,GTCTCCAGAG,TACGATCATC
ScaleMethyl_GATGATCGTA,TCGAAGCGCG,TACGATCATC
ScaleMethyl_GATGATCGTA,GCCGTCGCGT,TACGATCATC
ScaleMethyl_GATGATCGTA,TATTCCGTTA,TACGATCATC
ScaleMethyl_GATGATCGTA,ACTGGTAGAT,TACGATCATC
ScaleMethyl_GATGATCGTA,AATTGAGAGA,TACGATCATC
ScaleMethyl_GATGATCGTA,ATGAGTTCTC,TACGATCATC""" > ${runDir}/samplesheet.csv


#spot check indexes
#zcat Undetermined_S0_L001_I1_001.fastq.gz | sed -n '2~4p' - | head #
#zcat Undetermined_S0_L001_I2_001.fastq.gz | sed -n '2~4p' - | head #i5 looks right (matches with https://github.com/ScaleBio/ScaleMethyl/blob/main/references/i5.txt)

#Set up samples.csv for ScaleMethyl
echo """sample,barcodes,libName
DCIS_41T,3A01-3H08,ScaleMethyl
DCIS_66T,3A09-3H11,ScaleMethyl
TNBC,3A12-3H12,ScaleMethyl""" > ${runDir}/samples.csv

```


240523_VH00219_594_AAFLYGNM5_bclconvert.lsf
```bash
#BSUB -W 6:00
#BSUB -q e80medium
#BSUB -n 40
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/250526_RMMM_scalebio_dcis/bsub.log
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/250526_RMMM_scalebio_dcis
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 250526_RMMM_scalebio_dcis_bclconvert

#set dir variables
export proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export runDir="${proj_dir}/250526_RMMM_scalebio_dcis"
export fastqDir="${runDir}/240523_VH00219_594_AAFLYGNM5"

#run bcl-convert
singularity exec --bind ${runDir} \
${proj_dir}/src/bclconvert_3.9.3.sif \
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

## Now to submit a ScaleMethyl job to a seadragon2 node.

Consider max mem per queue
http://hpcweb/using_the_cluster.dir/lsf.dir/MaxMemoryChart.pdf
and max cores and wall time
http://hpcweb/using_the_cluster.dir/lsf.dir/QueueParameters_2023-12-08.pdf

Generate the lsf file for bsub submission.

```bash
#BSUB -W 240:00
#BSUB -q e80long
#BSUB -n 80
#BSUB -M 300
#BSUB -R rusage[mem=300]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/250526_RMMM_scalebio_dcis/%J_scalemet_processing.log
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/250526_RMMM_scalebio_dcis
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J 250526_RMMM_scalebio_dcis

#load modules
module load nextflow/23.04.3
module load singularity/3.7.0
module load glib/2.51.5 #might not need this

#set up directories and variables
export proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export scalebio_nf="${proj_dir}/tools/ScaleMethyl"
export runDir="${proj_dir}/250526_RMMM_scalebio_dcis"
export genome="${proj_dir}/ref/genome.json"
export fastqDir="${runDir}/240523_VH00219_594_AAFLYGNM5"
export samples="${runDir}/samples.csv"

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
--bamDedupOut true \
--covOut true \
--max_memory 300.GB \
--max_cpus 80 \
--profile docker,singularity \
-w ${SCRATCH}/scalemet_work \
-resume

```

```bash
bsub < scalemet.lsf
```

## Submit job.

```bash
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240526_RMMM_scalebio_dcis
rm -rf bsub.log
bsub < scalemet.lsf

mkdir -p transfer_dat
find ./cg_sort_cov -type l -exec bash -c 'cp -R "$(readlink -m "$0")" ./transfer_dat' {} \; #scale pipeline makes empty files, this throws errors for empty files (but can be ignored)
cp -R ./report ./transfer_dat
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240526_RMMM_scalebio_dcis/cg_sort_cov
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

premethyst for DCIS samples (run on geo)
```bash
singularity shell \
--bind $HOME \
--bind /volumes/seq/projects/metACT/ \
~/singularity/amethyst.sif

export cg_sort="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/cg_sort_cov"
export task_cpus=50
mkdir -p ${cg_sort}/h5_files
cd $cg_sort
find $cg_sort -maxdepth 1 -type d -printf '%f\n' | grep "sort$" > cg_cov_folders.txt

premethyst() {
outname=$(echo $1 | awk '{ gsub(".CG.chroms.sort", ""); print }')
indir=$1
echo $indir" "$outname
python ~/src/premethyst_calls2h5.py $cg_sort/$indir ${cg_sort}/h5_files/${outname}
}

export -f premethyst
parallel -j ${task_cpus} -a cg_cov_folders.txt premethyst


```