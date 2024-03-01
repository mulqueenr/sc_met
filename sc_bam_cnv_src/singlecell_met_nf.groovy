// Declare syntax version
//https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf
//Runs entirely in singularity container

nextflow.enable.dsl=2
// Script parameters

// Input parameters, user specified defaults
params.scalemethylout = "/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2"
params.readcountfilter = 100000

// Reference files
params.projectdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
params.refdir="${params.projectdir}/ref"
params.genes_bed="${params.refdir}/grch38/filteredGTF/GRCh38_transcriptsOnly.gtf"
params.genome_length="${params.refdir}/grch38/genome.txt"
params.src_dir="${params.projectdir}/src/"
params.out_dir="${params.scalemethylout}"
log.info """

		================================================
		    SCALEBIO CONTINUED PROCESSING v0.1
		================================================
		ScaleMethyl Pipeline Output : ${params.scalemethylout}
		Cell Readcount Filter : ${params.readcountfilter}
		NF Working Dir : ${workflow.launchDir}
		Src directory : ${params.src_dir}
		Out directory : ${params.out_dir}
		================================================

""".stripIndent()

// CNV BLOCK //
process COUNT_PER_GROUPEDBAM { 
	// Generate a count per grouped bam and pass list.
	input:
		path bams
	output:
		path("*.readcount.tsv")
	script:
	"""
	bam_name="${bams}"
	cellline=\$(echo \$bam_name | cut -d \'.\' -f 1 )
	well=\$(echo \$bam_name | cut -d \'.\' -f 2 )
	samtools view ${bams} \\
	| awk -v b=\$1 \'{split(\$1,a,\":\"); print a[1],b}\' \\
	| sort \\
	| uniq -c \\
	| sort -k1,1n \\
	| awk \'\$1>${params.readcountfilter} {print \$0}\' > \${cellline}.\${well}.readcount.tsv
        """
}

process SPLIT_GROUPED_BAM { 
	//Generate single-cell bams from those that pass filter
	input:
		path readcount
	output:
		path("*sorted.bam")
	script:
	"""
	test=${readcount}
	idx=\$(echo \$test | cut -d \' \' -f 2 )
	outidx=\$(echo \$idx | sed -e \'s/+/_/g\' -)
	bam=\$(echo \$test | cut -d \' \' -f 3)
	outprefix=\$(echo \$bam | cut -d \'/\' -f 2)

	((samtools view -H \$bam) && (samtools view \$bam | awk -v i=\$idx \'{split(\$1,a,\":\"); if(a[1]==i); print \$0}\')) \\
	| samtools view -bS - \\
	| samtools sort -T . -O BAM -o \${outprefix}.\${outidx}.sorted.bam -

	"""
}

process RUN_COPYKIT {
	//Run Copykit on single-cell bams
	publishDir "${params.scalemethylout}/postprocessing", mode: 'copy', overwrite: true

	input:
		path sc_bam_path
	output:
		path("*.{rds,pdf}")
	script:
	"""
	Rscript ${params.projectdir}/src/copykit_run.R \\
	-i "."
	"""
}


//process RUN_RIDDLER

// CELL CLUSTERING AND TYPING BLOCK //
process MAKE_TRANSCRIPT_BED {
	//Generate transcript bed from ScaleMethyl GTF file, filter to unique transcripts and the longest one.
	//requires bedtools and parallel
	cpus 10
	input:
		path ref
	output:
		path("GRCh38_transcripts.longest.bed")
	script:
	"""	
	awk \'OFS=\"\\t\" {split(\$10,a,\"\\\"\"); split(\$14,b,\"\\\"\"); print \"chr\"\$1,\$4,\$5,a[2],b[2]}\' ${ref} > GRCh38_transcripts.bed

	select_longest() { 
		awk -v gene=\"\$1\" \'{if(\$5==gene){print \$0,\$3-\$2}}\' GRCh38_transcripts.bed \\
		| sort -k6,6n - \\
		| tail -n 1
	}

	export -f select_longest

	parallel -j ${task.cpus} select_longest ::: \$(awk \'{print \$5}\' GRCh38_transcripts.bed | uniq) > GRCh38_transcripts.longest.bed

	"""
}

process MAKE_100KB_BED {
	//Generate transcript bed from ScaleMethyl genome.txt file
	//Add Blacklist filter??
	//requires bedtools
	input:
		path ref
	output:
		path("genome_windows.100kb.bed")

	script:
	"""
	grep -v \"^K\" ${ref} \\
	| grep -v \"^G\" \\
	| awk \'OFS=\"\\t\" {print \"chr\"\$1,\$2}\' > genome.filt.txt
	bedtools makewindows -w 100000 -g genome.filt.txt | awk \'OFS=\"\\t\" {print \$1,\$2,\$3,\$1\"_\"\$2\"_\"\$3}\' > genome_windows.100kb.bed
	"""
}	


//process MAKE_TF_MOTIF_BED {}
process SUMMARIZE_CG_OVER_BEDFEATURES {
	publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true, pattern: "*.tsv.gz"
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	//Output data frames of 1. coverage 2. methylation rate 3. posterior estimate (epiclomal)
	//requires from pybedtools import BedTool
	//requires import multiprocessing
	//requires import pandas as pd
	//requires import numpy as np
	//requires import os, glob
	//requires import scipy.io, scipy.sparse
	//requires import numpy as np
	//requires import argparse, pathlib
	input:
		path cov_folder
		path feat
	output:
		path("*tsv.gz")
	script:
	"""
	python /src/mc_cov_feature_summary.py --features ${feat} --cov_folder ${cov_folder}
	"""
}

//process MAKE_METHYLATION_SEURAT_OBJ {
	//Make two assays, 100kb for clustering, and gene met for identification
	//Perform dim reduction
	//Perform label transfer with RNA object????
//}

workflow {
	/* RUN CNV PROFILING ON CNVS */	
		def bams = Channel.fromPath("${params.scalemethylout}/bamDeDup/*/*bam")
		bams.view()

		scbam_dir= COUNT_PER_GROUPEDBAM(bams) \
		| SPLIT_GROUPED_BAM \
		| collect

		RUN_COPYKIT(scbam_dir)
		//RIDDLER_ON_SPLIT_BAMS(scbam_dir)


	/* GENERATE CLUSTERS AND GENE SUMMARY WINDOWS FOR CELL TYPE ANALYSIS */
		def covs = Channel.fromPath("${params.scalemethylout}/cg_sort_cov/**/**")

		gene_bed = MAKE_TRANSCRIPT_BED(${params.genes_bed})
		100kb_bed = MAKE_100KB_BED(${params.genome_length})
		//tf_bed = MAKE_TF_BED


		//#grab bed.gz, subset window bed to just that chr, make window x cell matrix, merge window x cell matrix column wise (add other cells), merge window x cell matrix row wise (add subset window beds), save it as a mtx file format

}


/*
bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash



#load modules
module load nextflow/23.04.3
module load singularity

#set up environment variables 
export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export srcDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src"
export sif="${srcDir}/copykit.sif"

#set up nextflow variables 
#export APPTAINER_BIND="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
#export APPTAINER_BINDPATH="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2

nextflow ${srcDir}/singlecell_met_nf.groovy \
--refdir $refdir \
-with-singularity $sif \
-w ${SCRATCH}/met_work \
--scalemethylout ${projDir}/240205_RMMM_scalebiotest2 \
-resume



cd /rsrch4/scratch/genetics/rmulqueen/met_work/72/390e90fff9625c529aa3d993044841
singularity shell --bind /rsrch4/scratch/genetics/rmulqueen/met_work --bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2 $sif
*/



/*



-make riddler container
-make epiclomal container
-make met seurat container


https://github.com/yardimcilab/RIDDLER for CNV calls??


https://github.com/molonc/Epiclomal
-cluster on large windows
-generate gene body met signals to highlight clusters

https://www.nature.com/articles/s41586-020-03182-8#Sec12

Methylation features were calculated as fractions of methylcytosine over total cytosine across gene bodies ± 2kb flanking regions and 100kb bins spanning the entire genome. Methylation features were further split into CG and CH methylation types. Features overlapping our methylation mm10 blacklist were removed. 100kb bin features were then filtered on minimum mean coverage >500 and maximum mean coverage <3000. Gene body features were filtered on minimum coverage >5 and all remaining features were normalized per cell using the beta binomial normalization technique in allcools.16

https://www.sciencedirect.com/science/article/pii/S2666979X23002987#sec5

*/
