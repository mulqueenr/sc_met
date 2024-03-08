nextflow.enable.dsl=2
// Script parameters
// Input parameters, user specified defaults
params.scalemethylout = "/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2"
params.readcountfilter = 100000

// Reference files
params.projectdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
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
		NF Working Dir : ${workflow.workDir}
		Src directory : ${params.src_dir}
		Out directory : ${params.out_dir}
		================================================

""".stripIndent()

// CNV BLOCK //
process SPLIT_GROUPED_BAM { 
	// Generate a count per grouped bam and pass list.
	maxForks 10
	cpus 5
	label 'cnv'

	input:
		path bams
	output:
		path("*.sorted.bam"), optional: true
	script:
	"""
	bam_name="${bams}"
	export cellline=\$(echo \$bam_name | cut -d \'.\' -f 1 )
	export well=\$(echo \$bam_name | cut -d \'.\' -f 2 )
	samtools view -@ ${task.cpus} ${bams} \\
	| awk -v b=${bams} \'{split(\$1,a,\":\"); print a[1],b}\' \\
	| sort \\
	| uniq -c \\
	| sort -k1,1n \\
	| awk \'\$1>${params.readcountfilter} {print \$0}\' > \${cellline}.\${well}.readcount.tsv


	split_bam() {
	test=\$1
	idx=\$(echo \$test | cut -d \' \' -f 2 )
	outidx=\$(echo \$idx | sed -e \'s/+/_/g\' -)
	bam=\$(echo \$test | cut -d \' \' -f 3)

	((samtools view -H \$bam) && (samtools view \$bam | awk -v i=\$idx \'{split(\$1,a,\":\"); if(a[1]==i) print \$0}\')) \\
	| samtools view -bS - \\
	| samtools sort -T . -O BAM -o \${cellline}.\${well}.\${outidx}.sorted.bam - 
	}
    
    export -f split_bam
	parallel -j ${task.cpus} -a \${cellline}.\${well}.readcount.tsv split_bam 
    """
}

process RUN_COPYKIT {
	//Run Copykit on single-cell bams
	publishDir "${params.scalemethylout}/postprocessing", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	cpus 10
	label 'cnv'

	input:
		path sc_bam_path
	output:
		path("*.{rds,pdf}")
	script:
	"""
	Rscript /src/copykit_run.R \\
	-i "." \\
	-c ${task.cpus}
	"""
}

// CELL CLUSTERING AND TYPING BLOCK //
process MAKE_TRANSCRIPT_BED {
	//Generate transcript bed from ScaleMethyl GTF file, filter to unique transcripts and the longest one.
	//requires bedtools and parallel
	cpus 10
	label 'cnv'

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

	parallel -j ${task.cpus} select_longest ::: \$(awk \'OFS =\"\\t\" {print \$5}\' GRCh38_transcripts.bed | uniq) |  awk \'OFS=\"\\t\" {print \$1,\$2,\$3,\$4}\' > GRCh38_transcripts.longest.bed

	"""
}

process MAKE_100KB_BED {
	//Generate transcript bed from ScaleMethyl genome.txt file
	//Add Blacklist filter??
	//requires bedtools
	label 'cnv'

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
process SUMMARIZE_CG_OVER_BED_100KB {
	//publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true, pattern: "*.tsv.gz"
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'

	input:
		path cov_folder
		path hundokb_bed
	output:
		path("*.tsv.gz")
	script:
	"""
	python /src/mc_cov_feature_summary.py \\
	--features ${hundokb_bed} \\
	--feat_name 100kb \\
	--cov_folder ${cov_folder} # --cpus 5
	"""
}


//process MAKE_TF_MOTIF_BED {}
process SUMMARIZE_CG_OVER_BED_GENES {
	//publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'

	input:
		path cov_folder
		path gene_bed
	output:
		path("*.tsv.gz")
	script:
	"""
	python /src/mc_cov_feature_summary.py \\
	--features ${gene_bed} \\
	--feat_name genebody \\
	--cov_folder ${cov_folder} #--cpus 5
	"""
}

process MERGED_100KB_SUMMARIES {
	publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'

	input:
		path cg_summaries
	output:
		path("*.tsv.gz")
	script:
	"""
	#concatenate total count
	set -- *100kb.total_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > total_count.100kb.merged.tsv.gz

	#concatentate mc count
	set -- *100kb.mc_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_count.100kb.merged.tsv.gz

	#concatenate rate estimates
	set -- *100kb.mc_posteriorest.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_posteriorest.100kb.merged.tsv.gz
	"""
}	


process MERGED_GENE_SUMMARIES {
	publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'

	input:
		path cg_summaries
	output:
		tuple path("total_count.genebody.merged.tsv.gz"), path("mc_count.genebody.merged.tsv.gz"), path("mc_posteriorest.genebody.merged.tsv.gz")
	script:
	"""
	#concatenate total count
	set -- *genebody.total_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > total_count.genebody.merged.tsv.gz

	#concatentate mc count
	set -- *genebody.mc_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_count.genebody.merged.tsv.gz

	#concatenate rate estimates
	set -- *genebody.mc_posteriorest.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_posteriorest.genebody.merged.tsv.gz
	"""
}	


///USE BC_MULTIOME.SIF TEMPORARILY FOR THIS
process MAKE_FINAL_SEURATOBJ {
	publishDir "${params.scalemethylout}/seurat_out", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype_seurat'

	input:
		tuple path(hundokb_total),path(hundokb_mc),path(hundokb_rate)
		tuple path(genebody_total),path(genebody_mc),path(genebody_rate)
		path(meta)
	output:
		path("met.SeuratObject.rds")
	script:
	"""
	Rscript /src/seurat_setup.R
	"""
}	
	//Perform dim reduction
	//Perform label transfer with RNA object????
//}

workflow {
	/* RUN CNV PROFILING ON CNVS */	
		def bams = 
		Channel.fromPath("${params.scalemethylout}/bamDeDup/*/*bam")

		scbam_dir = 
		SPLIT_GROUPED_BAM(bams) \
		| collect

		RUN_COPYKIT(scbam_dir)


	/* GENERATE CLUSTERS AND GENE SUMMARY WINDOWS FOR CELL TYPE ANALYSIS */
		def covs = 
		Channel.fromPath("${params.scalemethylout}/cg_sort_cov/*/*chroms.sort/", type: 'dir')

		def meta_data = 
		Channel.fromPath("${params.scalemethylout}/report/*/csv/*.passingCellsMapMethylStats.csv")
		met_in = meta_data | collect

		hundokb_bed = MAKE_100KB_BED("${params.genome_length}")
		hundokb_out =
		SUMMARIZE_CG_OVER_BED_100KB(covs, hundokb_bed) \
		| collect \
		| MERGED_100KB_SUMMARIES

		gene_bed = MAKE_TRANSCRIPT_BED("${params.genes_bed}")
		genebody_out = 
		SUMMARIZE_CG_OVER_BED_GENES(covs, gene_bed) \
		| collect \
		| MERGED_GENE_SUMMARIES

		MAKE_FINAL_SEURATOBJ(hundokb_out,genebody_out,met_in)

		//tf_bed = MAKE_TF_BED

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
export refDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref"
export sif="${srcDir}/copykit.sif"

#call nextflow
nextflow ${srcDir}/singlecell_met_nf.groovy \
--refdir $refDir \
-with-singularity $sif \
-w ${SCRATCH}/met_work \
--scalemethylout ${projDir}/240205_RMMM_scalebiotest2 \
-resume

cd /rsrch4/scratch/genetics/rmulqueen/met_work/c4/2ef5236657a57879d9e2c811c4b7aa
export sif="${srcDir}/scmetR.sif"

singularity shell \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src:/src/ \
--bind /rsrch4/scratch/genetics/rmulqueen/met_work \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2 \
--bind /rsrch4/scratch/genetics/rmulqueen/work \
$sif



*/



/*


ADD OVERDISPERSION
MAD SCORE
BREADTH OF COVERAGE

-make epiclomal container




https://github.com/molonc/Epiclomal
-cluster on large windows
-generate gene body met signals to highlight clusters

https://www.nature.com/articles/s41586-020-03182-8#Sec12

Methylation features were calculated as fractions of methylcytosine over total cytosine across gene bodies ± 2kb flanking regions and 100kb bins spanning the entire genome. Methylation features were further split into CG and CH methylation types. Features overlapping our methylation mm10 blacklist were removed. 100kb bin features were then filtered on minimum mean coverage >500 and maximum mean coverage <3000. Gene body features were filtered on minimum coverage >5 and all remaining features were normalized per cell using the beta binomial normalization technique in allcools.16

https://www.sciencedirect.com/science/article/pii/S2666979X23002987#sec5

*/
