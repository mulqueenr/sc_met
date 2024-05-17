nextflow.enable.dsl=2
// Script parameters
// Input parameters, user specified defaults
params.scalemethylout = "/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2"
params.readcountfilter = 100000

// Reference files
params.projectdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
params.genes_bed="${params.refdir}/grch38/filteredGTF/GRCh38_transcriptsOnly.gtf"
params.metatlas_bed="${params.refdir}/met_atlas.hg38.feat.bed"
params.multiome_atac_bed="${params.refdir}/multiome_bc.marker_atac_peaks.bed"
params.genome_length="${params.refdir}/grch38/genome.txt"
params.bismark_ref="${params.refdir}/grch38/"
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
		Genes Bed : ${params.genes_bed}
		Met Atlas Bed : ${params.metatlas_bed}
		Multiome ATAC Bed : ${params.multiome_atac_bed}
		Feature to summarize over : ${params.feat}
		CNV Calling : ${params.cnv_call}
		================================================

""".stripIndent()

// CNV BLOCK //
process SPLIT_GROUPED_BAM { 
	publishDir "${params.scalemethylout}/postprocessing/sc_bams", mode: 'copy', overwrite: true
	// Generate a count per grouped bam and pass list.
	maxForks 10
	cpus 5
	label 'cnv'

	input:
		path bams
	output:
		tuple path("*.readcount.tsv"),path("*.sorted.bam"), optional: true
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

process SPLIT_GROUPED_FASTQ { 
	publishDir "${params.scalemethylout}/postprocessing/sc_fq", mode: 'copy', overwrite: true
	// Generate a count per grouped fastq and pass list to filter down to single cells.
	// Note this requires both R1 and R2 to pass thresholds for each cell, which shouldn't be a problem really.
	maxForks 40
	cpus 2
	label 'cnv'

	input:
		tuple val(sampleid),path(fq1),path(fq2)
	output:
		tuple path("*-R1.fq.gz"),path("*-R2.fq.gz"), optional: true

	script:
	"""
	fq_name="${sampleid}"
	export cellline=\$(echo \$fq_name | cut -d \'.\' -f 1 )
	export well=\$(echo \$fq_name | cut -d \'.\' -f 2 | cut -d \'_\' -f 1)

	#count fq1 and use to filter possible cell ids
	zcat ${fq1} \\
	| grep "^@" \\
	| awk '{a=substr(\$1,2);split(a,b,":");print b[1]}' \\
	| sort \\
	| uniq -c \\
	| sort -k1,1n \\
	| awk \'\$1>${params.readcountfilter} {print \$0}\' > \${cellline}.\${well}.readcount.tsv

	#split fastq by cell ides
	split_fastq() {
	test=\$1
	idx=\$(echo \$test | cut -d \' \' -f 2 )
	outidx=\$(echo \$idx | sed -e \'s/+/_/g\' -)
	zcat ${fq1} | grep -A 3 "^@\${idx}" | grep -v "^\\-\\-\$"| gzip > \${cellline}.\${well}.\${outidx}-R1.fq.gz
	zcat ${fq2} | grep -A 3 "^@\${idx}" | grep -v "^\\-\\-\$"| gzip > \${cellline}.\${well}.\${outidx}-R2.fq.gz
	}
	export -f split_fastq
    
	parallel -j ${task.cpus} -a \${cellline}.\${well}.readcount.tsv split_fastq
    """
}


process BISMARK_REALIGN { 
	publishDir "${params.scalemethylout}/postprocessing/bismark/sc_bam", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.bismark_ref}:/ref/,${workflow.workDir}"
	maxForks 40
	cpus 1
	label 'allcool'

	input:
		tuple path(fq1),path(fq2)
	output:
		path("*deduplicated.bam")
	script:
	"""
	source activate bismark_env

	basename=\$(echo "${fq1.baseName}" | sed 's/.\\{6\\}\$//')
	cutadapt \
	-a CTGCGACGGCTGC \
	-A CTGTCTCTTAT \
	-u 8 \
	-U 8 \
	-o out1.fq.gz \
	-p out2.fq.gz \
	${fq1} \
	${fq2}

	bismark \\
	--pbat \\
	--temp_dir . \\
	--basename \${basename} \\
	--genome /ref/ \\
	-1 out1.fq.gz \\
	-2 out2.fq.gz

	deduplicate_bismark \\
	--paired \\
	--bam \\
	--outfile \${basename} \\
	\${basename}_pe.bam
    """
}

process BISMARK_TO_ALLCOOL {
	containerOptions "--bind ${params.bismark_ref}:/ref/"
	label 'allcool'
	
	input:
		path(bam)
	output:
		path("*{.gz,.gz.tbi}")
	script:
	"""
	source activate allcool_env 

	bam_name="${bam}"
	export cellline=\$(echo \$bam_name | cut -d \'.\' -f 1 )
	export well=\$(echo \$bam_name | cut -d \'.\' -f 2 )
	export idx=\$(echo \$bam_name | cut -d \'.\' -f 3 )
	export out_name=\"\${cellline}.\${well}.\${idx}\"

	samtools sort \\
	-o \${out_name}.sorted.bam \\
	${bam}

	allcools bam-to-allc \\
	--reference_fasta /ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
	--bam_path \${out_name}.sorted.bam  \\
	--output_path \${out_name}
	"""
}


process ALLCOOL_TO_DATASET {
	publishDir "${params.scalemethylout}/postprocessing/bismark/sc_allcool/", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.refdir}:/ref/"
	label 'allcool'
	cpus 50

	input:
		path(allc)
	output:
		path("*")
	script:
	"""
	source activate allcool_env 

	#move each sample into subdir and make a table
	for i in *.gz; do 
	out_name=\${i::-3};
	mkdir -p \${out_name};
	mv \${out_name}.gz \${out_name}.gz.tbi \${out_name}/;
	echo \${out_name}\$'\t'./\${out_name}/\${out_name}.gz; done > allc_table.tsv

	allcools generate-dataset  \\
	--allc_table allc_table.tsv \\
	--output_path data.mcds \\
	--chrom_size_path /ref/allcool_regions/genome.txt \\
	--obs_dim cell  \\
	--cpu ${task.cpus} \\
	--chunk_size 30 \\
	--regions chrom100k 100000 \\
	--quantifiers chrom100k count CGN \\
	--regions chrom5k 5000 \\
	--quantifiers chrom5k count CGN \\
	--regions transcripts /ref/allcool_regions/GRCh38_transcripts.longest.bed.gz \\
	--quantifiers transcripts count CGN \\
	--regions CGI /ref/allcool_regions/grch38.cgi.bed.gz \\
	--quantifiers CGI count CGN \\
	--regions METATLAS /ref/allcool_regions/met_atlas.hg38.feat.bed.gz \\
	--quantifiers METATLAS count CGN \\
	--regions ATAC /ref/allcool_regions/multiome_bc.marker_atac_peaks.bed.gz \\
	--quantifiers ATAC count CGN \\
	--quantifiers CGI hypo-score CGN cutoff=0.9 \\
	--quantifiers ATAC hypo-score CGN cutoff=0.9 \\
	--quantifiers METATLAS hypo-score CGN cutoff=0.9
	
	"""

}

process BAM_TO_BED_MET { 
	// Generate a bed file with read start-end mC count and methylation level
	//Theory here is that cell types have CGI which define differences, from the methylation atlas paper, so it should be at read level
	//Filter by XC:i:0 to ensure read is bs converted (based on CHH methylation)
	//Then count X (met CG) and x (nonmet CG) per read
	publishDir "${params.scalemethylout}/postprocessing/read_bed", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'cnv'

	input:
		path sorted_bams
	output:
		path("*.CGread.bed.gz")
	script:
	"""
	cellline=\$(echo "${sorted_bams}" | cut -d '.' -f 1 )
	well=\$(echo "${sorted_bams}" | cut -d '.' -f 2 )
	outidx=\$(echo "${sorted_bams}" | cut -d '.' -f 3 )

	#generate read bed 
	samtools view ${sorted_bams} \\
	| awk 'OFS=\"\\t\"{split(\$1,a,\":\");if(\$15==\"XC:i:0\")print \$3,\$4,\$4+\$9,a[1]}' \\
	| gzip > \${cellline}.\${well}.\${outidx}.read.bed.gz

	#generate x|X count
	samtools view ${sorted_bams} \\
	| awk 'OFS=\"\\t\"{if(\$15==\"XC:i:0\")print substr(\$14,6)}' \\
	|  sed 's/[^X|x]//g' \\
	| awk '{ print length }' \\
	| gzip > \${cellline}.\${well}.\${outidx}.CG.bed.gz

	#generate X count
	samtools view ${sorted_bams} \\
	| awk 'OFS=\"\\t\"{if(\$15==\"XC:i:0\")print substr(\$14,6)}' \\
	|  sed 's/[^X]//g' \\
	| awk '{ print length }' \\
	| gzip > \${cellline}.\${well}.\${outidx}.mCG.bed.gz

	command=paste
	for i in \\
	\${cellline}.\${well}.\${outidx}.read.bed.gz \\
	\${cellline}.\${well}.\${outidx}.CG.bed.gz \\
	\${cellline}.\${well}.\${outidx}.mCG.bed.gz; do
	    command=\"\$command <(gzip -cd \$i)\"
	done

	eval \$command \\
	| awk 'OFS=\"\\t\"{if(int(\$3)<int(\$2)){print \$1,\$3,\$2,\$4,\$5,\$6}else{print \$0}}' \\
	| sort -k1,1 -k2,2n -k3,3n -T . - \\
	| gzip > \${cellline}.\${well}.\${outidx}.CGread.bed.gz


    """
    //	| awk 'OFS=\"\\t\"{if(\$5>0)print\$0}' \\

}

process BED_MET_ATLAS_OVERLAP { 
	// Take per cell bed file of reads, binarize to methylated or not, and overlap with feature set
	publishDir "${params.scalemethylout}/postprocessing/read_bed", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'cnv'

	input:
		path cg_read
		path ref
	output:
		path("*.metatlas.bed.gz")
	script:
	"""
	feat_name=\"metatlas\"

	bedtools intersect -a ${ref} -b ${cg_read} -wao \\
	| awk '{ gsub(\"[.]\", \"0\"); print \$0}' \\
	| bedtools groupby -g 4 -c 9,10 -o sum \\
	| sort -k1,1 \\
	| awk 'OFS=\"\\t\"{if(\$3/\$2>0.9){print \$0,1}else{print \$0,0}}' \\
	| awk 'OFS=\"\\t\"{if(int(\$2)==0){print \$1,\$2,\$3,\".\"}else{print \$0}}' \\
	| gzip > ${cg_read.baseName}.\${feat_name}.bed.gz

	"""
}


process BED_MET_ATLAS_SUMMARY { 
	// Take per cell bed file of reads, binarize to methylated or not, and overlap with feature set
	publishDir "${params.scalemethylout}/postprocessing/read_bed", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'cnv'

	input:
		path metatlas_bed
	output:
		path("*read_met.tsv.gz")
	script:
	"""
	python /src/bed_merge.py \\
	--feat_name met_atlas \\
	--working_dir "."
	
	"""
}

// CNV BLOCK //
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


process MAKE_5KB_BED {
	//Generate transcript bed from ScaleMethyl genome.txt file
	//Add Blacklist filter??
	//requires bedtools
	
	label 'cnv'

	input:
		path ref
	output:
		path("genome_windows.5kb.bed")

	script:
	"""
	grep -v \"^K\" ${ref} \\
	| grep -v \"^G\" \\
	| awk \'OFS=\"\\t\" {print \"chr\"\$1,\$2}\' > genome.filt.txt
	
	bedtools makewindows -w 5000 -g genome.filt.txt \\
	| awk \'OFS=\"\\t\" {print \$1,\$2,\$3,\$1\"_\"\$2\"_\"\$3}\' > genome_windows.5kb.bed
	"""
}	


process MAKE_METATLAS_BED {
	label 'cnv'

	input:
		path ref
	output:
		path("${ref}")

	script:
	"""
	touch ${ref}

	"""
}	



//process MAKE_TF_MOTIF_BED {}
process SUMMARIZE_CG_OVER_BED {
	//publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true, pattern: "*.tsv.gz"
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'
	maxForks 20 //limiting maxForks for memory overlap (especially for atac peaks/ 5kb windows)

	input:
		path cov_folder
		path hundokb_bed
		val feat_name
		val min_cg
	output:
		path("*.tsv.gz")
	script:
	"""
	python /src/mc_cov_feature_summary.py \\
	--features ${hundokb_bed} \\
	--feat_name ${feat_name} \\
	--min_cg ${min_cg} \\
	--cov_folder ${cov_folder} \\
	--cpus 1
	"""
}



process MERGED_SUMMARIES {
	publishDir "${params.scalemethylout}/cg_dataframes", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.src_dir}:/src/,${params.scalemethylout}"
	label 'celltype'

	input:
		path cg_summaries
		val feat_name
	output:
		path("*merged.tsv.gz")
	script:
	"""
	#concatenate total count
	set -- *${feat_name}.total_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > total_count.${feat_name}.merged.tsv.gz

	#concatentate mc count
	set -- *${feat_name}.mc_count.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_count.${feat_name}.merged.tsv.gz

	#concatentate mc rate
	set -- *${feat_name}.mc_simplerate.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_simplerate.${feat_name}.merged.tsv.gz

	#concatenate rate estimates
	set -- *${feat_name}.mc_posteriorest.tsv.gz
	{
	zcat "\$1"; shift
	for file do
	    zcat "\$file" | sed '1d'
	done
	} | gzip > mc_posteriorest.${feat_name}.merged.tsv.gz
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
	/* SET CHANNELS */
		def bams = 
		Channel.fromPath("${params.scalemethylout}/bamDeDup/*/*bam")

		def fastqs = 
		Channel.fromFilePairs("${params.scalemethylout}/fastq/*/*_R{1,2}_*.fastq.gz",flat: true)

		def covs = 
		Channel.fromPath("${params.scalemethylout}/cg_sort_cov/*/*chroms.sort/", type: 'dir')

		def meta_data = 
		Channel.fromPath("${params.scalemethylout}/report/*/csv/*.passingCellsMapMethylStats.csv")
		met_in = meta_data | collect

	/* RUN CNV PROFILING ON CNVS */	
		fastqs \
		| SPLIT_GROUPED_FASTQ \
		| flatten | toSortedList | flatten | collate( 2 ) \
		| BISMARK_REALIGN \
		| BISMARK_TO_ALLCOOL \
		| collect \
		| ALLCOOL_TO_DATASET

		

	// if( "${params.cnv_call}" == "true" ) {
	//  		sorted_bams = 
	//  		SPLIT_GROUPED_BAM(bams)

	//  		scbam_dir = sorted_bams | collect
	//  		RUN_COPYKIT(scbam_dir)


	//  		//Summarize CG over reads
	//  		cg_beds=sorted_bams | flatten | BAM_TO_BED_MET
	//  		BED_MET_ATLAS_OVERLAP(cg_beds,"${params.metatlas_bed}") | collect | BED_MET_ATLAS_SUMMARY

	//  	}

	//  /* GENERATE CLUSTERS AND GENE SUMMARY WINDOWS FOR CELL TYPE ANALYSIS */

	//  	//Summarize CG at bp specificity
	//  	if( "${params.feat}" == "5kb" ) {
	//  		fiftykb_bed = MAKE_5KB_BED("${params.genome_length}")
	//  		fiftykb_out =
	//  		SUMMARIZE_CG_OVER_BED(covs, fiftykb_bed,"${params.feat}",1) \
	//  		| collect
	//  		MERGED_SUMMARIES(fiftykb_out,"${params.feat}")
	//  	}

	//  	else if( "${params.feat}" == "genebody" ) {
	//  		gene_bed = MAKE_TRANSCRIPT_BED("${params.genes_bed}")
	//  		genebody_out = 
	//  		SUMMARIZE_CG_OVER_BED(covs, gene_bed,"${params.feat}",5) \
	//  		| collect
	//  		MERGED_SUMMARIES(genebody_out,"${params.feat}")
	//  	}

	//  	else if( "${params.feat}" == "metatlas" ) {
	//  		metatlas_bed = MAKE_METATLAS_BED("${params.metatlas_bed}")
	//  		metatlas_out = 
	//  		SUMMARIZE_CG_OVER_BED(covs, metatlas_bed,"${params.feat}",1) \
	//  		| collect
	//  		MERGED_SUMMARIES(metatlas_out,"${params.feat}")
	//  	}

	//  	else if( "${params.feat}" == "atacpeaks" ) {
	//  		atac_bed = MAKE_METATLAS_BED("${params.multiome_atac_bed}")
	//  		atac_out = 
	//  		SUMMARIZE_CG_OVER_BED(covs, atac_bed,"${params.feat}",1) \
	//  		| collect
	//  		MERGED_SUMMARIES(atac_out,"${params.feat}")
	//  	}
	//  	//MAKE_FINAL_SEURATOBJ(hundokb_out,genebody_out,met_in)

	//  	//tf_bed = MAKE_TF_BED

}


/*


RUNNING MANUALLY
bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash

#to remove metatlas process scratch
rf -rf $(find /rsrch4/scratch/genetics/rmulqueen/met_work/ -type f -name "*read.bed.gz" | awk 'FS="/" {print $1"/"$2"/"$3}' -)

#set up environment variables 
export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact"
export srcDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src"
export refDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref"
export sif="${srcDir}/copykit.sif"

module load singularity 
#call nextflow
#nextflow ${srcDir}/singlecell_met_nf.groovy \
#--refdir $refDir \
#-w ${SCRATCH}/met_work \
#--scalemethylout ${projDir}/240205_RMMM_scalebiotest2 \
#-resume

export sif="${srcDir}/scmetR.sif"
export sif="${srcDir}/copykit.sif"
export sif="${srcDir}/allcool.sif"


cd /rsrch4/scratch/genetics/rmulqueen/met_work/97/9842c46c89915a3a739d9f2bf5ec29
singularity shell \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src:/src/ \
--bind /rsrch4/scratch/genetics/rmulqueen/met_work \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2 \
--bind /rsrch4/scratch/genetics/rmulqueen/work \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref/:/ref/ \
$sif

bsub -Is -W 4:00 -q transfer -n 10 -M 10 -R rusage[mem=10] /bin/bash
module load bedtools
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM5_ESM.zip
unzip 41586_2022_5580_MOESM5_ESM.zip -d met_atlas

conda activate bedtools-2.30.0
merge_bed () {
	in=$1
	celltype=$(printf '%s\n' "$in" | cut -d. -f1)
	awk -v outcell="$celltype" 'OFS ="\t" {print $1,$2,$3,$6,$7,outcell}' $in | tail -n +2
}
# add header and sort bed file, merge overlapping sites
export -f merge_bed

#header=$(echo chr$'\t'start$'\t'end$'\t'feat_description$'\t'gene$'\t'celltype)
(parallel merge_bed ::: *bed | sort -k 1,1 -k2,2n -k3,3n - | grep "^chr[1-9|X]" | bedtools merge ) > ../met_atlas.hg19.bed

#perform liftOver
#set up
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+wrx liftOver
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

#run
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref
/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src/liftOver \
met_atlas.hg19.bed \
hg19ToHg38.over.chain.gz \
met_atlas.hg38.bed \
met_atlas.hg19.unmapped.bed

cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref
awk 'OFS="\t" {print $1,$2,$3,$1"_"$2"_"$3}' met_atlas.hg38.bed > met_atlas.hg38.feat.bed


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
