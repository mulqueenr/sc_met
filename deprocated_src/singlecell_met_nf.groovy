nextflow.enable.dsl=2
// Script parameters
// Input parameters, user specified defaults


//TODO Correct blacklist regions for met clustering
//TODO Add outname parameter for CNV clones and clustering output

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
params.clustering_dim="chrom5k"
log.info """

		================================================
		    SCALEBIO CONTINUED PROCESSING v0.2
		================================================
		ScaleMethyl Pipeline Output : ${params.scalemethylout}
		Cell Readcount Filter : ${params.readcountfilter}
		NF Working Dir : ${workflow.workDir}
		Src directory : ${params.src_dir}
		Out directory : ${params.out_dir}
		Genes Bed : ${params.genes_bed}
		Met Atlas Bed : ${params.metatlas_bed}
		Multiome ATAC Bed : ${params.multiome_atac_bed}
		Clustering Dimension : ${params.clustering_dim}
		CNV Calling : ${params.cnv_call}
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
	publishDir "${params.scalemethylout}/postprocessing/clones", mode: 'copy', overwrite: true
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

// CLUSTERING BLOCK //
process SPLIT_GROUPED_FASTQ { 
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
	containerOptions "--bind ${params.bismark_ref}:/ref/"
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
	unset PYTHONPATH

	basename=\$(echo "${fq1.baseName}" | sed 's/.\\{6\\}\$//')

	cutadapt -a CTGCGACGGCTGC -A CTGTCTCTTAT -u 8 -U 8 -o out1.fq.gz -p out2.fq.gz ${fq1} ${fq2}

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

process ALLCOOL_TO_MCDS {
	publishDir "${params.scalemethylout}/postprocessing/bismark/sc_allcool", mode: 'copy', overwrite: true
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
	--regions transcripts /ref/allcool_regions/GRCh38_transcripts.longest.bed.gz \\
	--quantifiers transcripts count CGN \\
	--regions chrom5k 5000 \\
	--quantifiers chrom5k count CGN \\
	--quantifiers chrom5k hypo-score CGN cutoff=0.9 \\
	--regions chrom2k 2000 \\
	--quantifiers chrom2k count CGN \\
	--quantifiers chrom2k hypo-score CGN cutoff=0.9 
	"""
}

process MCDS_CLUSTERING {
	publishDir "${params.scalemethylout}/postprocessing/clusters", pattern:"*", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.refdir}:/ref/,${params.src_dir}:/src/"
	label 'allcool'
	cpus 50

	input:
		path(metadata)
		path(mcds_out)
		val var_dim
	output:
		path("*")
	script:
	"""
	#combine metadata, get header from first one, cat each one without header
	metafiles=\$( ls *passingCellsMapMethylStats.csv )
	header=\$(head -n 1 "\${metafiles[0]}")  
	( echo \$header ) && ( for i in *passingCellsMapMethylStats.csv ; do cat \$i | grep "pass" | tail -n +1 ; done) > cell_metadata.csv

	python /src/allcool_clustering.py \\
	--metadata cell_metadata.csv \\
	--mcds ./data.mcds \\
	--allcool_path . \\
	--var_dim ${var_dim} \\
	--outname sample
	"""
}


process ALLCOOL_CLUSTER_MERGE {
	publishDir "${params.scalemethylout}/postprocessing/clusters", pattern:"*", mode: 'copy', overwrite: true
	containerOptions "--bind ${params.refdir}:/ref/,${params.src_dir}:/src/"
	label 'allcool'
	cpus 50

	input:
		path(mcds_out)
	output:
		path("*")
	script:
	"""
	cluster_list=\$(zcat *ClusteringResults.csv.gz | tail -n +2 | awk '{print \$2}' | sort -k1,1n | uniq)
	
	for i in \$cluster_list; do
	cat zcat *ClusteringResults.csv.gz \\
	| awk -v clus=\$i '{if(\$2==clus) print \$3}' > tmp.txt
	
	xargs zcat < tmp.txt \\
	| awk 'OFS=\"\\t\" {print \$1,\$2-1,\$2,0,\$3,\$4,\$5,\$6,\$7}' \\
	| sort -T . --parallel=${task.cpus} -k1,1 -k2,2n -k3,3 \\
	| bedtools merge -i - -s -c 4,5,6,7 -o first,first,sum,sum \\
	| awk 'OFS=\"\\t\" {print \$1,\$3,\$4,\$5,\$6,\$7,1}' \\
	| bgzip > sample.Cluster_\${i}.allcool.gz

	tabix -@${task.cpus} -b 2 -e 2 -s 1 sample.Cluster_\${i}.allcool.gz; done
	"""
}
/*
process METATLAS_ALLCOOL_TO_MCDS {
	publishDir "${params.scalemethylout}/postprocessing/bismark/sc_allcool/metatlas_mcds", mode: 'copy', overwrite: true
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

	grep "^HBCA*" allc_table.tsv > allc_table2.tsv #only run HBCA samples

	allcools generate-dataset  \\
	--allc_table allc_table2.tsv \\
	--output_path data.mcds \\
	--chrom_size_path /ref/allcool_regions/genome.txt \\
	--obs_dim cell  \\
	--cpu ${task.cpus} \\
	--chunk_size 30 \\
	--regions METATLAS /ref/allcool_regions/met_atlas.hg38.feat.bed.gz \\
	--quantifiers METATLAS count CGN \\
	--quantifiers METATLAS hypo-score CGN cutoff=0.9 
	
	"""

	//--regions ATAC /ref/allcool_regions/multiome_bc.marker_atac_peaks.bed.gz \\
	//--quantifiers ATAC count CGN \\
	//--quantifiers ATAC hypo-score CGN cutoff=0.9
}


process METATLAS_ALLCOOL_TO_MCDS {
	publishDir "${params.scalemethylout}/postprocessing/bismark/sc_allcool/atac_mcds", mode: 'copy', overwrite: true
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

	grep "^HBCA*" allc_table.tsv > allc_table2.tsv #only run HBCA samples

	allcools generate-dataset  \\
	--allc_table allc_table2.tsv \\
	--output_path data.mcds \\
	--chrom_size_path /ref/allcool_regions/genome.txt \\
	--obs_dim cell  \\
	--cpu ${task.cpus} \\
	--chunk_size 30 \\
	//--regions ATAC /ref/allcool_regions/multiome_bc.marker_atac_peaks.bed.gz \\
	//--quantifiers ATAC count CGN \\
	//--quantifiers ATAC hypo-score CGN cutoff=0.9
	"""


}
*/

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

		def var_dim = Channel.of("${params.clustering_dim}")
		
	/* RUN CLUSTERING ON METHYLATION DATA */	
		allc=
		fastqs \
		| SPLIT_GROUPED_FASTQ \
		| flatten | toSortedList | flatten | collate( 2 ) \
		| BISMARK_REALIGN \
		| BISMARK_TO_ALLCOOL \
		| collect
		| ALLCOOL_TO_MCDS
		
		allc_out=
		MCDS_CLUSTERING(met_in,allc,var_dim) \
		| ALLCOOL_CLUSTER_MERGE

		
	/* RUN CNV PROFILING ON CNVS */	
		if( "${params.cnv_call}" == "true" ) {
	  		sorted_bams = SPLIT_GROUPED_BAM(bams) | collect | RUN_COPYKIT
	  	}


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

ADD OVERDISPERSION
MAD SCORE
BREADTH OF COVERAGE

*/
