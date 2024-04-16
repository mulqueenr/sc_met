```bash
https://lhqing.github.io/ALLCools/start/installation.html
mamba env create -f allcools_env.yaml
conda activate allcools
mamba install -n allcools rpy2
mamba install -n allcools tpot xgboost dask-ml scikit-mdr skrebate


mamba env create -f mapping_env.yaml #yaml in /volumes/seq/projects/metACT/src
conda activate mapping
pip install cemba-data

bismark_genome_preparation --verbose /volumes/seq/projects/metACT/ref/grch38


```

Generate allc files from single cell bam files
```bash
conda activate allcools
cd /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/sc_bam
ref="/volumes/seq/projects/metACT/ref/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
allc_out="/volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/postprocessing/sc_allcool"
mkdir -p $allc_out

bam_convert() { 
	in=$1
	allcools bam-to-allc --bam_path ${in} --output_path $allc_out --reference_fasta $ref 
}

#correct cat /volumes/seq/projects/metACT/ref/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai so col1 has chr??
export -f bam_convert

parallel -j 1 bam_convert ::: $(ls *sorted.bam)

