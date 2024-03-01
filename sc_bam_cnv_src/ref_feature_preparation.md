make transcript bed
```bash
gtf="/volumes/seq/projects/metACT/ref/grch38/filteredGTF/GRCh38_transcriptsOnly.gtf"

awk 'OFS="\t" {split($10,a,"\""); split($14,b,"\""); print "chr"$1,$4,$5,a[2],b[2]}' $gtf > GRCh38_transcripts.bed


select_longest() { 
	awk -v gene="$1" '{if($5==gene){print $0,$3-$2}}' GRCh38_transcripts.bed \
	| sort -k6,6n - \
	| tail -n 1
}

export -f select_longest

parallel -j 20 select_longest ::: $(awk '{print $5}' GRCh38_transcripts.bed | uniq) > GRCh38_transcripts.longest.bed

```

make 100kb bed
```bash
ref="/volumes/seq/projects/metACT/ref/grch38/genome.txt"
grep -v "^K" ${ref} \
| grep -v "^G" \
| awk 'OFS="\t" {print "chr"$1,$2}' > genome.filt.txt

bedtools makewindows -w 100000 -g genome.filt.txt | awk 'OFS="\t" {print $1,$2,$3,$1"_"$2"_"$3}' > genome_windows.100kb.bed
```
