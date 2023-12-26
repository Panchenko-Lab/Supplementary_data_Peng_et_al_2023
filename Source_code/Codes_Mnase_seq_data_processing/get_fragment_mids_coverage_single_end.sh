#!/bin/bash
module load samtools
module load bedtools
bam_file_name=$1
bedtools bamtobed -i $bam_file_name   > ${bam_file_name%.*}.bed
awk 'BEGIN{FS="\t";OFS="\t"}{if ($6=="+") print "chr"$1,$2+73,$2+74,$5; else if($6=="-") print "chr"$1,$3-73, $3-72, $5}' ${bam_file_name%.*}.bed | awk 'BEGIN{FS="\t";OFS="\t"}{if($2>0) print $1,$2,$3,$4}' > ${bam_file_name%.*}_fragment_mids.bed

samtools view -H $bam_file_name | grep @SQ|sed 's/@SQ\tSN:\|LN://g' | awk 'BEGIN{FS="\t"; OFS="\t"} {print "chr"$1, $2}' > ${bam_file_name%.*}_genome.txt

bedtools genomecov  -dz -i ${bam_file_name%.*}_fragment_mids.bed -g ${bam_file_name%.*}_genome.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2+1,$3}' > ${bam_file_name%.*}_fragment_mids_coverage.bed

gzip -f ${bam_file_name%.*}_fragment_mids_coverage.bed
