#!/bin/bash
module load samtools
module load bedtools
bam_file_name=$1
samtools sort -n $bam_file_name -o ${bam_file_name%.*}_samtools.bam
bedtools bamtobed -i ${bam_file_name%.*}_samtools.bam -bedpe  > ${bam_file_name%.*}.bed
awk 'BEGIN{FS="\t";OFS="\t"} {if (($6-$2)<=180 && ($6-$2)>=120) print "chr"$1, $2+int(($6-$2)/2), $2+int(($6-$2)/2)+1}' ${bam_file_name%.*}.bed | sort -k 1,1 > ${bam_file_name%.*}_fragment_mids_120_180.bed

awk 'BEGIN{FS="\t";OFS="\t"} {if (($6-$2)<=152 && ($6-$2)>=142) print "chr"$1, $2+int(($6-$2)/2), $2+int(($6-$2)/2)+1}' ${bam_file_name%.*}.bed | sort -k 1,1 > ${bam_file_name%.*}_fragment_mids_142_152.bed

awk 'BEGIN{FS="\t";OFS="\t"} {if (($6-$2)<=148 && ($6-$2)>=146) print "chr"$1, $2+int(($6-$2)/2), $2+int(($6-$2)/2)+1}' ${bam_file_name%.*}.bed | sort -k 1,1 > ${bam_file_name%.*}_fragment_mids_146_148.bed

#awk 'BEGIN{FS="\t";OFS="\t"} {if (($6-$2)==147) print "chr"$1, $2+int(($6-$2)/2), $2+int(($6-$2)/2)+1}' ${bam_file_name%.*}.bed | sort -k 1,1 > ${bam_file_name%.*}_fragment_mids_147.bed

samtools view -H ${bam_file_name%.*}_samtools.bam | grep @SQ|sed 's/@SQ\tSN:\|LN://g' | awk 'BEGIN{FS="\t"; OFS="\t"} {print "chr"$1, $2}'> ${bam_file_name%.*}_genome.txt

bedtools genomecov  -dz -i ${bam_file_name%.*}_fragment_mids_120_180.bed -g ${bam_file_name%.*}_genome.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2+1,$3}'  > ${bam_file_name%.*}_fragment_mids_coverage_120_180.bed

bedtools genomecov  -dz -i ${bam_file_name%.*}_fragment_mids_142_152.bed -g ${bam_file_name%.*}_genome.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2+1,$3}' > ${bam_file_name%.*}_fragment_mids_coverage_142_152.bed

bedtools genomecov  -dz -i ${bam_file_name%.*}_fragment_mids_146_148.bed -g ${bam_file_name%.*}_genome.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2+1,$3}' > ${bam_file_name%.*}_fragment_mids_coverage_146_148.bed

gzip -f ${bam_file_name%.*}_fragment_mids_coverage_120_180.bed
gzip -f ${bam_file_name%.*}_fragment_mids_coverage_142_152.bed
gzip -f ${bam_file_name%.*}_fragment_mids_coverage_146_148.bed

