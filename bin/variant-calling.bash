#!/usr/bin/env bash

set -e
mkdir /workspace/nextflow_tutorial/bash-results
cd /workspace/nextflow_tutorial/bash-results

genome=/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in /workspace/nextflow_tutorial/data/trimmed_fastq/*_1.trim.fastq.gz;
do
    echo "working with file $fq1"; base=$(basename $fq1 _1.trim.fastq.gz);  echo "base name is $base" \

    fq1=/workspace/nextflow_tutorial/data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2=/workspace/nextflow_tutorial/data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam=/workspace/nextflow_tutorial/results/sam/${base}.aligned.sam
    bam=/workspace/nextflow_tutorial/results/bam/${base}.aligned.bam
    sorted_bam=/workspace/nextflow_tutorial/results/bam/${base}.aligned.sorted.bam
    raw_bcf=/workspace/nextflow_tutorial/results/bcf/${base}_raw.bcf
    variants=/workspace/nextflow_tutorial/results/vcf/${base}_variants.vcf
    final_variants=/workspace/nextflow_tutorial/results/vcf/${base}_final_variants.vcf

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
    vcfutils.pl varFilter $variants > $final_variants

done
