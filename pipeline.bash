#!/bin/zsh/
# ChIP-Seq Bash Pipeline

# Quality assessment of the fastq files
fastqc /path/to/fastq/files/*.fastq.gz -t 16 --outdir=/path/to/output/files

# Make Bowtie2 Index from Mus_musculus.GRCm39.dna.primary_assembly.fa.gz (Ensemble)
bowtie2-build /path/to/Mus_musculus.GRCm39.dna.primary_assembly.fa /path/to/Mus_musculus.GRCm39.dna.primary_assembly

# Alignment to Mus_musculus.GRCm39.dna.primary_assembly.fa.gz (Ensemble) using Bowtie2
bowtie2 -p 16 -q --local -x /path/to/Mus_musculus.GRCm39.dna.primary_assembly -U /path/to/KI-Ischemia_CEBP-beta_mm_i43.fastq.gz -S /path/to/KI-Ischemia_CEBP-beta_mm_i43.sam &> /path/to/KI-Ischemia_CEBP-beta_mm_i43.log
bowtie2 -p 16 -q --local -x /path/to/Mus_musculus.GRCm39.dna.primary_assembly -U /path/to/KI-Sham_CEBP-beta_mm_i40.fastq.gz -S /path/to/KI-Sham_CEBP-beta_mm_i40.sam &> /path/to/KI-Sham_CEBP-beta_mm_i40.log
bowtie2 -p 16 -q --local -x /path/to/Mus_musculus.GRCm39.dna.primary_assembly -U /path/to/Pooled_Input_mm_i64.fastq.gz -S /path/to/Pooled_Input_mm_i64.sam &> /path/to/Pooled_Input_mm_i64.log
bowtie2 -p 16 -q --local -x /path/to/Mus_musculus.GRCm39.dna.primary_assembly -U /path/to/WT-Ischemia_CEBP-beta_mm_i96.fastq.gz -S /path/to/WT-Ischemia_CEBP-beta_mm_i96.sam &> /path/to/WT-Ischemia_CEBP-beta_mm_i96.log
bowtie2 -p 16 -q --local -x /path/to/Mus_musculus.GRCm39.dna.primary_assembly -U /path/to/WT-Sham_CEBP-beta_mm_i95.fastq.gz -S /path/to/WT-Sham_CEBP-beta_mm_i95.sam &> /path/to/WT-Sham_CEBP-beta_mm_i95.log

# Sam to Bam Conversion
samtools view -h -S -b -o /path/to/KI-Ischemia_CEBP-beta_mm_i43.bam /path/to/KI-Ischemia_CEBP-beta_mm_i43.sam
samtools view -h -S -b -o /path/to/KI-Sham_CEBP-beta_mm_i40.bam /path/to/KI-Sham_CEBP-beta_mm_i40.sam
samtools view -h -S -b -o /path/to/Pooled_Input_mm_i64.bam /path/to/Pooled_Input_mm_i64.sam
samtools view -h -S -b -o /path/to/WT-Ischemia_CEBP-beta_mm_i96.bam /path/to/WT-Ischemia_CEBP-beta_mm_i96.sam
samtools view -h -S -b -o /path/to/WT-Sham_CEBP-beta_mm_i95.bam /path/to/WT-Sham_CEBP-beta_mm_i95.sam 

# Sorting BAM files by genomics coordinates
sambamba sort -t 16 -o /path/to/KI-Ischemia_CEBP-beta_mm_i43_sorted.bam /path/to/KI-Ischemia_CEBP-beta_mm_i43.bam
sambamba sort -t 16 -o /path/to/KI-Sham_CEBP-beta_mm_i40_sorted.bam /path/to/KI-Sham_CEBP-beta_mm_i40.bam
sambamba sort -t 16 -o /path/to/Pooled_Input_mm_i64_sorted.bam /path/to/Pooled_Input_mm_i64.bam
sambamba sort -t 16 -o /path/to/WT-Ischemia_CEBP-beta_mm_i96_sorted.bam /path/to/WT-Ischemia_CEBP-beta_mm_i96.bam
sambamba sort -t 16 -o /path/to/WT-Sham_CEBP-beta_mm_i95_sorted.bam /path/to/WT-Sham_CEBP-beta_mm_i95.bam

# Filtering Uniquely Mapping Reads
sambamba view -h -t 16 -f bam -F "[XS] == null and not unmapped  and not duplicate" /path/to/KI-Ischemia_CEBP-beta_mm_i43_sorted.bam > /path/to/KI-Ischemia_CEBP-beta_mm_i43_sorted_filtered.bam
sambamba view -h -t 16 -f bam -F "[XS] == null and not unmapped  and not duplicate" /path/to/KI-Sham_CEBP-beta_mm_i40_sorted.bam > /path/to/KI-Sham_CEBP-beta_mm_i40_sorted_filtered.bam
sambamba view -h -t 16 -f bam -F "[XS] == null and not unmapped  and not duplicate" /path/to/Pooled_Input_mm_i64_sorted.bam > /path/to/Pooled_Input_mm_i64_sorted_filtered.bam
sambamba view -h -t 16 -f bam -F "[XS] == null and not unmapped  and not duplicate" /path/to/WT-Ischemia_CEBP-beta_mm_i96_sorted.bam > /path/to/WT-Ischemia_CEBP-beta_mm_i96_sorted_filtered.bam
sambamba view -h -t 16 -f bam -F "[XS] == null and not unmapped  and not duplicate" /path/to/WT-Sham_CEBP-beta_mm_i95_sorted.bam > /path/to/WT-Sham_CEBP-beta_mm_i95_sorted_filtered.bam

# Check alignment
multiqc -d /path/to/ChIP-Seq/directory -o /path/to/multiqc/output
