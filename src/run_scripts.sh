#!/bin/bash

##1. FASTQC run for quality of sequencing:
fastqc -f fastq Sample_file.fastq## Step 1 (Python) - Calculating all the CREs between the organs and their overlaps

##2. Trimming the poor-quality reads and the adapters: 
trim_galore -q 28 --paired --gzip --phred33 --stringency 3 Sample_fastq1.fq.gz Sample_fastq2.fq.gz

##3. Alignment using STAR:
STAR --runThreadN 8 --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat -- genomeDir /path/to/genome/star/ --readFilesIn infile.fq.gz --outFileNamePrefix mapped_

##4. Convert sam to bam format:
samtools view -bSo OutFile.bam mapped_Infile.sam

##5. Sorting the aligned coordinates:
samtools sort InFlie.bam OutFile.sorted

##6. Indexing the bam file:
samtools index InFile.sorted.bam

##7. Remove duplicates:
samtools rmdup -s InFile.sorted.bam OutFile.sorted_dedup.bam

##8. Peak calling using macs2:
macs2 callpeak -t InFileH3K4me1.sorted_dedup.bam -c InFileInput.sorted_dedup.bam --format B