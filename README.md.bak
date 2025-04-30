# ZebrafishHeartEnhancers
 Zebrafish Heart Enhancers Discovery
## Description
Pipeline for ChIP-Seq analyses of the Zebrafish Heart Enhancers datasets.

Code for the manuscript: Markus Tondl^, Gulrez Chahal1^, Michael P. Eichenlaub1^, Michał Pawlak, Monika Mohenska, Lin Grimm, Lauren Bottrell, Mark Drvodelic, Sara Alaei, Jeannette Hallab, Lisa N. Waylen, Jose M. Polo, Cédric Blanpain, Eileen Furlong, Nathan Palpant, Ben Hogan, Cecilia Winata, Ekaterina Salimova, Hieu T. Nim*, Mirana Ramialison*. "An in vivo repertoire of zebrafish cardiomyocyte-specific cis-regulatory elements". 

Languages: Bash. Operating systems: Windows, Linux, Mac OSX. Fully tested on Linux Ubuntu 20.04. 

## Getting the Source Code

To get the source code, please click the "fork" button in the upper-right and then add this repo as an upstream source:

````
$ git clone <your_fork_of_the_repo> ppds
$ cd ppds
$ REPO=https://github.com/Ramialison-Lab/ZebrafishHeartEnhancers.git
$ git remote add upstream $REPO
````

To get new content, use 
````
$ git pull upstream master 
````

Usage:

```text
ChIP-seq analysis pipeline

Here we outline the steps involved in processing the ChIP-seq raw data to call out the peaks of interest and perform gene ontology. The following packages were installed to run the analysis:
Reference genome: Danio rerio (Genome assembly:GRCz10) Zv10
FastQC  (Version 0.11.8)
Cutadapt (version 1.16)
STAR (version 2.5.0a)
MACS2 (version 2.1.0.20140616)

1. FASTQC run for quality of sequencing:

fastqc -f fastq Sample_file.fastq

2. Trimming the poor-quality reads and the adapters: 
A q-cutoff of 28 and an adapter overlap (stringency) of 3nt are reasonable values for paired end files (here fwd reads are fastq1 and rev reads are fastq2) with the default illumina adapter:

trim_galore -q 28 --paired --gzip --phred33 --stringency 3 Sample_fastq1.fq.gz Sample_fastq2.fq.gz

3. Alignment using STAR:


STAR --runThreadN 8 --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat -- genomeDir /path/to/genome/star/ --readFilesIn infile.fq.gz --outFileNamePrefix mapped_


4. Convert sam to bam format:

samtools view -bSo OutFile.bam mapped_Infile.sam

5. Sorting the aligned coordinates:


samtools sort InFlie.bam OutFile.sorted

6. Indexing the bam file:

samtools index InFile.sorted.bam


7. Remove duplicates:


samtools rmdup -s InFile.sorted.bam OutFile.sorted_dedup.bam


8. Peak calling using macs2:


macs2 callpeak -t InFileH3K4me1.sorted_dedup.bam -c InFileInput.sorted_dedup.bam --format B

The output peaks (potential cREs) were obtained in the BED file format for zebrafish genome ZV10 version. These were lifted over to zv9 to run gene ontology analysis using GREAT (version 3.0.0).
