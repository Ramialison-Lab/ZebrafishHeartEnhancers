# ZebrafishHeartEnhancers
 Zebrafish Heart Enhancers Discovery
## Description
Pipeline for ChIP-Seq analyses of the Zebrafish Heart Enhancers datasets.

Code for the manuscript: Gulrez Chahal^, Michael P. Eichenlaub^, Markus Tondl^, Michal Pawlak, Monika Mohenska, Lin Grimm, Lauren Bottrell, Mark Drvodelic, Sara Alaei, Jeannette Hallab, Lisa N. Waylen, Jose M. Polo, Cédric Blanpain, Nathan Palpant, Fernando Rossello, Minna-Liisa Änkö, Peter D. Currie, Benjamin M. Hogan, Cecilia Winata, Ekaterina Salimova, Hieu T. Nim*, Mirana Ramialison*. "An in vivo repertoire of zebrafish cardiomyocyte-specific cis-regulatory elements". 

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
## ChIP-seq analysis pipeline

Here we outline the steps involved in processing the ChIP-seq raw data to call out the peaks of interest and perform gene ontology. The following packages were installed to run the analysis:
Reference genome: Danio rerio (Genome assembly:GRCz10) Zv10
- [Quality check for sequencing data] (#fastqc)  
- [Trimming the poor-quality reads and the adapters] (#TrimGalore)
- [Alignment using STAR/BWA] (#star)
- [Convert sam to bam format, sorting, indexing and removing duplicates] (#samtools)
- [Calling peaks] (#macs2)


## <a name="fastqc">Quality check for sequencing data </a>
- Fastqc (Version 0.11.8)
````
fastqc -f fastq Sample_file.fastq
````
## Trimming the poor-quality reads and the adapters

- A q-cutoff of 28 and an adapter overlap (stringency) of 3nt are reasonable values for paired end files (here fwd reads are fastq1 and rev reads are fastq2) with the default illumina adapter, output will be a trimmed version of the file (eg. infile_trimmed.fq):
````
trim_galore -q 28 --paired --gzip --phred33 --stringency 3 Sample_fastq1.fq.gz Sample_fastq2.fq.gz
````
## <a name="star">Alignment using STAR/BWA </a>
- STAR (version 2.7.10b)
-Parameters:
  Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022
  genome build: zv10 , 
  indexing tool: STAR
  parameters: default
  masking options: none
````
wget http://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
mkdir -p STAR_index
STAR --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles Danio_rerio.GRCz10.dna.toplevel.fa \
     --runThreadN 8

STAR --runThreadN 8 --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat -- genomeDir /path/to/genome/star/ --readFilesIn infile_trimmed.fq.gz --outFileNamePrefix mapped_
```` 

- Alternatively BWA (version 0.7.18-r1243-dirt)
    genome build: zv10 
    indexing tool: bwa index
    parameters: default
    masking options: none 
````
wget http://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
bwa index Danio_rerio.GRCz10.dna.toplevel.fa


bwa mem -t 8 Danio_rerio.GRCz10.dna.toplevel.fa infile_trimmed.fq.gz > infile_aligned.sam
````

## <a name="samtools">Convert sam to bam format, sorting the aligned reads and indexing</a>

- sam file to bam file conversion:
````
samtools view -bSo OutFile.bam infile_aligned.sam
````
- Sorting the aligned coordinates:
````
samtools sort OutFile.bam OutFile.sorted.bam
````
- Indexing the bam file:
````
samtools index OutFile.sorted.bam
````
- Remove duplicates:
````
samtools rmdup -s OutFile.sorted.bam OutFile.sorted_dedup.bam
````

## <a name="macs2">Calling peaks</a>
Run the above steps for both the H3Kme1 and contol sample file.
-macs2 (version 2.1.0.20140616)
````
macs2 callpeak -t OutFileH3K4me1.sorted_dedup.bam -c OutFileInput.sorted_dedup.bam --format BAM -g 1.5e9 --broad –name OutFileH3K4me1_sorted_dedup 
````
The output peaks (potential cREs) were obtained in the BED file format for zebrafish genome ZV10 version. These were lifted over to zv9 to run gene ontology analysis using GREAT (version 3.0.0).
