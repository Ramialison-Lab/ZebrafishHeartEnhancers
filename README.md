# ZebrafishHeartEnhancers
 Zebrafish Heart Enhancers Discovery
## Description
Pipeline for analysis of Zebrafish Heart [ChIP-Seq raw data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1057907) to obtain the Zebrafish Heart Enhancers [datasets](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252150).

Code for the manuscript: Gulrez Chahal^, Michael P. Eichenlaub^, Markus Tondl^, Michal Pawlak, Monika Mohenska, Lin Grimm, Lauren Bottrell, Mark Drvodelic, Sara Alaei, Jeannette Hallab, Lisa N. Waylen, Jose M. Polo, Cédric Blanpain, Nathan Palpant, Fernando Rossello, Minna-Liisa Änkö, Peter D. Currie, Benjamin M. Hogan, Cecilia Winata, Ekaterina Salimova, Hieu T. Nim*, Mirana Ramialison*. "An in vivo repertoire of zebrafish cardiomyocyte-specific cis-regulatory elements". 

Languages: Bash. Operating systems: Windows, Linux, Mac OSX. Fully tested on Linux Ubuntu 20.04. 

## Downloading the raw data
The following fastq files can be downloaded from GEO using the code below:
- SRR27368649: GFP positive Input
- SRR27368650: GFP negative H3K4me1
- SRR27368651: GFP positive H3K4me1
- SRR27368652: GFP negative input
```
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR27368649 -o SRR27368649.fastq.gz
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR27368650 -o SRR27368650.fastq.gz
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR27368651 -o SRR27368651.fastq.gz
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR27368652 -o SRR27368652.fastq.gz
gunzip *.gz
```
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


## ChIP-seq analysis pipeline

Here we outline the steps involved in processing the ChIP-seq raw data to call out the peaks of interest and perform gene ontology. The following packages were installed to run the analysis:
Reference genome: Danio rerio (Genome assembly:GRCz10) Zv10
- [Quality check for sequencing data](#fastqc)
- [Trimming the poor-quality reads and the adapters](#TrimGalore)
- [Alignment using STAR/BWA](#star)
- [Convert sam to bam format, sorting, indexing and removing duplicates](#samtools)
- [Calling peaks](#macs2)


### <a name="fastqc">Quality check for sequencing data </a>
- Fastqc (Version 0.11.8)
````
fastqc -f fastq SRR27368651.fastq
````
### <a name="TrimGalore">Trimming the poor-quality reads and the adapter </a>

- A q-cutoff of 28 and an adapter overlap (stringency) of 3nt are reasonable values for paired end files (here fwd reads are fastq1 and rev reads are fastq2) with the default illumina adapter, output will be a trimmed compressed version of the file (eg. SRR27368651_trimmed.fq.gz):
````
trim_galore -q 28 --phred33 --fastqc --gzip --stringency 3 SRR27368651.fastq
````
### <a name="star">Alignment using STAR/BWA </a>
- STAR (version 2.7.10b)
-Parameters:
  --Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022
  --genome build: zv10 , 
  --indexing tool: STAR
  --parameters: default
  --masking options: none
````
wget http://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
mkdir -p STAR_index
STAR --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles Danio_rerio.GRCz10.dna.toplevel.fa \
     --runThreadN 8

STAR --runThreadN 8 --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat -- genomeDir /path/to/genome/star/ --readFilesIn SRR27368651_trimmed.fq.gz --outFileNamePrefix aligned_
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


bwa mem -t 8 Danio_rerio.GRCz10.dna.toplevel.fa SRR27368651_trimmed.fq.gz > SRR27368651_trimmed_aligned.sam
````

### <a name="samtools">Convert sam to bam format, sorting the aligned reads and indexing</a>

- sam file to bam file conversion:
````
samtools view -bSo SRR27368651_trimmed_aligned.bam SRR27368651_trimmed_aligned.sam
````
- Sorting the aligned coordinates:
````
samtools sort SRR27368651_trimmed_aligned.bam SRR27368651_trimmed_aligned_sorted.bam
````
- Indexing the bam file:
````
samtools index SRR27368651_trimmed_aligned_sorted.bam
````
- Remove duplicates:
````
samtools rmdup -s SRR27368651_trimmed_aligned_sorted.bam SRR27368651_trimmed_aligned_sorted_dedup.bam
````

### <a name="macs2">Calling peaks</a>
Run the above steps for both the H3Kme1 and the corresponding contol (Input) file.
-macs2 (version 2.1.0.20140616)
````
macs2 callpeak -t SRR27368651_trimmed_aligned_sorted_dedup.bam -c SRR27368649_trimmed_aligned_sorted_dedup.bam --format BAM -g 1.5e9 --broad –name OutFileH3K4me1_sorted_dedup 
````
The output peaks (potential cREs) were obtained in the BED file format for zebrafish genome ZV10 version. These were lifted over to zv9 to run gene ontology analysis using GREAT (version 3.0.0).
