# ChIP-seq Quality Assessment using ChIPQC

# Load required libraries
library(BiocManager)
BiocManager::install("ChIPQC")
BiocManager::install("GenomeInfoDbData")

# Install them all
BiocManager::install(txdb_packages)
BiocManager::install("TxDb.Hsapiens.UCSC.hg18.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene") 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm9.knownGene",
                       "TxDb.Mmusculus.UCSC.mm10.knownGene",
                       "TxDb.Mmusculus.UCSC.mm39.knownGene"))
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene")

install.packages("RCurl")
library(ChIPQC)


rm(list=ls())
register(SerialParam())

# Create directory structure if needed
dir.create("results", showWarnings = FALSE)
dir.create("ChIPQCreport", showWarnings = FALSE)

samples <- data.frame(
  SampleID = c("30kneg_Rep1", "30kpos_Rep1"),
  Tissue = c("sample", "sample"), 
  Factor = c("30kneg", "30kpos"),
  Condition = c("negative", "positive"),
  Replicate = c(1, 1),
  bamReads = c("mapped_GFPnegH3K4me1Aligned.out.sorted.bam", 
               "mapped_GFPposH3K4me1Aligned.out.sorted.bam"),
  ControlID = c("30kneginput", "30kposinput"),
  bamControl = c("mapped_30KnegInputAligned.out.sorted.bam",
                 "mapped_GFPposInputAligned.out.sorted.bam"),
  Peaks = c(NA, NA),  # You'll need to add peak files if available
  PeakCaller = c(NA, NA),
  stringsAsFactors = FALSE
)

print("Sample sheet:")
print(samples)

write.csv(samples, "samplesheet.csv", row.names = FALSE)

# Check if all BAM files exist
bam_files <- c(samples$bamReads, samples$bamControl)
missing_files <- bam_files[!file.exists(bam_files)]
if(length(missing_files) > 0) {
  cat("Missing files:\n")
  print(missing_files)
  stop("Please check file paths and ensure all BAM files are present")
}

# Check for BAM index files
bai_files <- paste0(bam_files, ".bai")
missing_bai <- bai_files[!file.exists(bai_files)]
if(length(missing_bai) > 0) {
  cat("Missing BAM index files:\n")
  print(missing_bai)
  cat("Creating missing index files...\n")
  
  # Create missing index files
  for(bam in bam_files[!file.exists(paste0(bam_files, ".bai"))]) {
    system(paste("samtools index", bam))
  }
}

cat("Creating ChIPQC object... This may take several minutes.\n")

# Create ChIPQC object
chipObj <- ChIPQC(samples)

# Generate summary statistics
cat("ChIPQC Summary:\n")
print(chipObj)

# Create ChIPQC report
cat("Generating ChIPQC report...\n")
ChIPQCreport(chipObj, 
             reportName="ChIP QC report: 30kneg vs 30kpos Analysis", 
             reportFolder="ChIPQCreport")

cat("Analysis complete! Check the ChIPQCreport folder for the HTML report.\n")

# Save the ChIPQC object for later use
saveRDS(chipObj, "chipqc_object.rds")

# Print some basic metrics
cat("\n=== Basic Quality Metrics ===\n")
cat("Read depths:\n")
print(QCmetrics(chipObj)$Reads)

if(ncol(QCmetrics(chipObj)) > 1) {
  cat("\nFragment lengths:\n")
  print(QCmetrics(chipObj)$FragLen)
  
  cat("\nRelative cross-correlation (RelCC):\n")
  print(QCmetrics(chipObj)$RelCC)
}

