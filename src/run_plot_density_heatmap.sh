#!/bin/bash

# Additional heatmaps for GFPpos_H3K4me1_peaks and GFPneg_H3K4me1_peaks

# Step 1: Generate matrix for GFPpos_H3K4me1_peaks using both GFPpos and GFPneg signals
echo "Generating matrix for GFPpos_H3K4me1_peaks..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R GFPpos_H3K4me1_peaks.bed \
    -S GFPpos_norm.bw GFPneg_norm.bw \
    --skipZeros \
    --missingDataAsZero \
    -o matrix_GFPpos_peaks.gz \
    --samplesLabel "GFP+" "GFP-" \
    -p 4

# Step 2: Generate matrix for GFPneg_H3K4me1_peaks using both GFPpos and GFPneg signals
echo "Generating matrix for GFPneg_H3K4me1_peaks..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R GFPneg_H3K4me1_peaks.bed \
    -S GFPpos_norm.bw GFPneg_norm.bw \
    --skipZeros \
    --missingDataAsZero \
    -o matrix_GFPneg_peaks.gz \
    --samplesLabel "GFP+" "GFP-" \
    -p 4

# Step 3: Create heatmap for GFPpos_H3K4me1_peaks
echo "Creating heatmap for GFPpos_H3K4me1_peaks..."
plotHeatmap -m matrix_GFPpos_peaks.gz \
    -out heatmap_GFPpos_H3K4me1_peaks.png \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --zMin -4 --zMax 4 \
    --xAxisLabel "Distance from peak center (bp)" \
    --refPointLabel "Center" \
    --regionsLabel "GFP+ H3K4me1 peaks" \
    --plotTitle "H3K4me1 Signal at GFP+ Peaks" \
    --heatmapHeight 12 \
    --heatmapWidth 6

# Step 4: Create heatmap for GFPneg_H3K4me1_peaks
echo "Creating heatmap for GFPneg_H3K4me1_peaks..."
plotHeatmap -m matrix_GFPneg_peaks.gz \
    -out heatmap_GFPneg_H3K4me1_peaks.png \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --zMin -4 --zMax 4 \
    --xAxisLabel "Distance from peak center (bp)" \
    --refPointLabel "Center" \
    --regionsLabel "GFP- H3K4me1 peaks" \
    --plotTitle "H3K4me1 Signal at GFP- Peaks" \
    --heatmapHeight 12 \
    --heatmapWidth 6

# Step 5: Generate individual heatmaps for each sample
# GFP+ signal at GFP+ peaks
echo "Creating heatmap for GFP+ signal at GFP+ peaks..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R GFPpos_H3K4me1_peaks.bed \
    -S GFPpos_norm.bw \
    --skipZeros \
    --missingDataAsZero \
    -o matrix_GFPpos_signal_at_GFPpos_peaks.gz \
    --samplesLabel "GFP+" \
    -p 4

plotHeatmap -m matrix_GFPpos_signal_at_GFPpos_peaks.gz \
    -out heatmap_GFPpos_signal_at_GFPpos_peaks.png \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --zMin -4 --zMax 4 \
    --xAxisLabel "Distance from peak center (bp)" \
    --refPointLabel "Center" \
    --regionsLabel "GFP+ H3K4me1 peaks" \
    --plotTitle "GFP+ H3K4me1 Signal at GFP+ Peaks" \
    --heatmapHeight 12 \
    --heatmapWidth 3

# GFP- signal at GFP- peaks
echo "Creating heatmap for GFP- signal at GFP- peaks..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R GFPneg_H3K4me1_peaks.bed \
    -S GFPneg_norm.bw \
    --skipZeros \
    --missingDataAsZero \
    -o matrix_GFPneg_signal_at_GFPneg_peaks.gz \
    --samplesLabel "GFP-" \
    -p 4

plotHeatmap -m matrix_GFPneg_signal_at_GFPneg_peaks.gz \
    -out heatmap_GFPneg_signal_at_GFPneg_peaks.png \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --zMin -4 --zMax 4 \
    --xAxisLabel "Distance from peak center (bp)" \
    --refPointLabel "Center" \
    --regionsLabel "GFP- H3K4me1 peaks" \
    --plotTitle "GFP- H3K4me1 Signal at GFP- Peaks" \
    --heatmapHeight 12 \
    --heatmapWidth 3

echo "Additional heatmaps generated successfully!"