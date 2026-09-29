#!/bin/bash
# tEPOR gene editing study: RNA-seq processing
# Analysis parameters are retained from the original study script.
# Set the path variables below before running; requires Bash on Linux.
set -euo pipefail

timestamp() {
    date +"[%Y-%m-%d %H:%M:%S]"
}

# Sample identifiers
samples=(K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18 K19 K20 K21 K22 K23 K24)

# Configurable input, reference, and output paths
BASE_DIR="${BASE_DIR:-/path/to/RNAseqK1-24}"
TRIM_DIR="$BASE_DIR/trimmed"
FASTQC_DIR="$BASE_DIR/fastqc_trimmed"
STAR_DIR="$BASE_DIR/star_aligned"
COUNTS_OUTPUT="$BASE_DIR/gene_counts.txt"

ADAPTERS="${ADAPTERS:-/path/to/TruSeq3-PE.fa}"
STAR_INDEX="${STAR_INDEX:-/path/to/star_index}"
GTF="${GTF:-/path/to/gencode.v40.annotation.gtf}"

# Create output directories
mkdir -p "$TRIM_DIR" "$FASTQC_DIR" "$STAR_DIR" "$BASE_DIR/logs"

echo "$(timestamp) Starting Trimmomatic trimming..."
for sample in "${samples[@]}"
do
    echo "$(timestamp) Trimming $sample..."
    trimmomatic PE -threads 4         "$BASE_DIR/${sample}_1.fq.gz" "$BASE_DIR/${sample}_2.fq.gz"         "$TRIM_DIR/${sample}_1_paired.fq.gz" "$TRIM_DIR/${sample}_1_unpaired.fq.gz"         "$TRIM_DIR/${sample}_2_paired.fq.gz" "$TRIM_DIR/${sample}_2_unpaired.fq.gz"         ILLUMINACLIP:${ADAPTERS}:2:30:10         LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36         &> "$BASE_DIR/logs/${sample}_trimmomatic.log"
done
echo "$(timestamp) Trimmomatic step completed."

echo "$(timestamp) Running FastQC..."
for sample in "${samples[@]}"
do
    fastqc -t 4         "$TRIM_DIR/${sample}_1_paired.fq.gz"         "$TRIM_DIR/${sample}_2_paired.fq.gz"         -o "$FASTQC_DIR"
done

echo "$(timestamp) Running MultiQC..."
multiqc "$FASTQC_DIR" -o "$FASTQC_DIR"

echo "$(timestamp) Starting STAR alignment..."
for sample in "${samples[@]}"
do
    echo "$(timestamp) Aligning $sample..."
    STAR --runThreadN 8         --genomeDir "$STAR_INDEX"         --readFilesIn "$TRIM_DIR/${sample}_1_paired.fq.gz" "$TRIM_DIR/${sample}_2_paired.fq.gz"         --readFilesCommand zcat         --outFileNamePrefix "$STAR_DIR/${sample}_"         --outSAMtype BAM SortedByCoordinate         --outSAMunmapped Within         --quantMode GeneCounts         &> "$BASE_DIR/logs/${sample}_star.log"
done
echo "$(timestamp) STAR alignment completed."

echo "$(timestamp) Running featureCounts..."
BAM_FILES=$(ls $STAR_DIR/*_Aligned.sortedByCoord.out.bam)

featureCounts -T 8 -p -t exon -g gene_id   -a "$GTF"   -o "$COUNTS_OUTPUT"   $BAM_FILES   &> "$BASE_DIR/logs/featurecounts.log"

echo "$(timestamp) All steps completed successfully!"
