#!/bin/bash
# RNA-seq Pipeline Script
# Author: Ajinkya Khilari
#
# This script processes RNA-seq data using a CSV file input.
# It supports quality control, trimming, alignment, SAM-to-BAM conversion,
# feature counting, and quantification.
#
# Pipeline Overview:
#
#     +---------+    +-------------+    +---------+    +------------+    +--------------+    +-------------+
#     | FastQC  | -> | Trimming    | -> | Alignment| -> | SAM-to-BAM | -> | FeatureCount | -> | Quantification |
#     +---------+    +-------------+    +---------+    +------------+    +--------------+    +-------------+
#
# Usage:
#   ./run_pipeline.sh -c samples.csv -g genome.fa -a annotation.gtf|annotation.gff -t transcripts.fa \
#                      -T [fastp|trim_galore] -l [STAR|HISAT2] -q [salmon|rsem] -p threads
#
# Example:
#   ./run_pipeline.sh -c samples.csv -g /path/to/genome.fa -a /path/to/annotation.gff \
#                      -t /path/to/transcripts.fa -T fastp -l STAR -q salmon -p 4

# Function to display usage information
usage() {
  echo "Usage: $0 -c samples.csv -g genome.fa -a annotation.gtf|annotation.gff -t transcripts.fa -T [fastp|trim_galore] -l [STAR|HISAT2] -q [salmon|rsem] -p threads"
  exit 1
}

# Parse command-line options
while getopts "c:g:a:t:T:l:q:p:" opt; do
    case $opt in
        c) CSV_FILE="$OPTARG" ;;
        g) GENOME_FASTA="$OPTARG" ;;
        a) ANNOTATION="$OPTARG" ;;
        t) TRANSCRIPTS_FA="$OPTARG" ;;
        T) TRIMMER="$OPTARG" ;;
        l) ALIGNER="$OPTARG" ;;
        q) QUANT="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check for required parameters
if [ -z "$CSV_FILE" ] || [ -z "$GENOME_FASTA" ] || [ -z "$ANNOTATION" ]; then
    echo "Error: Missing required parameters."
    usage
fi

# Set defaults if options not provided
: ${TRIMMER:="fastp"}
: ${ALIGNER:="STAR"}
: ${QUANT:="salmon"}
: ${THREADS:=4}

# Display chosen parameters
echo "--------------------------------------------"
echo "RNA-seq Pipeline"
echo "Author: Ajinkya Khilari"
echo "--------------------------------------------"
echo "Pipeline Configuration:"
echo "  CSV File:          $CSV_FILE"
echo "  Genome FASTA:      $GENOME_FASTA"
echo "  Annotation File:   $ANNOTATION"
echo "  Transcripts FASTA: $TRANSCRIPTS_FA"
echo "  Trimming Tool:     $TRIMMER"
echo "  Aligner:           $ALIGNER"
echo "  Quantification:    $QUANT"
echo "  Threads:           $THREADS"
echo "--------------------------------------------"
echo ""
echo "Pipeline Workflow:"
echo "  [1] FastQC  ->  [2] Trimming  ->  [3] Alignment  ->  [4] SAM-to-BAM  ->  [5] Feature Counting  ->  [6] Quantification"
echo "--------------------------------------------"
echo ""

# ---------------------------
# Define output directories and index directories (only for selected tools)
# ---------------------------
mkdir -p fastqc_reports

if [ "$TRIMMER" == "fastp" ] || [ "$TRIMMER" == "trim_galore" ]; then
    mkdir -p trimmed_reads
fi

if [ "$ALIGNER" == "STAR" ]; then
    mkdir -p star_alignments
elif [ "$ALIGNER" == "HISAT2" ]; then
    mkdir -p hisat2_alignments
fi

if [ "$QUANT" == "salmon" ]; then
    mkdir -p salmon_quant
elif [ "$QUANT" == "rsem" ]; then
    mkdir -p rsem_quant
fi

mkdir -p featureCounts_output

# ---------------------------
# Annotation Conversion (if provided file is GFF, convert to GTF)
# ---------------------------
if [[ "$ANNOTATION" == *.gff* ]]; then
    echo "[Annotation Conversion] Converting GFF to GTF..."
    GTF_FILE="${ANNOTATION%.*}.gtf"
    if [ ! -f "$GTF_FILE" ]; then
        gffread "$ANNOTATION" -T -o "$GTF_FILE"
    fi
    ANNOTATION_USED="$GTF_FILE"
else
    ANNOTATION_USED="$ANNOTATION"
fi
echo "Using annotation file: $ANNOTATION_USED"
echo "--------------------------------------------"

# ---------------------------
# Index Creation (only for tools specified)
# ---------------------------
if [ "$ALIGNER" == "STAR" ]; then
    STAR_INDEX="STAR_genome_index"
    if [ ! -d "$STAR_INDEX" ]; then
        echo "[Index Creation] Creating STAR index..."
        mkdir -p "$STAR_INDEX"
        STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir "$STAR_INDEX" \
             --genomeFastaFiles "$GENOME_FASTA" --sjdbGTFfile "$ANNOTATION_USED" --sjdbOverhang 100 --genomeSAindexNbases 11
    else
        echo "[Index Creation] STAR index found."
    fi
fi

if [ "$ALIGNER" == "HISAT2" ]; then
    HISAT2_INDEX="HISAT2_index"
    if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
        echo "[Index Creation] Creating HISAT2 index..."
        hisat2-build "$GENOME_FASTA" "$HISAT2_INDEX"
    else
        echo "[Index Creation] HISAT2 index found."
    fi
fi

if [ "$QUANT" == "salmon" ]; then
    SALMON_INDEX="salmon_index"
    if [ ! -d "$SALMON_INDEX" ]; then
        echo "[Index Creation] Creating Salmon index..."
        salmon index -t "$TRANSCRIPTS_FA" -i "$SALMON_INDEX" -k 31
    else
        echo "[Index Creation] Salmon index found."
    fi
fi

if [ "$QUANT" == "rsem" ]; then
    RSEM_REF="rsem_reference"
    if [ ! -f "${RSEM_REF}.grp" ]; then
        echo "[Index Creation] Creating RSEM reference..."
        rsem-prepare-reference --gtf "$ANNOTATION_USED" "$GENOME_FASTA" "$RSEM_REF"
    else
        echo "[Index Creation] RSEM reference found."
    fi
fi
echo "--------------------------------------------"

# ---------------------------
# Pipeline Processing for Each Sample
# ---------------------------
BAM_LIST="bam_list.txt"
rm -f "$BAM_LIST"

# Process each sample from the CSV file (CSV header: sample,fastq1,fastq2,strandness)
tail -n +2 "$CSV_FILE" | while IFS=, read -r SAMPLE FASTQ1 FASTQ2 STRANDNESS; do
    echo "====================================================="
    echo "Processing Sample: $SAMPLE"
    echo "-----------------------------------------------------"

    # Step 1: FastQC
    echo "[Step 1] Running FastQC on raw FASTQ files..."
    fastqc -o fastqc_reports "$FASTQ1" "$FASTQ2"

    # Step 2: Adapter Trimming
    echo "[Step 2] Adapter Trimming using $TRIMMER..."
    if [ "$TRIMMER" == "fastp" ]; then
        TRIM_R1="trimmed_reads/${SAMPLE}_R1_trimmed.fastq.gz"
        TRIM_R2="trimmed_reads/${SAMPLE}_R2_trimmed.fastq.gz"
        fastp -i "$FASTQ1" -I "$FASTQ2" \
              -o "$TRIM_R1" -O "$TRIM_R2" \
              --detect_adapter_for_pe --trim_poly_g --cut_front --cut_tail \
              --length_required 36 --thread $THREADS \
              --html "trimmed_reads/${SAMPLE}_fastp.html" --json "trimmed_reads/${SAMPLE}_fastp.json"
    elif [ "$TRIMMER" == "trim_galore" ]; then
        trim_galore --paired --cores $THREADS -o trimmed_reads "$FASTQ1" "$FASTQ2"
        TRIM_R1="trimmed_reads/$(basename ${FASTQ1%%.*})_val_1.fq.gz"
        TRIM_R2="trimmed_reads/$(basename ${FASTQ2%%.*})_val_2.fq.gz"
    else
        echo "Unknown trimming tool: $TRIMMER" && exit 1
    fi

    # Step 3: Alignment
    echo "[Step 3] Aligning reads using $ALIGNER..."
    if [ "$ALIGNER" == "STAR" ]; then
        ALIGN_OUT_PREFIX="star_alignments/${SAMPLE}_"
        STAR --runThreadN $THREADS --genomeDir "$STAR_INDEX" \
             --readFilesIn "$TRIM_R1" "$TRIM_R2" --readFilesCommand zcat \
             --outFileNamePrefix "$ALIGN_OUT_PREFIX"
        SAM_FILE="${ALIGN_OUT_PREFIX}Aligned.out.sam"
    elif [ "$ALIGNER" == "HISAT2" ]; then
        SAM_FILE="hisat2_alignments/${SAMPLE}.sam"
        hisat2 -p $THREADS -x "$HISAT2_INDEX" -1 "$TRIM_R1" -2 "$TRIM_R2" -S "$SAM_FILE"
    else
        echo "Unknown aligner: $ALIGNER" && exit 1
    fi

    # Step 4: SAM-to-BAM Conversion
    echo "[Step 4] Converting SAM to BAM..."
    BAM_FILE="${SAMPLE}.bam"
    samtools view -Sb "$SAM_FILE" > "$BAM_FILE"
    echo "$BAM_FILE" >> "$BAM_LIST"

    # Step 5: Set strandness for featureCounts
    echo "[Step 5] Setting strandness for feature counting..."
    if [ "$STRANDNESS" == "unstranded" ]; then
        FC_STRAND=0
    elif [ "$STRANDNESS" == "forward" ]; then
        FC_STRAND=1
    elif [ "$STRANDNESS" == "reverse" ]; then
        FC_STRAND=2
    else
        FC_STRAND=0
    fi

    # Step 6: Quantification
    echo "[Step 6] Quantification using $QUANT..."
    if [ "$QUANT" == "salmon" ]; then
        if [ "$STRANDNESS" == "unstranded" ]; then
            SALMON_LIB="A"
        elif [ "$STRANDNESS" == "forward" ]; then
            SALMON_LIB="ISF"
        elif [ "$STRANDNESS" == "reverse" ]; then
            SALMON_LIB="ISR"
        else
            SALMON_LIB="A"
        fi
        salmon quant -i "$SALMON_INDEX" -l "$SALMON_LIB" -1 "$TRIM_R1" -2 "$TRIM_R2" \
             -p $THREADS -o "salmon_quant/${SAMPLE}_salmon"
    elif [ "$QUANT" == "rsem" ]; then
        rsem-calculate-expression --paired-end --bam "$BAM_FILE" "$RSEM_REF" "rsem_quant/${SAMPLE}_rsem"
    else
        echo "Unknown quantification tool: $QUANT" && exit 1
    fi

    echo "Finished processing sample: $SAMPLE"
    echo "====================================================="
    echo ""
done

# ---------------------------
# Step 7: Run featureCounts on all BAM files (batch mode)
# ---------------------------
echo "[Step 7] Running featureCounts on all BAM files..."
mapfile -t BAM_FILES < "$BAM_LIST"
featureCounts -T $THREADS -s $FC_STRAND -a "$ANNOTATION_USED" -o featureCounts_output/all_counts.txt "${BAM_FILES[@]}"

echo "--------------------------------------------"
echo "RNA-seq Pipeline completed successfully."

