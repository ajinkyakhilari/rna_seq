#!/bin/bash
# Usage:
#   ./run_pipeline.sh -c samples.csv -g genome.fa -a annotation.gtf|annotation.gff -t transcripts.fa -T [fastp|trim_galore] -l [STAR|HISAT2] -q [salmon|rsem] -p threads
#
# Example:
#   ./run_pipeline.sh -c samples.csv -g /path/to/genome.fa -a /path/to/annotation.gff -t /path/to/transcripts.fa -T fastp -l STAR -q salmon -p 4

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
echo "CSV File:          $CSV_FILE"
echo "Genome FASTA:      $GENOME_FASTA"
echo "Annotation File:   $ANNOTATION"
echo "Transcripts FASTA: $TRANSCRIPTS_FA"
echo "Trimming Tool:     $TRIMMER"
echo "Aligner:           $ALIGNER"
echo "Quantification:    $QUANT"
echo "Threads:           $THREADS"

# ---------------------------
# Define directories for indices and outputs
# ---------------------------
# Define index directories (modify these paths as needed)
STAR_INDEX="STAR_genome_index"
HISAT2_INDEX="HISAT2_index"         # HISAT2 will create files like HISAT2_index.1.ht2 etc.
SALMON_INDEX="salmon_index"         # Salmon index directory
RSEM_REF="rsem_reference"           # RSEM reference prefix

# ---------------------------
# Annotation Conversion (if provided file is GFF, convert to GTF)
# ---------------------------
if [[ "$ANNOTATION" == *.gff* ]]; then
    echo "Annotation appears to be GFF. Converting to GTF using gffread..."
    GTF_FILE="${ANNOTATION%.*}.gtf"
    if [ ! -f "$GTF_FILE" ]; then
        gffread "$ANNOTATION" -T -o "$GTF_FILE"
    fi
    ANNOTATION_USED="$GTF_FILE"
else
    ANNOTATION_USED="$ANNOTATION"
fi
echo "Using annotation file: $ANNOTATION_USED"

# ---------------------------
# Create output directories
# ---------------------------
mkdir -p fastqc_reports trimmed_reads star_alignments hisat2_alignments salmon_quant rsem_quant

# ---------------------------
# Index Creation
# ---------------------------
echo "Checking and creating indices if necessary..."

# STAR index creation (requires annotation in GTF format)
if [ ! -d "$STAR_INDEX" ]; then
    echo "STAR index not found. Creating STAR index..."
    mkdir -p "$STAR_INDEX"
    STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir "$STAR_INDEX" \
         --genomeFastaFiles "$GENOME_FASTA" --sjdbGTFfile "$ANNOTATION_USED" --sjdbOverhang 100 --genomeSAindexNbases 11
else
    echo "STAR index found."
fi

# HISAT2 index creation
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "HISAT2 index not found. Creating HISAT2 index..."
    hisat2-build "$GENOME_FASTA" "$HISAT2_INDEX"
else
    echo "HISAT2 index found."
fi

# Salmon index creation (only if using Salmon quantification)
if [ "$QUANT" == "salmon" ]; then
    if [ ! -d "$SALMON_INDEX" ]; then
        echo "Salmon index not found. Creating Salmon index..."
        salmon index -t "$TRANSCRIPTS_FA" -i "$SALMON_INDEX" -k 31
    else
        echo "Salmon index found."
    fi
fi

# RSEM reference creation (only if using RSEM quantification)
if [ "$QUANT" == "rsem" ]; then
    if [ ! -f "${RSEM_REF}.grp" ]; then
        echo "RSEM reference not found. Creating RSEM reference..."
        rsem-prepare-reference --gtf "$ANNOTATION_USED" "$GENOME_FASTA" "$RSEM_REF"
    else
        echo "RSEM reference found."
    fi
fi

# ---------------------------
# Pipeline Processing for Each Sample
# ---------------------------
# Temporary file to collect BAM filenames for featureCounts
BAM_LIST="bam_list.txt"
rm -f "$BAM_LIST"

# Process each sample from the CSV file (assumes header exists; columns: sample,fastq1,fastq2,strandness)
tail -n +2 "$CSV_FILE" | while IFS=, read -r SAMPLE FASTQ1 FASTQ2 STRANDNESS; do
    echo "-----------------------------------------------------"
    echo "Processing sample: $SAMPLE"

    ## Step 1: FastQC on raw FASTQ files
    fastqc -o fastqc_reports "$FASTQ1" "$FASTQ2"

    ## Step 2: Adapter Trimming
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

    ## Step 3: Alignment
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

    ## Step 4: SAM-to-BAM Conversion
    BAM_FILE="${SAMPLE}.bam"
    samtools view -Sb "$SAM_FILE" > "$BAM_FILE"
    echo "$BAM_FILE" >> "$BAM_LIST"

    ## Step 5: Feature Counting (strandness mapping: unstranded=0, forward=1, reverse=2)
    if [ "$STRANDNESS" == "unstranded" ]; then
        FC_STRAND=0
    elif [ "$STRANDNESS" == "forward" ]; then
        FC_STRAND=1
    elif [ "$STRANDNESS" == "reverse" ]; then
        FC_STRAND=2
    else
        FC_STRAND=0
    fi

    ## Step 6: Quantification Options
    if [ "$QUANT" == "salmon" ]; then
        # Define library type based on strandness: unstranded -> A, forward -> ISF, reverse -> ISR
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
done

## Step 7: Run featureCounts on all BAM files (batch mode)
mapfile -t BAM_FILES < "$BAM_LIST"
featureCounts -T $THREADS -s $FC_STRAND -a "$ANNOTATION_USED" -o counts.txt "${BAM_FILES[@]}"

echo "Pipeline completed successfully."
