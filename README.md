# rna_seq
RNA seq analysis for illumina data


# RNA-seq Pipeline for Transcriptome Analysis

This repository contains an automated RNA-seq analysis pipeline designed for processing Illumina transcriptome data. The pipeline supports quality control, adapter trimming, alignment, SAM-to-BAM conversion, feature counting, and quantification using tools such as FastQC, fastp/Trim Galore, STAR/HISAT2, featureCounts, Salmon/RSEM, and more.


## Features

- **CSV-Driven Sample Input:**  
  Provide a CSV file with columns for sample name, R1 FASTQ, R2 FASTQ, and strandness.

- **Quality Control:**  
  Uses FastQC to assess raw read quality.

- **Adapter Trimming:**  
  Options for trimming with either fastp or Trim Galore.

- **Alignment Options:**  
  Aligns reads with STAR or HISAT2.  
  Automatically creates indices if they do not exist.

- **SAM-to-BAM Conversion:**  
  Converts SAM output from aligners to BAM format using samtools.

- **Feature Counting:**  
  Uses featureCounts (from the Subread package) to count reads per gene.

- **Quantification:**  
  Quantification options include Salmon (pseudoalignment) and RSEM (alignment-based).

- **Annotation Flexibility:**  
  Supports GTF or GFF annotations. GFF files are automatically converted to GTF using gffread.


## Installation

The pipeline is designed to run in a Conda environment. An example `environment.yml` file is provided to install all necessary software.

1. **Clone this repository:**
   
   git clone https://github.com/ajinkyakhilari/rna_seq.git
   
   cd rnaseq-pipeline

**Create the Conda Environment:**

conda env create -f environment.yml

conda activate rnaseq_env


**Usage**

Command-Line Input

The main pipeline script (run_pipeline.sh) accepts input parameters from the command line. The required options include:

-c: CSV file with sample information.

-g: Genome FASTA file.

-a: Annotation file (GTF or GFF).

-t: Transcripts FASTA file (for Salmon indexing).

-T: Trimming tool (fastp or trim_galore).

-l: Aligner (STAR or HISAT2).

-q: Quantification tool (salmon or rsem).

-p: Number of threads to use.

#### Example:

***bash run_pipeline.sh -c samples.csv -g /path/to/genome.fa -a /path/to/annotation.gff -t /path/to/transcripts.fa -T fastp -l STAR -q salmon -p 4***


**CSV File Format**

Your CSV file should contain a header with the following columns:

sample: Unique sample name

fastq1: Path to the first FASTQ file (R1)

fastq2: Path to the second FASTQ file (R2)

strandness: Strandness of the data (e.g., unstranded, forward, or reverse)

#### Example:

sample,fastq1,fastq2,strandness

Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz,unstranded

Sample2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz,reverse


**Directory Structure**

After running the pipeline, the following directories/files will be generated:

fastqc_reports/: Contains FastQC reports for raw FASTQ files

trimmed_reads/: Contains trimmed FASTQ files and QC reports from fastp/Trim Galore

star_alignments/ or hisat2_alignments/: Contains alignment outputs (SAM files)

salmon_quant/: Output directory for Salmon quantification (if selected)

rsem_quant/: Output directory for RSEM quantification (if selected)

counts.txt: Gene-level count matrix from featureCounts

bam_list.txt: Temporary file listing all BAM files processed



**Pipeline Workflow**

- Index Creation:

Builds STAR, HISAT2, Salmon, and RSEM indices/references if not already present

- Converts GFF annotation to GTF if necessary

- Quality Control:

Runs FastQC on the input FASTQ files

- Adapter Trimming:

Trims adapters and low-quality bases using the selected trimming tool

- Alignment:

Aligns trimmed reads using the selected aligner (STAR/HISAT2)

- Conversion:

Converts SAM files to BAM format using samtools

- Feature Counting:

Counts reads using featureCounts

- Quantification:

Performs transcript quantification with Salmon or RSEM


**Contributing**

Contributions to improve this pipeline are welcome! Feel free to open issues or pull requests.


**License**

This project is licensed under the GPL3.0 License. See the LICENSE file for details.


**Contact**

For questions or support, please open an issue in this repository.


