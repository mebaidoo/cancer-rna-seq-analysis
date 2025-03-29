import os
import subprocess

# Set directories
fastq_dir = "fastq_files"  # Path to FASTQ files
qc_dir = "qc_reports"  # Output directory for FastQC and MultiQC reports
trimmed_dir = "trimmed_reads"  # Directory for trimmed reads
trimmed_qc_dir = "trimmed_qc_reports"
aligned_dir = "aligned"  # Directory for aligned BAM files
genome_dir = "genome_indices/STAR_R64-1-1" # Genome index path

# Ensure output directories exist
for directory in [trimmed_dir, trimmed_qc_dir]:
    os.makedirs(directory, exist_ok=True)

# Function to run shell commands
def run_command(command):
    try:
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

# 1. Quality Control with FastQC
def run_fastqc():
    run_command(f"fastqc {fastq_dir}/*.fastq -o {qc_dir}")

# 2. Summarize Quality Reports with MultiQC
def run_multiqc():
    run_command(f"multiqc {qc_dir} -o {qc_dir}")

# 3. Trim Reads (if necessary) using Trim Galore
def trim_reads():
    run_command(f"trim_galore --quality 20 --phred33 --output_dir {trimmed_dir} {fastq_dir}/*.fastq")

# 4. Function to perform QC on trimmed reads
def qc_trimmed_reads():
    # Step 1: Run FastQC on trimmed reads
    run_command(f"fastqc {trimmed_dir}/*.fq -o {trimmed_qc_dir}")
    
    # Step 2: Summarize QC reports with MultiQC
    run_command(f"multiqc {trimmed_qc_dir} -o {trimmed_qc_dir}")
            
# 5. Align Reads using STAR (for single-end data)
def align_reads():
    for fastq_file in os.listdir(trimmed_dir):
        if fastq_file.endswith("_1_trimmed.fq"):
            sample_id = fastq_file.split("_1_trimmed.fq")[0]
            read1 = os.path.join(trimmed_dir, fastq_file)
            output_prefix = os.path.join(aligned_dir, sample_id)

            run_command(f"STAR --runThreadN 4 --genomeDir {genome_dir} --readFilesIn {read1} "
                        f"--outFileNamePrefix {output_prefix}_ --outSAMtype BAM SortedByCoordinate")
            # If two reads (paired-end), add:
            # read2 = os.path.join(trimmed_dir, f"{sample_id}_2.fastq")
            # run_command(f"STAR --runThreadN 4 --genomeDir {genome_dir} --readFilesIn {read1} {read2} "
            #             f"--outFileNamePrefix {output_prefix}_ --outSAMtype BAM SortedByCoordinate")


# Main function to execute the pipeline
def main():
    print("Starting RNA-seq analysis...")

    # Step 1: Run FastQC on FASTQ files
    run_fastqc()

    # Step 2: Summarize QC reports with MultiQC
    run_multiqc()

    # Step 3: Trim reads if needed
    trim_reads()

    # Step 4: Run FastQC and MultiQC on trimmed reads to check their quality after trimming
    qc_trimmed_reads()

    # Step 5: Align reads using STAR
    align_reads()

    print("RNA-seq analysis complete.")

if __name__ == "__main__":
    main()