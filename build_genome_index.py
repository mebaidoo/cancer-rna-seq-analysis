import os
import subprocess

# Set paths
genome_dir = "genome_indices/STAR_R64-1-1"
fasta_url = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-60/fungi/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gtf_url = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-60/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.60.gtf.gz"
fastq_file = os.path.join("fastq_files", "SRR002058_1.fastq")

# Ensure the genome directory exists
os.makedirs(os.path.expanduser(genome_dir), exist_ok=True)

# Download the FASTA and GTF files
def download_file(url, output_dir):
    filename = os.path.join(output_dir, os.path.basename(url))
    if not os.path.exists(filename):
        print(f"Downloading: {url}")
        subprocess.run(f"wget {url} -O {filename}", shell=True, check=True)
    else:
        print(f"File already exists: {filename}")
    return filename

fasta_file = download_file(fasta_url, genome_dir)
gtf_file = download_file(gtf_url, genome_dir)

# Unzip the files (if not already unzipped)
def unzip_file(file_path):
    if file_path.endswith(".gz"):
        print(f"Unzipping: {file_path}")
        subprocess.run(f"gunzip -f {file_path}", shell=True, check=True)
        return file_path.replace(".gz", "")
    return file_path

fasta_file = unzip_file(fasta_file)
gtf_file = unzip_file(gtf_file)

# Get the read length from the first sequence in a FASTQ file (to use in sjdbOverhang calculation)
def get_read_length(fastq_file):
    
    with open(fastq_file, "r") as f:
        # FASTQ format: 4 lines per read, the sequence is on the 2nd line
        for i, line in enumerate(f):
            if i == 1:
                return len(line.strip())  # Return length of the sequence
    return None

# Build the STAR index with automatic sjdbOverhang calculation
def build_star_index(genome_dir, fasta_file, gtf_file, fastq_file):
    print("Determining read length from FASTQ...")
    read_length = get_read_length(fastq_file)
    
    if read_length is None:
        raise ValueError("Could not determine read length from FASTQ file.")
    
    # Calculate sjdbOverhang (STAR recommends sjdbOverhang = read_length - 1)
    sjdb_overhang = read_length - 1
    print(f"Detected read length: {read_length}, using sjdbOverhang: {sjdb_overhang}")

    print("Building STAR index...")
    command = f"""
    STAR --runThreadN 4 --runMode genomeGenerate \
         --genomeDir {genome_dir} \
         --genomeFastaFiles {fasta_file} \
         --sjdbGTFfile {gtf_file} \
         --sjdbOverhang {sjdb_overhang} \
         --genomeSAindexNbases 10
    """
    subprocess.run(command, shell=True, check=True)

build_star_index(genome_dir, fasta_file, gtf_file, fastq_file)
print("Index building complete.")