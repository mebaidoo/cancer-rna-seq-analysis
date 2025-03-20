import os
import pandas as pd
from Bio import Entrez
import re
import subprocess

# Retrieve fastq data files from SRA using GEO accession ID, querying with NCBI's Entrez

# Set NCBI email
Entrez.email = "baidoo394@gmail.com"

# Function to fetch GSM (Geo Sample) IDs from a GSE (Geo Series) accession ID
def get_gsm_ids_from_gse(gse_id):
    gsm_ids = []
    try:
        # Search gds (GEO datasets) for the given GSE accession
        handle = Entrez.esearch(db="gds", term=gse_id, retmax=1000)
        record = Entrez.read(handle)
        handle.close()

        # Fetch detailed information
        if record["IdList"]:
            gse_summary = Entrez.esummary(db="gds", id=record["IdList"][0])
            gse_info = Entrez.read(gse_summary)
            gse_summary.close()

            # Extract GSM IDs from the samples list
            for sample in gse_info[0]["Samples"]:
                gsm_ids.append(sample["Accession"])
    except Exception as e:
        print(f"Error: {e}")
    return gsm_ids

# Function to fetch SRA (Sequence Read Archive) IDs from a GSM ID
def get_sra_ids_from_gsm(gsm_ids):
    gsm_to_sra = {}

    for gsm_id in gsm_ids:
        sra_ids = []
        try:
            # Search for the GSM record in the SRA database
            handle = Entrez.esearch(db="sra", term=gsm_id)
            record = Entrez.read(handle)
            handle.close()

            # Get SRA accession numbers from the search results
            for sra_id in record['IdList']:
                summary_handle = Entrez.esummary(db="sra", id=sra_id)
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()

                # Extract SRA IDs from the run information
                runs = summary_record[0].get('Runs', '').split(',')
                # Extract the SRA ID from the runs list
                runs_sra_ids = [re.search(r'SRR\d+', run).group() for run in runs]
                sra_ids.extend(runs_sra_ids)

            # Store the SRA IDs for the current GSM ID
            gsm_to_sra[gsm_id] = sra_ids

        except Exception as e:
            print(f"Error fetching SRA for {gsm_id}: {e}")
            gsm_to_sra[gsm_id] = []

    return gsm_to_sra

# Download SRA using prefetch from SRA Toolkit
def download_sra(sra_id, output_dir="sra_files"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        subprocess.run(["prefetch", sra_id, "--output-directory", output_dir], check=True)
        print(f"Downloaded: {sra_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {sra_id}: {e}")

# Convert SRA to FASTQ using fastq-dump from SRA Toolkit
def sra_to_fastq(sra_id, input_dir="sra_files", output_dir="fastq_files"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # sra_path = os.path.join(input_dir, f"{sra_id}.sra")
    sra_path = os.path.join(input_dir, sra_id, f"{sra_id}.sra")

    
    if os.path.exists(sra_path):
        try:
            subprocess.run(["fastq-dump", "--split-files", "--outdir", output_dir, sra_path], check=True)
            print(f"Converted {sra_id} to FASTQ.")
        except subprocess.CalledProcessError as e:
            print(f"Error converting {sra_id}: {e}")
    else:
        print(f"SRA file not found for {sra_id} at {sra_path}")

# Run both download_sra and sra_to_fastq functions for each SRA ID in gsm_to_sra dict
def process_sra_from_gsm(gsm_to_sra, sra_output_dir="sra_files", fastq_output_dir="fastq_files"):
    for gsm_id, sra_ids in gsm_to_sra.items():
        print(f"Processing GSM: {gsm_id}")
        for sra_id in sra_ids:
            # Download SRA
            download_sra(sra_id, output_dir=sra_output_dir)
            
            # Convert SRA to FASTQ
            sra_to_fastq(sra_id, input_dir=sra_output_dir, output_dir=fastq_output_dir)

# Get GSM IDs for GSE accession ID
gse_accession = "GSE264108"
gsm_ids = get_gsm_ids_from_gse(gse_accession)
# print returned GSM IDs
print("GSM IDs:", gsm_ids)

# Get SRA IDs for each GSM ID
gsm_to_sra = get_sra_ids_from_gsm(gsm_ids)
print("SRA IDs for each GSM ID:", gsm_to_sra)

# Download all SRA files and convert to FASTQ
process_sra_from_gsm(gsm_to_sra)