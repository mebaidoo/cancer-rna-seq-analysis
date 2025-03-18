import requests
import os

# Specify the file ID
file_id = "ccd692f2-ae9d-48ad-8d6d-70385f97420a"

# Get the current working directory
current_dir = os.getcwd()

# Construct the URL for the GDC API
url = f"https://api.gdc.cancer.gov/data/{file_id}"

# Download the file
print(f"Downloading: {file_id}")
response = requests.get(url, stream=True)  # stream to download the file in chunks, useful for large files

# Save the file in the current working directory
output_file = os.path.join(current_dir, f"{file_id}.tar.gz")
with open(output_file, "wb") as output:
    for chunk in response.iter_content(chunk_size=1024):
        if chunk:
            output.write(chunk)
print(f"Saved to {output_file}")

