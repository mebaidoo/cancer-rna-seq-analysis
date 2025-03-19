import requests
import pandas as pd
import io

# TCGA file ID (replace with your actual ID)
file_id = "ccd692f2-ae9d-48ad-8d6d-70385f97420a"

# Step 1: Download the raw data
url = f"https://api.gdc.cancer.gov/data/{file_id}"
response = requests.get(url)

if response.status_code == 200:
    print(f"Downloading TCGA data for: {file_id}")
    
    # Step 2: Read data into a DataFrame (assuming it's tabular data)
    data = pd.read_csv(io.BytesIO(response.content), skiprows=1, sep="\t")
    
    # Step 3: Save to an Excel file
    output_file = "tcga_data.xlsx"
    data.to_excel(output_file, index=False)

    print(f"Data saved as: {output_file}")

else:
    print(f"Error: Unable to download file (status code: {response.status_code})")