# Amplicon Coverage Analyzer

A web-based tool for analyzing sequencing coverage from nanopore amplicon data.

## Author
Dr. Karthick Vasudevan  
Institute of Bioinformatics  
karthick@ibioinformatics.org  

## Features
- Upload reference FASTA and amplicon FASTQ files
- Real-time alignment using Minimap2
- Coverage calculation using Samtools
- Excel report generation
- Visual progress bar and on-screen preview

## Requirements
Install these tools before running the app:
```bash
conda install -c bioconda minimap2 samtools
```

## Run Locally
```bash
streamlit run app.py
```

## Deploy on Streamlit Cloud
You may need to use a custom Docker image that includes `minimap2` and `samtools`.
