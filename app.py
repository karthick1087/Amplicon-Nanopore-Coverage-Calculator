import streamlit as st
import os
import subprocess
import pandas as pd
from multiprocessing import Pool
import tempfile
import logging

logging.basicConfig(level=logging.INFO)

st.markdown("""
<style>
    .stApp { background-color: #ffffff; }
    .stMarkdown h1, .stMarkdown h4, .stMarkdown p, .stMarkdown span { color: #000000 !important; }
    .stFileUploader section { border: 2px dashed #000000; }
    .stButton button { background-color: #000000 !important; color: white !important; }
</style>
""", unsafe_allow_html=True)

def main():
    st.markdown("<h1>Amplicon Coverage Analyzer</h1>", unsafe_allow_html=True)
    st.markdown("<h4>Upload FASTQ files and a reference FASTA to compute amplicon coverage.</h4>", unsafe_allow_html=True)

    if 'processed' not in st.session_state:
        st.session_state.processed = False
    if 'output_data' not in st.session_state:
        st.session_state.output_data = None
    if 'dataframe' not in st.session_state:
        st.session_state.dataframe = None

    with st.sidebar:
        st.header("⚙️ Settings")
        reference_file = st.file_uploader("Upload Reference FASTA", type=['fasta', 'fa'])
        fastq_files = st.file_uploader("Upload FASTQ Files", type=['fastq', 'fastq.gz'], accept_multiple_files=True)
        process_btn = st.button("🚀 Start Analysis", use_container_width=True)

    if process_btn:
        if not reference_file:
            st.error("❌ Please upload a reference FASTA file.")
            return
        if not fastq_files:
            st.error("❌ Please upload at least one FASTQ file.")
            return

        with st.spinner("🔬 Processing files..."):
            try:
                output_data, df = process_all_reads_with_progress(reference_file, fastq_files)
                st.session_state.output_data = output_data
                st.session_state.dataframe = df
                st.session_state.processed = True
                st.success("✅ Analysis completed!")
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                logging.error(f"Unexpected error: {str(e)}")

    if st.session_state.processed:
        st.subheader("📊 Coverage Summary Table")
        st.dataframe(st.session_state.dataframe.head(20), use_container_width=True)

        st.download_button(
            label="📥 Download Full Excel Report",
            data=st.session_state.output_data,
            file_name="final_coverage_report.xls",
            mime="application/vnd.ms-excel",
            use_container_width=True
        )

def process_all_reads_with_progress(reference_file, fastq_files):
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_path = os.path.join(tmpdir, "reference.fasta")
        with open(ref_path, "wb") as f:
            f.write(reference_file.getbuffer())

        reads_dir = os.path.join(tmpdir, "reads")
        os.makedirs(reads_dir, exist_ok=True)
        for fq in fastq_files:
            fq_path = os.path.join(reads_dir, fq.name)
            with open(fq_path, "wb") as f:
                f.write(fq.getbuffer())

        output_dir = os.path.join(tmpdir, "results")
        os.makedirs(output_dir, exist_ok=True)
        final_output = os.path.join(output_dir, "final_coverage_report.xls")

        barcode_files = [f for f in os.listdir(reads_dir) if f.endswith((".fastq", ".fastq.gz"))]

        # Progress bar setup
        progress = st.progress(0)
        total = len(barcode_files)
        completed = 0

        def process_and_update(barcode_file):
            nonlocal completed
            try:
                process_file(barcode_file, reads_dir, output_dir, ref_path)
            finally:
                completed += 1
                progress.progress(completed / total)

        for bf in barcode_files:
            process_and_update(bf)

        # Merge results
        merged_df = merge_coverage_matrix(output_dir, final_output)

        with open(final_output, "rb") as f:
            output_data = f.read()

        return output_data, merged_df

def process_file(barcode_file, reads_dir, output_dir, reference):
    barcode_path = os.path.join(reads_dir, barcode_file)
    barcode_name = os.path.splitext(os.path.splitext(barcode_file)[0])[0]
    output_prefix = os.path.join(output_dir, barcode_name)

    sorted_bam = align_reads(barcode_path, output_prefix, reference)
    calculate_summary_coverage(sorted_bam, output_prefix)

def align_reads(input_file, output_prefix, reference):
    bam_file = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}_sorted.bam"

    subprocess.run(
        f"minimap2 -ax map-ont {reference} {input_file} | samtools view -b -F 2308 - > {bam_file}",
        shell=True, check=True
    )
    subprocess.run(f"samtools sort -m 1G -o {sorted_bam} {bam_file}", shell=True, check=True)
    subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
    os.remove(bam_file)
    return sorted_bam

def calculate_summary_coverage(bam_file, output_prefix):
    coverage_tsv = f"{output_prefix}_coverage.tsv"
    coverage_csv = f"{output_prefix}_coverage.csv"

    subprocess.run(f"samtools coverage {bam_file} > {coverage_tsv}", shell=True, check=True)

    with open(coverage_tsv, 'r') as f:
        for line in f:
            if line.startswith('#rname'):
                header = line.strip('#').strip().split('\t')
                break

    df = pd.read_csv(coverage_tsv, sep='\t', comment='#', header=None)
    df.columns = header if header else ['rname', 'startpos', 'endpos', 'numreads']

    if all(col in df.columns for col in ['#rname', 'startpos', 'endpos', 'numreads']):
        df = df[['#rname', 'startpos', 'endpos', 'numreads']]
    elif all(col in df.columns for col in ['rname', 'startpos', 'endpos', 'numreads']):
        df.rename(columns={'rname': '#rname'}, inplace=True)
        df = df[['#rname', 'startpos', 'endpos', 'numreads']]
    else:
        raise ValueError("Expected coverage output columns not found.")

    df.to_csv(coverage_csv, index=False)
    os.remove(coverage_tsv)

def merge_coverage_matrix(output_dir, final_output):
    combined_entries = {}
    coverage_files = [f for f in os.listdir(output_dir) if f.endswith("_coverage.csv")]

    for f in coverage_files:
        barcode = os.path.splitext(f)[0].split('_')[0]
        df = pd.read_csv(os.path.join(output_dir, f))

        for _, row in df.iterrows():
            key = (row['startpos'], row['endpos'])
            if key not in combined_entries:
                combined_entries[key] = {
                    '#rname': row['#rname'],
                    'startpos': row['startpos'],
                    'endpos': row['endpos'],
                    f'numreads_{barcode}': row['numreads']
                }
            else:
                if not row['#rname'].endswith('_R'):
                    combined_entries[key]['#rname'] = row['#rname']
                combined_entries[key][f'numreads_{barcode}'] = row['numreads']

    merged_df = pd.DataFrame.from_dict(combined_entries, orient='index').reset_index(drop=True)
    merged_df = merged_df.sort_values(by=['startpos'])
    merged_df.to_excel(final_output, index=False)
    return merged_df

if __name__ == "__main__":
    main()

