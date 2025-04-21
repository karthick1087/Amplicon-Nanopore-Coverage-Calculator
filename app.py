import streamlit as st
import os
import subprocess
import pandas as pd
import tempfile
import zipfile
import logging
import plotly.express as px
import re

logging.basicConfig(level=logging.INFO)
st.set_page_config(page_title="Amplicon Coverage Analyzer", layout="wide")

# App styling
st.markdown("""
<style>
    .stApp { background-color: #607cbd; }
    .stMarkdown h1, .stMarkdown h4, .stMarkdown p, .stMarkdown span { color: #000000 !important; }
    .stFileUploader section { border: 2px dashed #000000; }
    .stButton button { background-color: #000000 !important; color: white !important; }
</style>
""", unsafe_allow_html=True)

def sanitize(text):
    """Remove emojis and surrogate characters from rname"""
    return re.sub(r'[^\w\-.]', '_', text)

def main():
    st.markdown("<h1>Nanopore Amplicon Coverage Analyzer</h1>", unsafe_allow_html=True)
    st.markdown("<h4>Upload FASTQ files or zipped barcode folders and a reference FASTA.</h4>", unsafe_allow_html=True)

    # Acknowledgment
    st.markdown("""
    <div style="text-align: center; padding: 10px; margin-bottom: 20px;">
        <p style="color: black; font-size: 15px;">
            🧬 <b>Developed by</b> Dr. Karthick Vasudevan<br>
            <i>Institute of Bioinformatics</i><br>
            📧 <a href="mailto:karthick@ibioinformatics.org" style="color: black; text-decoration: none;">
                karthick@ibioinformatics.org
            </a>
        </p>
    </div>
    """, unsafe_allow_html=True)

    reference_file = st.file_uploader("🧬 Upload Reference FASTA", type=("fasta", "fa"))
    input_mode = st.radio("Select Input Type", ["Upload .fastq.gz files directly", "Upload zipped barcode folders (.zip)"])
    fastq_files = zip_files = None

    if input_mode == "Upload .fastq.gz files directly":
        fastq_files = st.file_uploader("📁 Upload FASTQ.GZ Files (1 per barcode)", type="gz", accept_multiple_files=True)
    else:
        zip_files = st.file_uploader("📁 Upload ZIP files (1 per barcode folder)", type="zip", accept_multiple_files=True)

    if st.button("🚀 Start Analysis"):
        if not reference_file or (not fastq_files and not zip_files):
            st.error("❌ Please upload reference and input files.")
            return

        with st.spinner("🔬 Processing..."):
            try:
                output_data, df = process_inputs(reference_file, fastq_files, zip_files)
                st.session_state.output_data = output_data
                st.session_state.dataframe = df
                st.session_state.processed = True
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                logging.error(str(e))

    if st.session_state.get("processed"):
        df = st.session_state.dataframe
        st.markdown("<h3>📊 Coverage Summary Table</h3>", unsafe_allow_html=True)

        threshold = st.slider("Minimum Read Count to Display", 0, 50000, 0, step=1000)
        filtered_df = df.copy()
        for col in filtered_df.columns:
            if col.startswith("numreads_"):
                filtered_df = filtered_df[filtered_df[col] >= threshold]

        if not filtered_df.empty:
            st.dataframe(filtered_df, use_container_width=True)
        else:
            st.warning("⚠️ No data above selected threshold.")

        st.subheader("📈 Coverage Plot")
        melted = df.melt(id_vars=["#rname"],
                         value_vars=[col for col in df.columns if col.startswith("numreads_")],
                         var_name="Barcode", value_name="Read Count")
        melted["Barcode"] = melted["Barcode"].str.replace("numreads_", "")
        fig = px.bar(melted, x="#rname", y="Read Count", color="Barcode", barmode="group")
        st.plotly_chart(fig, use_container_width=True)

        st.download_button("📥 Download Excel Report", st.session_state.output_data,
                           file_name="coverage_report.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

def process_inputs(reference_file, fastq_files, zip_files):
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_path = os.path.join(tmpdir, "reference.fasta")
        with open(ref_path, "wb") as f:
            f.write(reference_file.getbuffer())

        reads_dir = os.path.join(tmpdir, "reads")
        os.makedirs(reads_dir, exist_ok=True)

        # Option A
        if fastq_files:
            for fq in fastq_files:
                with open(os.path.join(reads_dir, fq.name), "wb") as f:
                    f.write(fq.getbuffer())

        # Option B
        if zip_files:
            for zip_file in zip_files:
                barcode = os.path.splitext(zip_file.name)[0]
                unzip_path = os.path.join(reads_dir, barcode)
                os.makedirs(unzip_path, exist_ok=True)
                with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                    zip_ref.extractall(unzip_path)

                output_fq = os.path.join(reads_dir, f"{barcode}.fastq.gz")
                with open(output_fq, "wb") as out_f:
                    for root, _, files in os.walk(unzip_path):
                        for f in sorted(files):
                            if f.endswith(".fastq.gz"):
                                with open(os.path.join(root, f), "rb") as in_f:
                                    out_f.write(in_f.read())

        output_dir = os.path.join(tmpdir, "results")
        os.makedirs(output_dir, exist_ok=True)
        final_excel = os.path.join(output_dir, "coverage_report.xlsx")

        input_files = [f for f in os.listdir(reads_dir) if f.endswith(".fastq.gz")]
        progress = st.progress(0)
        status = st.empty()

        for idx, fq_file in enumerate(input_files, 1):
            sample = os.path.splitext(os.path.splitext(fq_file)[0])[0]
            fq_path = os.path.join(reads_dir, fq_file)
            output_prefix = os.path.join(output_dir, sample)

            bam = f"{output_prefix}.bam"
            sorted_bam = f"{output_prefix}_sorted.bam"

            status.info(f"🔗 Aligning: {fq_file}")
            subprocess.run(f"minimap2 -ax map-ont {ref_path} {fq_path} | samtools view -b -F 2308 - > {bam}", shell=True, check=True)
            subprocess.run(f"samtools sort -m 1G -o {sorted_bam} {bam}", shell=True, check=True)
            subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
            os.remove(bam)

            calculate_summary_coverage(sorted_bam, output_prefix)
            progress.progress(idx / len(input_files))

        merged_df = merge_coverage_matrix(output_dir, final_excel)
        with open(final_excel, "rb") as f:
            output_data = f.read()
        return output_data, merged_df

def calculate_summary_coverage(bam_file, output_prefix):
    coverage_tsv = f"{output_prefix}_coverage.tsv"
    coverage_csv = f"{output_prefix}_coverage.csv"
    subprocess.run(f"LC_ALL=C samtools coverage {bam_file} > {coverage_tsv}", shell=True, check=True, executable="/bin/bash")

    combined_entries = {}
    with open(coverage_tsv, 'r', encoding='utf-8', errors='ignore') as infile:
        for line in infile:
            if line.startswith('#rname'): continue
            fields = line.strip().split('\t')
            if len(fields) < 4: continue
            rname = sanitize(fields[0])
            startpos = int(fields[1])
            endpos = int(fields[2])
            numreads = int(float(fields[3]))
            key = (startpos, endpos)
            if key in combined_entries:
                combined_entries[key]['numreads'] += numreads
                if not rname.endswith('_R') and not rname.endswith('_R1'):
                    combined_entries[key]['rname'] = rname
            else:
                combined_entries[key] = {
                    'rname': rname,
                    'startpos': startpos,
                    'endpos': endpos,
                    'numreads': numreads
                }

    with open(coverage_csv, 'w', encoding='utf-8') as outfile:
        outfile.write('#rname,startpos,endpos,numreads\n')
        for entry in combined_entries.values():
            outfile.write(f"{entry['rname']},{entry['startpos']},{entry['endpos']},{entry['numreads']}\n")
    os.remove(coverage_tsv)

def merge_coverage_matrix(output_dir, final_output):
    coverage_files = sorted([f for f in os.listdir(output_dir) if f.endswith("_coverage.csv")])
    combined = {}
    barcodes = []

    for f in coverage_files:
        barcode = os.path.splitext(f)[0].split("_")[0]
        barcodes.append(barcode)
        df = pd.read_csv(os.path.join(output_dir, f))
        for _, row in df.iterrows():
            key = (row['#rname'], row['startpos'], row['endpos'])
            if key not in combined:
                combined[key] = {
                    '#rname': row['#rname'],
                    'startpos': row['startpos'],
                    'endpos': row['endpos'],
                    f'numreads_{barcode}': row['numreads']
                }
            else:
                combined[key][f'numreads_{barcode}'] = row['numreads']

    merged_df = pd.DataFrame.from_dict(combined, orient='index').reset_index(drop=True)
    merged_df = merged_df.sort_values(by=['startpos'])
    cols = ['#rname', 'startpos', 'endpos'] + [f'numreads_{b}' for b in barcodes]
    for col in cols:
        if col not in merged_df.columns:
            merged_df[col] = 0
    merged_df = merged_df[cols]
    merged_df.to_excel(final_output, index=False)
    return merged_df

if __name__ == "__main__":
    main()
