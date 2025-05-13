#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
import argparse
import requests
import re
from Bio import SeqIO


def sbatch(job_name, command, index):
    print(f"Running job: {index}_{job_name}")
    subprocess.run(command, shell=True, check=True)

def filter_gtf(gtf_file, gene_id, output_dir):
    # Filter the GTF for the given gene_id.
    command = f"grep 'gene_id \"{gene_id}\"' {gtf_file} > {output_dir}/gene_{gene_id}.gtf"
    return command
    
def get_mane_transcript_id(output_dir, gene_id):
    """
    Read the filtered GTF file and extract the transcript ID corresponding to the MANE transcript.
    If a MANE_Select tag is present in the attributes, use those rows to extract the transcript_id.
    Otherwise, fall back to using all rows (or adjust as needed).
    Returns the MANE transcript ID as a string or None if not found.
    """
    gtf_file = f"{output_dir}/gene_{gene_id}.gtf"
    columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    
    try:
        gtf = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=columns)
    except pd.errors.EmptyDataError:
        return None
    
    # If a MANE_Select tag is present, filter the rows
    if gtf['attributes'].str.contains("MANE_Select", case=False, na=False).any():
        mane_rows = gtf[gtf['attributes'].str.contains("MANE_Select", case=False, na=False)]
    else:
        mane_rows = gtf
    
    transcript_ids = set()
    
    # Iterate over the attributes to extract transcript_id values
    for attr in mane_rows['attributes']:
        # Split the attribute field by semicolons
        fields = attr.split(';')
        for field in fields:
            field = field.strip()
            # Look for the transcript_id key
            if field.startswith("transcript_id"):
                # Expecting a format like: transcript_id "ENST00000367770"
                parts = field.split(" ")
                if len(parts) >= 2:
                    tid = parts[1].replace('"', '').strip()
                    transcript_ids.add(tid)
    
    if transcript_ids:
        # If multiple IDs are found, return the first in sorted order (or adjust selection criteria as needed)
        return sorted(transcript_ids)[0]
    else:
        return None

def fetch_ensembl_transcript_sequence(ensembl_transcript_id):
    """
    Retrieve the cDNA (mRNA) sequence for the given Ensembl transcript ID using the Ensembl REST API.
    
    Parameters:
        ensembl_transcript_id (str): The Ensembl transcript ID (e.g., "ENST00000367770").
        
    Returns:
        str: The cDNA sequence in plain text if successful, None otherwise.
    """
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_transcript_id}?type=cdna"
    headers = {"Content-Type": "text/plain"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        print(f"Error retrieving sequence for {ensembl_transcript_id}: {response.text}")
        return None
    return response.text

def get_strand_from_gtf(gene_id, output_dir):
    """
    Retrieve the strand information for the gene from the filtered GTF.
    (Assumes that a line with the 'gene' feature exists.)
    """
    gtf_file = f"{output_dir}/gene_{gene_id}.gtf"
    command = f"grep -w 'gene' {gtf_file} | awk '$3 == \"gene\" {{print $7}}'"
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        strand = result.stdout.strip()
        if strand not in ['+', '-']:
            raise ValueError("Invalid strand information retrieved.")
        print(f"Strand for gene {gene_id}: {strand}")
        return strand
    except subprocess.CalledProcessError as e:
        print(f"Error retrieving strand for gene {gene_id}: {e}")
        return None
    

def write_primer3_input(gene_id, mRNA, output_file, TM=60.0, min_primer=18, max_primer=20, repeat_lib="/user/gent/446/vsc44685/ScratchVO_dir/NSQ2K_126_test/primer_design/rep_lib/humrep_and_simple.txt"):
    minTM = TM - 4
    maxTM = TM + 4
    with open(output_file, 'w') as f:
        f.write(f"SEQUENCE_ID={gene_id}\n")
        f.write(f"SEQUENCE_TEMPLATE={mRNA}\n")
        f.write("PRIMER_MIN_THREE_PRIME_DISTANCE=30\n")
        # f.write(f"PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{max_product_size}\n")
        f.write("PRIMER_TASK=pick_primer_list\n")
        f.write("PRIMER_NUM_RETURN=10\n")
        f.write("PRIMER_OPT_SIZE=20\n")
        f.write("PRIMER_PICK_LEFT_PRIMER=1\n")
        f.write("PRIMER_PICK_RIGHT_PRIMER=0\n")
        f.write(f"PRIMER_MIN_SIZE={min_primer}\n")
        f.write(f"PRIMER_MAX_SIZE={max_primer}\n")
        f.write(f"PRIMER_OPT_TM={TM}\n")
        f.write(f"PRIMER_MIN_TM={minTM}\n")
        f.write(f"PRIMER_MAX_TM={maxTM}\n")
        f.write("PRIMER_OPT_GC_PERCENT=50.0\n")
        f.write("PRIMER_MIN_GC=40.0\n")
        f.write("PRIMER_MAX_GC=60.0\n")
        f.write(f"PRIMER_MISPRIMING_LIBRARY={repeat_lib}\n") 
        f.write(f"PRIMER_MAX_LIBRARY_MISPRIMING=14.00\n") #The maximum allowed weighted similarity with any sequence in the mispriming library
        f.write("PRIMER_EXPLAIN_FLAG=1\n")
        f.write("=\n")
    return output_file


def Run_Primer3(input_file, output_dir):
    command = f"singularity run /scratch/gent/vo/000/gvo00027/singularity_containers/primer3_v2.5.0.sif --output={output_dir} {input_file}"
    return command


def parse_primer3_output(file_path):
    gene_id = None
    primers = {}   # to store primer sequences keyed by side and index (e.g. 'LEFT_0')
    positions = {} # to store corresponding positions
    TM={}
    GC={}

    # Regex patterns to capture primer sequence and position lines:
    # e.g., "PRIMER_LEFT_0_SEQUENCE=ATCG..."
    seq_pattern = re.compile(r"PRIMER_LEFT_([0-9]+)_SEQUENCE=(.+)")
    # e.g., "PRIMER_LEFT_0=123,20" (position and length)
    pos_pattern = re.compile(r"PRIMER_LEFT_([0-9]+)=(.+)")
    TM_pattern = re.compile(r"PRIMER_LEFT_([0-9])+_TM=(.+)")
    GC_pattern = re.compile(r"PRIMER_LEFT_([0-9])+_GC_PERCENT=(.+)")

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Get the gene ID
            if line.startswith("SEQUENCE_ID"):
                gene_id = line.split('=')[1]
            else:
                # Check for primer sequence lines
                seq_match = seq_pattern.match(line)
                if seq_match:
                    index = seq_match.group(1)
                    sequence = seq_match.group(2)
                    key = f"primer_{index}"
                    primers[key] = sequence
                    continue  # move to next line

                # Check for primer position lines (only if needed)
                pos_match = pos_pattern.match(line)
                if pos_match:
                    # To avoid capturing the sequence line again, only record if key doesn't include 'SEQUENCE'
                    key = f"primer_{pos_match.group(1)}"
                    # Ensure that this line isn’t a sequence line (we already matched those)
                    if key not in positions:
                        pos_and_size = pos_match.group(2)
                        start = int(pos_and_size.split(',')[0])
                        primer_size = int(pos_and_size.split(',')[1])
                        positions[key] = (start, primer_size)
                    continue
                TM_match = TM_pattern.match(line)
                if TM_match:
                    key = f"primer_{TM_match.group(1)}"
                    if key not in TM:
                        TM[key] = TM_match.group(2)
                    continue
                
                GC_match = GC_pattern.match(line)
                if GC_match:
                    key = f"primer_{GC_match.group(1)}"
                    # Check if the key exists in the primers dictionary
                    if key not in GC:
                        GC[key] = GC_match.group(2)
                    continue
    # Convert the primers dictionary to a list (or you can return the dict if you want the keys)
    # primer_list = list(primers.values())
    return gene_id, primers, positions, TM, GC
                        
                        
    
    # Convert the primers dictionary to a list (or you can return the dict if you want the keys)
    # primer_list = list(primers.values())
    return gene_id, primers, positions

def write_primers_to_fasta(gene_id, primers, output_path):
    """
    Writes primer sequences to a FASTA file.

    Parameters:
      gene_id (str): The gene identifier.
      primers (dict): A dictionary with primer indices as keys and sequences as values.
      output_path (str): The file path where the FASTA file will be written.

    Each FASTA header is in the format: >gene_id_LEFT_index
    """
    with open(output_path, 'a') as f:
        # Sort keys by numerical value if keys are numbers stored as strings.
        for index in sorted(primers):
            sequence = primers[index]
            header = f">{gene_id}_LEFT_{index}"
            f.write(header + "\n")
            # Optionally, wrap the sequence to 80 characters per line.
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")


def parse_transcriptome(transcriptome):
    transcript_info = {}
    with open(transcriptome) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            transcript_info[record.id] = record.description
    return transcript_info
#run alignment with bowtie it takes an index prefix and a fasta file with primers                
def run_bowtie2_command(primers_fasta,index_prefix, output_file="bowtie2_results.sam"):
    input_file = "query.fasta"
    command = f"ml purge; ml Bowtie2/2.5.4-GCC-13.2.0; bowtie2 -f -x {index_prefix} -U {primers_fasta} -S {output_file} -a"
    return command

def run_make_alignment_summary(transcriptome, sam_file, out_file):
    command = f"ml purge; ml Biopython/1.84-foss-2023b; ml Pysam/0.22.0-GCC-13.2.0; python3 make_alignment_summary.py --transcriptome {transcriptome} --sam_file {sam_file} --out_file {out_file}"
    return command

#creates a dataframe with primer3 ouputs and calulates the distance to the end of the transcript
def make_primer_dataframe(gene, primer_seqs, primer_positions, TM, GC, ampl_size, file_path):
    """
    Creates a pandas DataFrame from gene, primer sequences, and primer positions.

    Parameters:
        gene (str): The gene identifier.
        primer_seqs (dict): A dictionary with primer indices as keys and sequences as values.
        primer_positions (dict): A dictionary with primer indices as keys and positions as values.

    Returns:
        pd.DataFrame: A DataFrame containing the gene, primer index, sequence, and position.
    """
    data = []
    for index in primer_seqs:
        sequence = primer_seqs[index]
        position = primer_positions.get(index, None)
        start = int(position[0])
        primer_size = int(position[1])
        distance = int(ampl_size - start)
        TempMelt = TM.get(index, None)
        GCpercent = GC.get(index, None)
        data.append({'Gene': gene, 'Primer_Index': index, 'Sequence': sequence, 'Primer_size': primer_size, 'TempMelt': TempMelt, 'GC_percent':GCpercent, 'Distance_to_End': distance})
    df_new = pd.DataFrame(data)
    
    # Check if the CSV file already exists
    if os.path.exists(file_path):
        # Read the existing CSV file
        df_existing = pd.read_csv(file_path)
        # Append the new records to the existing DataFrame
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        # Write the combined DataFrame back to the CSV file
        df_combined.to_csv(file_path, index=False)
        return df_combined
    else:
        # File doesn't exist, so create a new CSV file with the new DataFrame
        df_new.to_csv(file_path, index=False)
        return df_new

    

def main(out_dir, gtf_file, gene_id, ampl_size, TM, lib, min_primer_size ,max_primer_size):
    output_dir = f"{out_dir}/{gene_id}_out"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filter_gtf_command = filter_gtf(gtf_file, gene_id, output_dir)
    sbatch("filter_gtf", filter_gtf_command, 1)
    #Step 1: Get the MANE transcript ID
    print("Running job: 2_Fetching_MANE_transcript_sequence")
    mane_tid = get_mane_transcript_id(output_dir, gene_id)
    if mane_tid:
        mRNA = fetch_ensembl_transcript_sequence(mane_tid)
        if mRNA:
            print(f"{mane_tid}: {mRNA}")
    else:
        print("No Transcript ID found for gene", gene_id)

    #step 2: Get the strand information
    print("Running job: 3_Strand_from_gtf")
    strand = get_strand_from_gtf(gene_id, output_dir)
    strand_label = 'pos' if strand == '+' else 'neg'
    #Step 3: Write the Primer3 input file
    print("Running job: 4_make_primer3_input")
    sequence= mRNA[-ampl_size:]
    primer3_input_file = f"{output_dir}/{gene_id}_primer3_input_{strand_label}.txt"
    primer3_input_path = write_primer3_input(gene_id, sequence, primer3_input_file, TM, min_primer=min_primer_size, max_primer=max_primer_size, repeat_lib=lib)
    print(f"Primer3 input file written to: {primer3_input_path}")
    primer3_output_file = f"{output_dir}/{gene_id}_primer3_out.txt"
    PrimerFasta_output_file = f"{out_dir}/primers.fa"
    #STep 4: Run Primer3
    primer3_command = Run_Primer3(primer3_input_file, primer3_output_file)
    sbatch("Run_Primer3", primer3_command, 5)
    gene, primer_seqs, primer_positions, TM, GC = parse_primer3_output(primer3_output_file)
    write_primers_to_fasta(gene, primer_seqs, PrimerFasta_output_file)
    print("Running job: 6_make_summary_files")
    summary_csv_file = f'{out_dir}/Primer3_summary.csv'
    df = make_primer_dataframe(gene, primer_seqs, primer_positions, TM, GC, ampl_size, summary_csv_file)
    print("All steps finished")
    return df