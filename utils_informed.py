#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
import argparse
import re

def sbatch(job_name, command, index):
    print(f"Running job: {index}_{job_name}")
    subprocess.run(command, shell=True, check=True)

def get_coverage(bam_file, output_dir):
    # Generate a BED file with coverage using bedtools genomecov.
    coverage_bed = os.path.join(output_dir, "coverage.bed")
    if os.path.exists(coverage_bed):
        print(f"Coverage file already exists: {coverage_bed}")
        return None
    else: 
        command = ["ml purge; ml BEDTools/2.31.0-GCC-12.3.0;",
        "bedtools genomecov -ibam", bam_file, f"-bg > {coverage_bed}"]
    return " ".join(command)

# def get_coverage(bam_file, output_dir):
#     # Generate a BED file with coverage using bedtools genomecov.
#     command = ["ml purge; ml BEDTools/2.31.0-GCC-12.3.0;", 
#                "bedtools genomecov -ibam", bam_file, f"-bg > {output_dir}/coverage.bed"]
#     return " ".join(command)

def filter_gtf(gtf_file, gene_id, output_dir):
    """Filter the GTF for the given gene_id"""
    gtf_filtered = os.path.join(output_dir, f"gene_{gene_id}.gtf")
    if os.path.exists(gtf_filtered):
        print(f"Filtered GTF file already exists: {gtf_filtered}")
        return None
    else:
        command = f"grep 'gene_id \"{gene_id}\"' {gtf_file} > {output_dir}/gene_{gene_id}.gtf"
    return command

def merge_exons(output_dir, gene_id):
    """
    Extract only exon features from the filtered GTF and merge them.
    The large merge distance (-d 1000000) is used so that all exons for the gene
    are merged into one continuous interval for primer design.
    """
    input_file = f"{output_dir}/gene_{gene_id}.gtf"
    output_file = f"{output_dir}/gene_{gene_id}_merged_exons.bed"
    if os.path.exists(output_file):
        print(f"Filtered Exon-merged Bed file already exists: {output_file}")
        return None
    else:
        command = f"ml purge; ml BEDTools/2.31.0-GCC-12.3.0; grep -w 'exon' {input_file} | bedtools sort -i - | bedtools merge -d 1000000 -i - > {output_file}"
    return command

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))


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

def get_gene_coverage(output_dir, gene_id):
    """
    Use the merged exons file to intersect with the genome coverage file.
    This avoids interpreting intronic gaps as zero coverage.
    """
    merged_exons_file = f"{output_dir}/gene_{gene_id}_merged_exons.bed"
    command = f"ml purge; ml BEDTools/2.31.0-GCC-12.3.0; bedtools intersect -a {output_dir}/coverage.bed -b {merged_exons_file} > {output_dir}/gene_{gene_id}_coverage.bed"
    return command

def get_peak_and_low_coverage_prct(window_size, output_dir, gene_id, strand, coverage_percentage):
    """
    Determines the highest coverage region and finds the first coordinate
    after the peak where the coverage falls below a specified percentage
    of the highest coverage. This is done in a strand-aware manner.
    
    Parameters:
      window_size (int): Number of bases to expand around the peak for extraction.
      output_dir (str): Directory where the gene's coverage BED file is located.
      gene_id (str): Identifier for the gene.
      strand (str): '+' or '-' indicating the transcriptional strand.
      coverage_percentage (float, optional): The percentage (of the highest coverage)
                                               to use as the threshold for identifying
                                               a low coverage region.
    
    Returns:
      tuple: (chrom, start_expanded, end_expanded, peak_coord, low_cov_coord)
             where low_cov_coord is the coordinate at which the coverage first 
             drops below the relative threshold.
    """
    import pandas as pd
    
    bed_file = f"{output_dir}/gene_{gene_id}_coverage.bed"
    
    # Load the BED file (assumed columns: chrom, start, end, coverage)
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'coverage'])
    bed_df = bed_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
    
    # Identify the interval with the maximum coverage and compute its midpoint.
    max_coverage_row = bed_df.loc[bed_df['coverage'].idxmax()]
    chrom = max_coverage_row['chrom']
    start = max_coverage_row['start']
    end = max_coverage_row['end']
    peak_coord = int((start + end) / 2)
    
    # Expand the region around the peak.
    start_expanded = max(0, peak_coord - window_size)
    end_expanded = peak_coord + window_size
    
    print(f"Highest coverage region: {chrom}:{start}-{end} on strand {strand} with coverage {max_coverage_row['coverage']}")
    print(f"Extracting expanded region: {chrom}:{start_expanded}-{end_expanded}")
    
    # Calculate the relative threshold based on the peak's coverage.
    max_cov = max_coverage_row['coverage']
    relative_threshold = max_cov * (coverage_percentage / 100.0)
    
    # Limit to intervals on the same chromosome.
    bed_chr = bed_df[bed_df['chrom'] == chrom].reset_index(drop=True)
    
    # Find the interval that contains the peak.
    matching = bed_chr[(bed_chr['start'] <= peak_coord) & (bed_chr['end'] > peak_coord)]
    if matching.empty:
        print(f"Peak coordinate {peak_coord} not found in any interval on {chrom}.")
        low_cov_coord = None
    else:
        idx = matching.index[0]
        # For the positive strand, iterate forward from the interval immediately after the peak.
        if strand == '+':
            low_cov_coord = None
            for i in range(idx + 1, len(bed_chr)):
                if bed_chr.loc[i, 'coverage'] < relative_threshold:
                    low_cov_coord = bed_chr.loc[i, 'start']
                    break
            if low_cov_coord is None:
                # If no low coverage interval was found, default to the end of the last interval.
                low_cov_coord = bed_chr.iloc[-1]['end']
        elif strand == '-':
            low_cov_coord = None
            for i in range(idx - 1, -1, -1):
                if bed_chr.loc[i, 'coverage'] < relative_threshold:
                    low_cov_coord = bed_chr.loc[i, 'end']
                    break
            if low_cov_coord is None:
                # If no low coverage interval was found, default to the start of the first interval.
                low_cov_coord = bed_chr.iloc[0]['start']
        else:
            print("Invalid strand specified. Use '+' or '-'.")
            low_cov_coord = None

    print(f"Peak coordinate: {peak_coord}")
    print(f"Low coverage coordinate (coverage < {coverage_percentage}% of {max_cov} i.e. {relative_threshold}) on strand {strand}: {low_cov_coord}")

    return chrom, start_expanded, end_expanded, peak_coord, low_cov_coord



def check_3prime_utr(output_dir, gene_id, region_chrom, region_start, region_end):
    # Check if the region overlaps a three_prime_utr feature.
    gtf = f"{output_dir}/gene_{gene_id}.gtf"
    gtf_columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    gtf = pd.read_csv(gtf, sep="\t", comment='#', header=None, names=gtf_columns)
    gtf['chrom'] = gtf['chrom'].astype(str)
    region_chrom = str(region_chrom)
    if ((gtf['chrom'] == region_chrom) &
        (gtf['start'] <= region_end) &
        (gtf['end'] >= region_start) &
        (gtf['feature'].str.lower() == 'three_prime_utr')).any():
        return True
    else:
        return False

def get_fasta_sequence(chrom, start_expanded, end_expanded, ref_genome, output_dir, gene_id):
    # Fetch the genomic sequence for the given coordinates.
    output_file = f"{output_dir}/{gene_id}_{chrom}_{start_expanded}_{end_expanded}.fa"
    command = f"ml purge; ml SAMtools/1.18-GCC-12.3.0; samtools faidx {ref_genome} {chrom}:{start_expanded}-{end_expanded} > {output_file}"
    return command

def write_primer3_input2(gene_id, fasta_file, output_file, strand, TM, min_primer, max_primer, min_product_size, max_product_size, repeat_lib="/rep_lib/humrep_and_simple.txt"):
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        sequence = "".join(lines[1:]).replace("\n", "")
        if strand == "-":
            sequence = reverse_complement(sequence)
    strand_label = 'pos' if strand == '+' else 'neg'
    output_file_with_strand = f"{output_file}_{strand_label}.txt"
    
    minTM = TM - 3
    maxTM = TM + 3
    
    with open(output_file_with_strand, 'w') as f:
        f.write(f"SEQUENCE_ID={gene_id}\n")
        f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
        f.write(f"PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{max_product_size}\n")
        f.write("PRIMER_TASK=pick_primer_list\n")
        f.write("PRIMER_NUM_RETURN=5\n")
        f.write("PRIMER_OPT_SIZE=20\n")
        f.write(f"PRIMER_MIN_SIZE={min_primer}\n")
        f.write(f"PRIMER_MAX_SIZE={max_primer}\n")
        f.write(f"PRIMER_OPT_TM={TM}\n")
        f.write(f"PRIMER_MIN_TM={minTM}\n")
        f.write(f"PRIMER_MAX_TM={maxTM}\n")
        f.write("PRIMER_OPT_GC_PERCENT=50.0\n")
        f.write("PRIMER_MIN_GC=40.0\n")
        f.write("PRIMER_MAX_GC=60.0\n")
        f.write(f"PRIMER_MISPRIMING_LIBRARY={repeat_lib}\n") 
        f.write(f"PRIMER_MAX_LIBRARY_MISPRIMING=14.00\n")
        f.write("=\n")
    return output_file_with_strand

# def Run_Primer3(input_file, output_dir):
#     command = f"singularity run ./primer3/primer3_v2.5.0.sif --output={output_dir} {input_file}"
#     return command

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
                
def get_exons(output_dir, gene_id):
    """
    Read the filtered GTF file and extract all unique exon intervals.
    Returns a sorted list of (start, end) tuples.
    """
    gtf_file = f"{output_dir}/gene_{gene_id}.gtf"
    columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    try:
        gtf = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=columns)
    except pd.errors.EmptyDataError:
        return []
    # Filter for exon features (case-insensitive) and select only the start and end columns.
    exons = gtf[gtf['feature'].str.lower() == 'exon'][['start', 'end']]
    # Remove duplicate rows.
    exons = exons.drop_duplicates()
    # Return list of (start, end) tuples sorted by start coordinate.
    exons_list = list(zip(exons['start'], exons['end']))
    return sorted(exons_list, key=lambda x: x[0])

def get_mane_exons(output_dir, gene_id):
    """
    Read the filtered GTF file and extract unique exon intervals for the MANE transcript.
    If a MANE_Select tag is present in the attributes, use only those exons.
    Otherwise, fall back to using all exons.
    Returns a sorted list of (start, end) tuples.
    """
    gtf_file = f"{output_dir}/gene_{gene_id}.gtf"
    columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    try:
        gtf = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=columns)
    except pd.errors.EmptyDataError:
        return []
    
    # Check for a MANE tag in the attributes (this may vary depending on your file)
    if gtf['attributes'].str.contains("MANE_Select", case=False, na=False).any():
        exons = gtf[(gtf['feature'].str.lower() == 'exon') &
                    (gtf['attributes'].str.contains("MANE_Select", case=False, na=False))]
    else:
        exons = gtf[gtf['feature'].str.lower() == 'exon']
    
    # Deduplicate and sort
    exons = exons[['start', 'end']].drop_duplicates()
    exons_list = list(zip(exons['start'], exons['end']))
    return sorted(exons_list, key=lambda x: x[0])

def get_mane_UTR(output_dir, gene_id):
    """
    Read the filtered GTF file and extract unique exon intervals for the MANE transcript.
    If a MANE_Select tag is present in the attributes, use only those exons.
    Otherwise, fall back to using all exons.
    Returns a sorted list of (start, end) tuples.
    """
    gtf_file = f"{output_dir}/gene_{gene_id}.gtf"
    columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    try:
        gtf = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=columns)
    except pd.errors.EmptyDataError:
        return []
    
    # Check for a MANE tag in the attributes (this may vary depending on your file)
    if gtf['attributes'].str.contains("MANE_Select", case=False, na=False).any():
        exons = gtf[(gtf['feature'].str.lower() == 'three_prime_utr') &
                    (gtf['attributes'].str.contains("MANE_Select", case=False, na=False))]
    else:
        exons = gtf[gtf['feature'].str.lower() == 'three_prime_utr']
    
    # Deduplicate and sort
    exons = exons[['chrom','start', 'end']].drop_duplicates()
    exons_list = list(zip(exons['chrom'],exons['start'], exons['end']))
    return sorted(exons_list, key=lambda x: x[0])


def adjust_zero_cov(zero_cov, exons, strand):
    """
    Adjust the zero_cov coordinate so that it falls within an exon.
    If zero_cov already lies in an exon, return it as is.
    Otherwise, return the nearest exon boundary.
    """
    # Check if zero_cov is within any exon.
    for start, end in exons:
        if start <= zero_cov < end:
            return zero_cov

    # If not, adjust it.
    if strand == '+':
        # If zero_cov is before the first exon, use its start.
        if zero_cov < exons[0][0]:
            return exons[0][0]
        # If after the last exon, use its end.
        if zero_cov > exons[-1][1]:
            return exons[-1][1]
        # Otherwise, zero_cov is in an intron between two exons.
        for i in range(len(exons) - 1):
            if exons[i][1] <= zero_cov < exons[i+1][0]:
                # Snap to the boundary that is closer.
                if (zero_cov - exons[i][1]) < (exons[i+1][0] - zero_cov):
                    return exons[i][1]
                else:
                    return exons[i+1][0]
        return zero_cov  # Fallback (should not reach here)
    else:
        # For minus strand, the genomic coordinates are the same,
        # so use the same logic.
        if zero_cov > exons[-1][1]:
            return exons[-1][1]
        if zero_cov < exons[0][0]:
            return exons[0][0]
        for i in range(len(exons) - 1):
            if exons[i][1] <= zero_cov < exons[i+1][0]:
                if (zero_cov - exons[i][1]) < (exons[i+1][0] - zero_cov):
                    return exons[i][1]
                else:
                    return exons[i+1][0]
        return zero_cov

def compute_spliced_distance(primer_genomic, zero_cov_adj, exons, strand):
    """
    Compute the spliced (introns excluded) distance between the primer genomic coordinate
    and the adjusted zero coverage coordinate by summing the overlaps with exons.
    
    For both strands, we simply consider the genomic interval between the two positions and
    sum up the portions that fall into exons.
    """
    # Determine the bounds of the genomic interval.
    left_bound = min(primer_genomic, zero_cov_adj)
    right_bound = max(primer_genomic, zero_cov_adj)
    
    spliced = 0
    for start, end in exons:
        # Skip exons that lie completely before or after the interval.
        if end <= left_bound:
            continue
        if start >= right_bound:
            break
        # Calculate overlap.
        overlap = max(0, min(end, right_bound) - max(start, left_bound))
        spliced += overlap
    return spliced

def calculate_amplicon_size(position_primer, strand, zero_cov_coord, start_expanded, end_expanded, output_dir, gene_id):
    """
    Calculate the amplicon size along the transcript (i.e. excluding introns)
    by (1) converting the primer position (reported by Primer3, relative to the FASTA)
    into a genomic coordinate, (2) adjusting the zero coverage coordinate to lie within an exon,
    and then (3) summing the exonic bases between them.
    """
    out_file_ampl = f"{output_dir}/{gene_id}_amplicon_size.txt"
    pos_start = int(position_primer.split(',')[0])
    
    # Map the primer position (relative to the FASTA region) to its genomic coordinate.
    # The FASTA region was extracted from [start_expanded, end_expanded].
    if strand == "-":
        primer_genomic = end_expanded - pos_start
    else:
        primer_genomic = start_expanded + pos_start

    # Retrieve the individual exon intervals from the filtered GTF.
    exons_list = get_mane_exons(output_dir, gene_id)  # Should return a sorted list of (start, end) tuples.
    if not exons_list:
        print(f"No exons found for {gene_id}. Cannot compute spliced amplicon size.")
        return

    # Adjust the zero coverage coordinate to fall within an exon.
    zero_cov_adj = adjust_zero_cov(zero_cov_coord, exons_list, strand)
    print(f"For {gene_id}, zero_cov_coord {zero_cov_coord} adjusted to {zero_cov_adj} based on exons: {exons_list}")

    # Compute the spliced distance by summing the overlapping exonic segments.
    spliced_distance = compute_spliced_distance(primer_genomic, zero_cov_adj, exons_list, strand)
    
    with open(out_file_ampl, "w") as f:
        f.write(f"{gene_id} {spliced_distance}\n")
    print(f"Calculated spliced amplicon size for {gene_id}: {spliced_distance} bp")


def run_blast_local_command2(sequence, database="core_nt", output_file="blast_results.xml"):
    input_file = "query.fasta"
    with open(input_file, "w") as f:
        f.write(f">query\n{sequence}")
    command = f"ml purge; ml BLAST+/2.14.1-gompi-2023a; blastn -query {input_file} -db {database} -out {output_file} -outfmt 5 -word_size 7 -reward 1 -penalty -3 -gapopen 5 -gapextend 2"
    return command

def run_blast_remote_command2(sequence, output_file="blast_results.xml"):
    """
    Generates a command to run a local BLAST search for a given sequence.
    """
    input_file = "query.fasta"
    with open(input_file, "w") as f:
        f.write(f">query\n{sequence}")
    command = f"ml purge; ml BLAST+/2.14.1-gompi-2023a; blastn -query {input_file} -db core_nt -out {output_file} -outfmt 5 -word_size 7 -reward 1 -penalty -3 -gapopen 5 -gapextend 2"
    return command



def main(bam_file, gtf_file, gene_id, TM, ref_genome, window, cov_thresh, max_primer, min_primer, max_product_size, min_product_size, repeat_lib, db, out_dir):
    # Create an output directory for this gene.
    output_dir = f"{out_dir}/{gene_id}_out"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Generate genome coverage.
    coverage_command = get_coverage(bam_file, output_dir)
    if coverage_command:
        sbatch("get_coverage", coverage_command, 1)

    # Step 2: Filter the GTF for the gene.
    filter_gtf_command = filter_gtf(gtf_file, gene_id, output_dir)
    
    if filter_gtf_command:
        sbatch("filter_gtf", filter_gtf_command, 2)

    # Step 2.5: Merge the exons for primer design.
    merge_exons_command = merge_exons(output_dir, gene_id)
    if merge_exons_command:
        sbatch("merge_exons", merge_exons_command, 2.5)

    # Get strand information from the filtered GTF.
    strand = get_strand_from_gtf(gene_id, output_dir)
    strand_label = "pos" if strand=="+" else "neg"

    # Step 3: Intersect coverage with the merged exons.
    gene_coverage_command = get_gene_coverage(output_dir, gene_id)
    sbatch("get_gene_coverage", gene_coverage_command, 3)

    # Step 4: Determine the peak coverage and the zero coverage coordinate.
    chrom, start_expanded, end_expanded, peak_coord, zero_cov_coord = get_peak_and_low_coverage_prct(window, output_dir, gene_id, strand, coverage_percentage=cov_thresh)

    # Check for overlap with 3' UTR.
    if check_3prime_utr(output_dir, gene_id, chrom, start_expanded, end_expanded):
        UTR = "utr"
        print(f"Gene {gene_id} has a 3' UTR in the region {chrom}:{start_expanded}-{end_expanded}")
    else:
        UTR = "no_utr"
        print(f"Gene {gene_id} does not have a 3' UTR in the region {chrom}:{start_expanded}-{end_expanded}")

    # Step 5: Retrieve the genomic sequence.
    fasta_command = get_fasta_sequence(chrom, start_expanded, end_expanded, ref_genome, output_dir, gene_id)
    sbatch("get_fasta_sequence", fasta_command, 5)

    # Step 6: Create Primer3 input.
    
    fasta_input = f"{output_dir}/{gene_id}_{chrom}_{start_expanded}_{end_expanded}.fa"
    primer3_input_file = f"{output_dir}/{gene_id}_{UTR}_primer3_input"
    primer3_input_path = write_primer3_input2(gene_id, fasta_input, primer3_input_file, strand, TM, min_primer, max_primer, min_product_size, max_product_size, repeat_lib)
    print(f"Primer3 input file written to: {primer3_input_path}")

    # Step 7: Run Primer3.
    primer3_output_file = f"{output_dir}/{gene_id}_{strand_label}_primer3_out.txt"
    primer3_command = Run_Primer3(primer3_input_path, primer3_output_file)
    sbatch("Primer3", primer3_command, 7)
    gene, primer_seqs, primer_positions, TM, GC = parse_primer3_output(primer3_output_file)
    PrimerFasta_output_file = f"{out_dir}/primers.fa"
    write_primers_to_fasta(gene, primer_seqs, PrimerFasta_output_file)

    # Step 8: Calculate the spliced amplicon size.
    #calculate_amplicon_size(primer_positions, strand, zero_cov_coord, start_expanded, end_expanded, output_dir, gene_id)

    # Step 9: Run local BLAST for each primer.
    # primer_seqs is a dictionary where keys are primer indices (e.g., 'primer_0', 'primer_1', ...)
    # and values are the corresponding primer sequences (strings).
    if not primer_seqs:
        print(f"No primers found in the Primer3 file {primer3_output_file}")
        return None
    else:
        for primer_name, primers in primer_seqs.items():
            print(f"Running BLAST for {primer_name} primer: {primers}")
            # db = "/user/gent/446/vsc44685/ScratchVO_dir/NSQ2K_126_test/primer_design/blast_corent/core_nt"
            out_file = f"{output_dir}/{gene_id}_{primer_name}_blast.xml"
            blast = run_blast_local_command2(primers, db, out_file)
            # blast = run_blast_remote_command2(primer_sequence, out_file)
            sbatch("blast", blast, 8)
