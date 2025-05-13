#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import pysam
import argparse

def make_alignment_summary(transcriptome, sam_file, out_file):
    # Create a lookup dictionary for transcript annotations from the FASTA file.
    transcript_info = {}
    with open(transcriptome) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # record.id is typically the transcript ID (e.g., "ENST00000456328")
            # record.description contains the full header (e.g., "ENST00000456328 gene=DDX11L2")
            transcript_info[record.id] = record.description

    samfile = pysam.AlignmentFile(sam_file, "r")
    ref_lengths = {sq['SN']: sq['LN'] for sq in samfile.header['SQ']}
    # Open an output file to write the tab-delimited summary.
    with open(out_file, "w") as outfile:
        # Write a header row including start and end positions, alignment length, and mismatches.
        outfile.write("query_name\tPrimer_Index\tmap_to_transcript\tmap_to_gene_name\tstart_align\tend_align\talign_length\tmismatch\tdirection\tDistance_to_end\n")
        
        # Iterate through each alignment record in the SAM file.
        for read in samfile:
            gene_id = read.query_name.split("_")[0]
            primer_id = "_".join(read.query_name.split("_")[2:])
            transcript_id = read.reference_name
            transcript_length = ref_lengths.get(transcript_id, "NA")
            # Get the transcript annotation from the dictionary.
            annotation = transcript_info.get(transcript_id, "Annotation not found")
            # Extract gene_name from the annotation.
            gene_name = "NA"
            if annotation != "Annotation not found" and "gene=" in annotation:
                gene_name = annotation.split("gene=")[1].split()[0]
            # Get the aligned portion length.
            alignment_length = read.query_alignment_length
            # Get the start and end positions of the alignment.
            start = read.reference_start
            end = read.reference_end if read.reference_end is not None else "NA"
            #Compute distance to the end of the transcript.
            if transcript_length != "NA":
                distance_to_end = transcript_length - start
            else:
                distance_to_end = "NA"
            # Determine the alignment direction.
            # If read.is_reverse is True, the read aligns to the reverse strand.
            direction = "reverse" if read.is_reverse else "forward"
            # Extract mismatches from the NM tag (if present).
            mismatches = read.get_tag("NM") if read.has_tag("NM") else "NA"
            
            # Write one row per alignment.
            outfile.write(f"{gene_id}\t{primer_id}\t{transcript_id}\t{gene_name}\t{start}\t{end}\t{alignment_length}\t{mismatches}\t{direction}\t{distance_to_end}\n")

    samfile.close()
    return

parser = argparse.ArgumentParser(description="Process Sam File and outputs a summary.")
parser.add_argument('--transcriptome', required=True, help='output of parse transcriptome')
parser.add_argument('--sam_file', help='path to save sam_file file')
parser.add_argument('--out_file', help='path to save summary file')

args = parser.parse_args()
transcriptome = args.transcriptome
sam_file = args.sam_file
out_file = args.out_file


if __name__ == "__main__":
    make_alignment_summary(transcriptome, sam_file, out_file)
