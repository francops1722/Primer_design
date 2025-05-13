#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 14:57:30 2025

@author: fpomasot
"""
#Run main pipeline

import os
import subprocess
import pandas as pd
import argparse
import requests
import re
from Bio import SeqIO

from utils import (
    main, sbatch, run_bowtie2_command, run_make_alignment_summary 
)


# Parse command line arguments.
parser = argparse.ArgumentParser(description="Process BAM and GTF files for a specific gene.")
parser.add_argument('--gene_ids', required=True, help='Path to the file containing gene IDs')
parser.add_argument('--Tm', type=float, default=60.0, help='Melting temperature for primers')
parser.add_argument('--min_primer_size', type=int, default=18, help='Min primer size')
parser.add_argument('--max_primer_size', type=int, default=18, help='Max primer size')
parser.add_argument('--gtf_file', required=True, help='Path to the GTF file')
parser.add_argument('--btw_index', default="/data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index", help='Path to the bowtie2 index prefix')
parser.add_argument('--transcriptome', required=True, default="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa", help="path to transcriptome fasta")
parser.add_argument('--dist', type=int, default=100, help='Distance to the end of the transcript')
#parser.add_argument('--min_product_size', type=int, default=100, help='Min product size for primer design')
parser.add_argument('--repeat_lib', default='.', help='path to repeat library for masking')
parser.add_argument('--out_dir', default='.', help='Output directory for results')

args = parser.parse_args()
TM = args.Tm
gene_ids_file = args.gene_ids
gtf_file = args.gtf_file
ampl_size = args.dist
out_dir = args.out_dir
index = args.btw_index
lib=args.repeat_lib
transcriptome = args.transcriptome
min_primer_size=args.min_primer_size
max_primer_size=args.max_primer_size

with open(gene_ids_file, 'r') as file:
    gene_ids = [line.strip() for line in file]
if __name__ == "__main__":
    for gene_id in gene_ids:
        sum_df = main(out_dir, gtf_file, gene_id, ampl_size, TM, lib, min_primer_size, max_primer_size)
    #run bowtie2
    PrimerFasta_output_file = f"{out_dir}/primers.fa"
    sam_file=f"{out_dir}/bowtie2_results.sam"
    bowtie2_command = run_bowtie2_command(PrimerFasta_output_file, index, output_file=sam_file)
    sbatch("bowtie2", bowtie2_command, 7)
    sam_summary=f"{out_dir}/Alignment_summary.tsv"
    bowtie_summary = run_make_alignment_summary(transcriptome, sam_file, sam_summary)
    sbatch('Bowtie_summary', bowtie_summary, 8)
    print("All steps finished.")


    