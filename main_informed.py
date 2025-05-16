
from utils_informed import main
from utils import (sbatch, run_bowtie2_command, run_make_alignment_summary)
import argparse
# Parse command line arguments.
parser = argparse.ArgumentParser(description="Process BAM and GTF files for a specific gene.")
parser.add_argument('--gene_ids', required=True, help='Path to the file containing gene IDs')
parser.add_argument('--bam_file', required=False, help='Path to the BAM file')
parser.add_argument('--Tm', type=float, default=60.0, help='Melting temperature for primers')
parser.add_argument('--min_primer_size', type=int, default=18, help='Min primer size')
parser.add_argument('--max_primer_size', type=int, default=18, help='Max primer size')
parser.add_argument('--min_product_size', type=int, default=18, help='Min product size')
parser.add_argument('--max_product_size', type=int, default=18, help='Max product size')
parser.add_argument('--gtf_file', required=True, help='Path to the GTF file')
parser.add_argument('--btw_index', default="/data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index", help='Path to the bowtie2 index prefix')
parser.add_argument('--transcriptome', default="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa", help='Path to transcriptome(to get annotation of each transcirpt)')
parser.add_argument('--blastDB', required=True, default="/data/gent/vo/000/gvo00027/resources/coreNT_db/blast_corent/core_nt", help='Path to the a local database for BLAST')
parser.add_argument('--ref_genome', required=True, help='Path to the reference genome file')
parser.add_argument('--window', type=int, default=100, help='Window size for peak coverage expansion')
parser.add_argument('--dist', type=int, default=100, help='Distance to the end of the UTR')
parser.add_argument('--repeat_lib', default="./rep_lib/humrep_and_simple.txt", help='path to repeat library for masking')
parser.add_argument('--out_dir', default='.', help='Output directory for results')
parser.add_argument('--cov_thresh', required=False, type=float, default=5, help='coverage threshold to compute the low coverage region')

args = parser.parse_args()
bam_file = args.bam_file
TM = args.Tm
gene_ids_file = args.gene_ids
gtf_file = args.gtf_file
ref_genome = args.ref_genome
window = args.window
distance = args.dist
output_dir = args.out_dir
cov_thresh = args.cov_thresh
db = args.blastDB
max_primer_size=args.max_primer_size
min_primer_size=args.min_primer_size
max_product_size=args.max_product_size
min_product_size=args.min_product_size
repeat_lib=args.repeat_lib
index = args.btw_index
transcriptome = args.transcriptome
# transcriptome = "/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa"
with open(gene_ids_file, 'r') as file:
    gene_ids = [line.strip() for line in file]

# if __name__ == "__main__":
    for gene_id in gene_ids:
        if bam_file:
            main(bam_file, gtf_file, gene_id, TM, ref_genome, window, cov_thresh, max_primer=max_primer_size, min_primer=min_primer_size, max_product_size=max_product_size, min_product_size=min_product_size,db=db, repeat_lib=repeat_lib, out_dir=output_dir)
    #run bowtie2
    PrimerFasta_output_file = f"{output_dir}/primers.fa"
    sam_file=f"{output_dir}/bowtie2_results.sam"
    bowtie2_command = run_bowtie2_command(PrimerFasta_output_file, index, output_file=sam_file)
    sbatch("bowtie2", bowtie2_command, 9)
    sam_summary=f"{output_dir}/Alignment_summary.tsv"
    bowtie_summary = run_make_alignment_summary(transcriptome, sam_file, sam_summary)
    sbatch('Bowtie_summary', bowtie_summary, 8)
    print("All steps finished.")
            
