
from utils_informed import main
import argparse
# Parse command line arguments.
parser = argparse.ArgumentParser(description="Process BAM and GTF files for a specific gene.")
parser.add_argument('--gene_ids', required=True, help='Path to the file containing gene IDs')
parser.add_argument('--bam_file', required=False, help='Path to the BAM file')
parser.add_argument('--Tm', type=float, default=60.0, help='Melting temperature for primers')
parser.add_argument('--min_primer_size', type=int, default=18, help='Min primer size')
parser.add_argument('--max_primer_size', type=int, default=20, help='Max primer size')
parser.add_argument('--min_product_size', type=int, default=18, help='Min product size')
parser.add_argument('--max_product_size', type=int, default=60, help='Max product size')
parser.add_argument('--gtf_file', required=True, help='Path to the GTF file')
parser.add_argument('--blastDB', required=True, default="/data/gent/vo/000/gvo00027/resources/Blast_databases/RefSeq_select/human_refseq_db", help='Path to the a local database for BLAST')
parser.add_argument('--ref_genome', required=True, help='Path to the reference genome file')
parser.add_argument('--window', type=int, default=100, help='Window size for peak coverage expansion')
parser.add_argument('--dist', type=int, default=100, help='Distance to the end of the UTR')
parser.add_argument('--out_dir', default='.', help='Output directory for results')
parser.add_argument('--cov_thresh', required=False, type=float, default=5, help='coverage threshold to compute the low coverage region')
parser.add_argument('--repeat_lib', default='"./rep_lib/humrep_and_simple.txt"', help='path to repeat library for masking')

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

with open(gene_ids_file, 'r') as file:
    gene_ids = [line.strip() for line in file]

if __name__ == "__main__":
    for gene_id in gene_ids:
        if bam_file:
            main(bam_file, gtf_file, gene_id, ref_genome, window, cov_thresh, min_primer_size, max_primer_size, max_product_size, min_product_size, TM, db, output_dir, repeat_lib)
