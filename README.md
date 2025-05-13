# Primer_design

The repository contains scripts to design only left primers using 3 modes: 

- **main.py:** script to design primers give a distance from the 3'end
- **main_informed.py:** script to design primers in an informed manner, using alignments in a bam file to target regions with high coverage
- **main_5end:** script to generate primers at 5'end, aimed for full-length amplification   

A python pipeline for high-throughput 3'end primer design

## Inputs:

For 3'end designing of primers:

- **gene_id:** A list of ensembl gene ids in a text file, example ./testdata/gene_ids.txt
- **gtf_file**: A GTF with gene annotation to extract the MANE transcript information
- **transcriptome**: A fasta file of the Transcriptome to be used for checking primer specificity
- **btw_index**: A bowtie index bulit on the same transcriptome, step not included in the pipeline, has to be created prior to run the pipeline
- **dist**: Parameter that indicates how many nucleotide are taken from the 3'end of each transcript to be used as template for primer3
- **Tm**: melting temperature for primer 3
- **--min_primer_size and --max_primer_size**: min and max Primer size
- **Repeat_lib**: a fasta file with repetitive sequences to exclude when designing primers.

For Informed designing of primers:

- **bam_file:** path to bam file to be used to define peaks of coverage 
- **Tm**: melting temperature for primer 3
- **--min_primer_size and --max_primer_size**: min and max Primer size
- **--min_product_size and --max_product_size**: min and max Product size
- **gtf_file:** A GTF with gene annotation to extract the MANE transcript information
- **blastDB:**  path to a local saved database, to run blast against
- **ref_genome:** path to reference genome fasta file
- **window:** given a window N, the pipeline detects the coordinate of highest coverage and extracts the genome sequence of length 2xN around the peak coordinate  
- **cov_thresh** parameter used to compute an expected product length. Finds the first low-coverage (coverage < cov_thresh) coordinate downstream/upstream the peak coordinate
- **Repeat_lib**: a fasta file with repetitive sequences to exclude when designing primers.
    

## Workflow:

### 3'end

1. For each gene ID get the MANE transcript ID
2. Fetch the https://rest.ensembl.org site and gets the transcript sequence
3. Write a text file with all needed parameters to run Primer3
4. Run `Primer3` and saves output files 
5. Run `Bowtie2` enabling multiple alignments report
6. Makes summary files. 

### Coverage Informed 

1. **Coverage Calculation**
2. **Gene Processing**

   * Filter GTF: Extract entries for the target gene
   * Merge exons: Combine exons into a continuous region
   * Determine transcriptional strand: Identify as (+) or (-)

3. **Peak Detection**

   * Identify peak: Find the maximum coverage interval
   * Expand region: Use `window_size` to extend around the peak
   * Dynamic threshold: Calculate `cov_thresh%` of the peak value
   * Boundary detection: Find the first low-coverage coordinate downstream/upstream (strand-aware)

4. **Primer Design**

   * Extract sequence: Retrieve genomic sequence using `SAMtools`
   * Strand handling: Reverse complement if on the negative strand
   * Prepare Primer3 input 
   * Run Primer3: Execute within a singularity container

5. **Amplicon Calculation**

   * Coordinate adjustment: Align coordinates with exons
   * Spliced length: Compute the length considering exons only
   * Report generation: Create an amplicon size report

6. **Specificity Check**

   * BLAST: Search primers against a local database

## How to Run

1. Clone the repository

        git clone https://github.com/francops1722/Primer_design.git

2. Pull Primer3 singularity container

        cd ./primer3

        ./PullSif.sh

2. Edit input parameters the job.pbs files:
        - Job_3end.pbs
        - Job_Informed.pbs
        - Job_5end.pbs
3. Submit job

        qsub Job.pbs
