# Primer_design
A python pipeline for high-throughput 3'end primer design

## Inputs:

- **gene_id:** A list of ensembl gene ids in a text file, example ./testdata/gene_ids.txt
- **gtf_file**: A GTF with gene annotation to extract the MANE transcript information
- **transcriptome**: A fasta file of the Transcriptome to be used for checking primer specificity
- **btw_index**: A bowtie index bulit on the same transcriptome, step not included in the pipeline, has to be created prior to run the pipeline
- **dist**: Parameter that indicates how many nucleotide are taken from the 3'end of each transcript to be used as template for primer3
- **Tm**: melting temperature for primer 3
- **--min_primer_size and --max_primer_size**: min and max Primer size
- **Repeat_lib**: a fasta file with repetitive sequences to exclude when designing primers.

## Workflow:

1. For each gene ID get the MANE transcript ID
2. Fetch the https://rest.ensembl.org site and gets the transcript sequence
3. Write a text file with all needed parameters to run Primer3
4. Run Primer3 and saves output files 
5. Run Bowtie2 enabling multiple alignments report
6. Makes summary files. 

## How to Run

1. Pull Primer3 singularity container

        `cd ./primer3`

        `./PullSif.sh`

2. Edit input parameters the job.pbs file
3. Submit job

        qsub Job.pbs