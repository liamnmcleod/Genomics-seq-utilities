# Genomics-seq-utilities
Python utilities for working with FASTA/FASTQ, BLAST, and ORF parsing.


Small, readable Python utilities for working with FASTA/FASTQ, running simple BLAST queries, and parsing ORFs.

What’s inside

bio_blast_ncbi.py – helper to submit sequences to NCBI BLAST and capture results (e.g., BLASTN/BLASTP).

fasta_seq_to_dictionary.py / fasta_to_dict.py – parse a FASTA into {header → sequence} and (optionally) filter by minimum length.

fastq_analyzer.ipynb – QC on FASTQ (read length distribution, base composition, quality histogram, GC by position).

multi_fasta_orf_parser.py – scan multi-FASTA files for putative open reading frames and report coordinates/sequences.
