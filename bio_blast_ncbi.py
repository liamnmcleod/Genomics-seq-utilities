from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

FASTA_FILE = r"C:\Users\liamm\Desktop\myseq.FASTA"
E_VALUE_THRESH = 0.01

try:
    with open(FASTA_FILE, "r") as f:
        fasta_string = f.read()
except FileNotFoundError:
    print(f"Error: The file '{FASTA_FILE}' was not found.")
    exit()


print(f"Performing BLASTN search for sequence in {FASTA_FILE}...")
# This can take a few moments depending on the server load.
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
print("Search complete. Parsing results...")


with open("blast_results.xml", "w") as out_file:
    out_file.write(result_handle.read())
result_handle.close()


blast_record = NCBIXML.read(open("blast_results.xml"))

print("\n--- Significant Alignments ---")
# Loop through alignments and High-scoring Segment Pairs (HSPs).
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        # Check if the E-value is below our threshold.
        if hsp.expect < E_VALUE_THRESH:
            print("\n**** Alignment ****")
            # Use f-strings for cleaner, labeled output.
            print(f"Sequence: {alignment.title}")
            print(f"Length: {alignment.length} bp")
            print(f"E-value: {hsp.expect}")
            print(f"Score: {hsp.score} bits")


            print("\nAlignment Details:")
            print(f"Query: {hsp.query[0:75]}...")
            print(f"Match: {hsp.match[0:75]}...")
            print(f"Sbjct: {hsp.sbjct[0:75]}...")
            print("-" * 30)

print("Analysis finished")