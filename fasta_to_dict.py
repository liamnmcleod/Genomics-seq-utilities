
"""Reads a FASTA file and builds a dictionary with all sequences,
optionally filtering by a minimum length."""

import argparse
import sys

def parse_fasta(filepath):
    sequences = {}
    current_header = None
    try:
        with open(filepath, 'r') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # First part of the header as the key
                    current_header = line[1:].split()[0]
                    sequences[current_header] = []
                elif current_header:
                    sequences[current_header].append(line)
    except IOError as e:
        print(f"Error: Cannot read file '{filepath}'.\n{e}", file=sys.stderr)
        return None


    final_sequences = {header: "".join(parts) for header, parts in sequences.items()}
    return final_sequences

def main():
    """Main function to parse arguments and print sequences."""
    parser = argparse.ArgumentParser(
        description="Reads a FASTA file into a dictionary, filtering by sequence length."
    )
    parser.add_argument(
        "filename",
        help="The FASTA file to be read."
    )
    parser.add_argument(
        "-l", "--length",
        type=int,
        default=0,
        help="Minimum length of sequences to keep. Default is 0 (keep all)."
    )
    args = parser.parse_args()

    all_seqs = parse_fasta(args.filename)

    if all_seqs is not None:
        # Filter dictionary based on the length argument using a comprehension
        filtered_seqs = {
            name: seq for name, seq in all_seqs.items() if len(seq) >= args.length
        }

        for name, sequence in filtered_seqs.items():
            print(f">{name}\n{sequence}")

if __name__ == "__main__":
    main()