from Bio import SeqIO
from collections import Counter
import argparse


def parse_fasta_to_dict(filepath):
    """
    Reads a FASTA file and returns its contents as a dictionary.
    Uses Bio.SeqIO for robust parsing.

    Args:
        filepath (str): The path to the FASTA file.

    Returns:
        dict: A dictionary of {sequence_id: sequence_string}.
    """
    try:
        # Use a dictionary comprehension for a concise and pythonic way to build the dictionary
        return {record.id: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None


def print_sequence_summary(sequences):
    """Prints the number of sequences and identifies the shortest and longest ones."""
    print(f"Number of sequences: {len(sequences)}")
    print("-" * 80)

    if not sequences:
        return

    # Find the IDs of the min and max length sequences
    min_id = min(sequences, key=lambda seq_id: len(sequences[seq_id]))
    max_id = max(sequences, key=lambda seq_id: len(sequences[seq_id]))

    print(f"Shortest sequence: {min_id} ({len(sequences[min_id])} bp)")
    print(f"Longest sequence:  {max_id} ({len(sequences[max_id])} bp)")
    print("-" * 80)


def find_longest_orfs(sequences):
    """
    Finds the longest Open Reading Frame (ORF) in reading frames 1, 2, and 3.

    Args:
        sequences (dict): A dictionary of sequences.

    Returns:
        dict: A dictionary mapping each sequence ID to its longest ORF info.
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    longest_orfs = {}

    for seq_id, sequence in sequences.items():
        max_orf_len = 0
        max_orf_info = (0, None, None)  # (length, frame, start_position)

        for frame in range(3):
            for i in range(frame, len(sequence), 3):
                codon = sequence[i:i + 3]

                # Search for a start codon to begin a potential ORF
                if codon == start_codon:
                    # Scan forward from the start codon for a stop codon
                    for j in range(i + 3, len(sequence), 3):
                        next_codon = sequence[j:j + 3]
                        if next_codon in stop_codons:
                            orf_len = j + 3 - i
                            if orf_len > max_orf_len:
                                max_orf_len = orf_len
                                max_orf_info = (orf_len, frame + 1, i)  # Store frame as 1,2,3
                            break  # Found the first stop codon, end this ORF
        longest_orfs[seq_id] = max_orf_info
    return longest_orfs


def count_sequence_repeats(filepath, repeat_length):
    """
    Counts all overlapping repeats of a given length in a FASTA file.

    Args:
        filepath (str): The path to the FASTA file.
        repeat_length (int): The length of the repeat to search for.

    Returns:
        tuple: (A Counter object with all repeats, a list of the most frequent repeats, the highest count).
    """
    repeat_counts = Counter()
    for record in SeqIO.parse(filepath, "fasta"):
        seq = str(record.seq)
        if len(seq) >= repeat_length:
            for i in range(len(seq) - repeat_length + 1):
                repeat = seq[i:i + repeat_length]
                repeat_counts[repeat] += 1

    if not repeat_counts:
        return repeat_counts, [], 0

    max_count = max(repeat_counts.values())
    most_frequent = [repeat for repeat, count in repeat_counts.items() if count == max_count]
    return repeat_counts, most_frequent, max_count


def main():
    """Main function to run the DNA sequence analysis."""
    # Setup argument parser for command-line usage
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("repeat_length", type=int, help="Length of the repeat sequence to search for.")
    args = parser.parse_args()

    sequences = parse_fasta_to_dict(args.fasta_file)
    if sequences is None:
        return  # Exit if file reading failed

    print_sequence_summary(sequences)

    orf_summary = find_longest_orfs(sequences)
    print("Longest ORF Analysis:")
    for seq_id, (length, frame, start) in orf_summary.items():
        if frame is not None:
            print(f"  {seq_id}: Length={length} bp, Frame={frame}, Start={start}")
        else:
            print(f"  {seq_id}: No ORF found")
    print("-" * 80)

    repeats, most_frequent, max_count = count_sequence_repeats(args.fasta_file, args.repeat_length)
    print(f"Most Frequent {args.repeat_length}-mer Repeat(s) (occurred {max_count} times):")
    print(f"  {', '.join(most_frequent)}")
    print("-" * 80)


if __name__ == "__main__":
    main()