from pyexpat.errors import messages

import sys
import getopt

def usage():
    print("""" fasta_to_dictionary.py reads a FASTA file and builds a dictionary with all sequences larger 
    than a specified length
    
    fasta_to_dictionary.py [-h] [-l <length>] <filename>
    -h                  print this message
    -l <length>         filter out all sequences smaller than <length>, default <length> = 0
    <filename>          the file has to be in FASTA format

    """)

o, a = getopt.getopt(sys.argv[1:], "l:h") # o = list of optional arguments, a = list of required arguments

seqs = {}

try:
    with open(r"C:\Users\liamm\Desktop\testfasta.fasta") as fasta_file:
        for line in fasta_file:
            line = line.rstrip() #removes any trailing whitespace
            if not line:
                continue
            if line.startswith(">"):    #checks for header line at position 0
                name = line[1:].split()[0]  #if header then remove ">"
                seqs[name] = ""
            else:
                if name is None:
                    raise ValueError("Sequence data found before a header.")
                seqs[name] += line  # add sequence to current name
except IOError:
    print("file not found")
for name, sequence in seqs.items():
    print(f"{name}: {sequence}")
