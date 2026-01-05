from Bio import SeqIO
from collections import Counter
import csv

def calculate_cpb(fasta_files, cps_file, output_file):
    # Load the CPS scores into a dictionary
    cps_scores = {}
    with open(cps_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cps_scores[row["CodonPair"]] = float(row["CPS"])

    # Open the output file
    with open(output_file, "w") as out_file:
        # Process each fasta file
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = str(record.seq)
                codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
                codon_pairs = [(codons[i], codons[i+1]) for i in range(len(codons)-1)]
                codon_pairs_counter = Counter(codon_pairs)

                # Calculate the CPB for the sequence
                cpb = sum(count * cps_scores.get(pair[0]+pair[1], 0) for pair, count in codon_pairs_counter.items()) / len(codon_pairs)

                # Write the result to the output file
                out_file.write(f"{fasta_file}, {record.id}, {record.description.split()[1]}, {cpb}\n")

# Usage
fasta_files = ["SARS2_S_WT.fasta"]
cps_file = "CPS_human.csv"
output_file = "output_SARS2_S_WT_CPB.txt"
calculate_cpb(fasta_files, cps_file, output_file)