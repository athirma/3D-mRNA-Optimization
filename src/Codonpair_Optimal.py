import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# Read codon pair replacement list from CSV file
# Format: Column 0 contains codon pairs to be replaced, Column 1 contains replacement codon pairs
df = pd.read_csv("replace_codonpair.csv", header=None)

# Create dictionary mapping original codon pairs to optimized replacements
replacement_dict = dict(zip(df[0], df[1]))

# Process FASTA files iteratively for codon pair optimization
for idx in range(len(replacement_dict)):
    # Determine input FASTA file name
    # For first iteration, use wild-type sequence; for subsequent iterations, use previous optimization result
    input_file = f"test_{idx}.fasta" if idx > 0 else "SARS2_S_WT.fasta"

    # Read and parse FASTA file containing nucleotide sequences
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Process each sequence in the FASTA file
    for seq_record in sequences:
        # Convert Bio.Seq object to string for codon-level manipulation
        seq_str = str(seq_record.seq)

        # Divide nucleotide sequence into consecutive codon triplets
        original_codons = [seq_str[i:i + 3] for i in range(0, len(seq_str), 3)]

        # Create copy of codon list for modification (preserve original for reference)
        optimized_codons = original_codons.copy()

        # Iterate through adjacent codon pairs in the sequence
        for i in range(len(original_codons) - 1):
            # Concatenate adjacent codons to form codon pair (hexanucleotide)
            codon_pair = original_codons[i] + original_codons[i + 1]

            # Check if current codon pair matches target for replacement at this iteration
            if codon_pair == df[0][idx]:
                # Replace with optimized codon pair: split 6-nt replacement into two 3-nt codons
                optimized_codons[i] = df[1][idx][:3]  # First codon of replacement pair
                optimized_codons[i + 1] = df[1][idx][3:]  # Second codon of replacement pair

        # Update sequence record with optimized codon sequence
        seq_record.seq = Seq(''.join(optimized_codons))

    # Write optimized sequences to new FASTA file for next iteration
    output_file = f"test_{idx + 1}.fasta"
    SeqIO.write(sequences, output_file, "fasta")

# Rename final optimization result to indicate completion
final_output = f"test_{len(replacement_dict)}.fasta"
os.rename(final_output, "SARS2_S_Opi_replace.fasta")

# Cleanup: Remove intermediate files generated during iterative optimization
for idx in range(len(replacement_dict)):
    temp_file = f"test_{idx}.fasta"
    if os.path.exists(temp_file):
        os.remove(temp_file)
