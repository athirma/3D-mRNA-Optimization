"""
Nucleotide and dinucleotide composition of mRNA coding sequences.

Reports, for each sequence: length, A/C/G/U content, GC and GC3 content,
CpG and UpA counts with per-kilobase rates and observed/expected ratios, and
optionally the codon pair bias (CPB) score.

The observed/expected ratio normalises a dinucleotide count for the sequence's
own mononucleotide composition and for its length, which is what makes counts
comparable between constructs of different GC content and between antigens of
different coding-sequence length:

    O/E(XY) = (count_XY * L) / (count_X * count_Y)

A value of 1 is the frequency expected if bases were arranged at random; values
below 1 indicate suppression of that dinucleotide.

Usage
-----
    python src/Composition_Analysis.py examples/SARS2_S_WT.fasta \\
        --cps data/CPS_human.csv --output composition.csv
"""

import argparse
import csv

from Bio import SeqIO

BASES = "ACGT"


def load_cps(path):
    """Read the codon pair score table used by CPB_Calculation.py."""
    scores = {}
    with open(path, "r") as handle:
        for row in csv.DictReader(handle):
            scores[row["CodonPair"].strip().upper()] = float(row["CPS"])
    return scores


def dinucleotide_oe(sequence, first, second):
    """Return the count and observed/expected ratio of a dinucleotide."""
    length = len(sequence)
    observed = sum(1 for i in range(length - 1)
                   if sequence[i] == first and sequence[i + 1] == second)
    n_first = sequence.count(first)
    n_second = sequence.count(second)
    if not n_first or not n_second:
        return observed, float("nan")
    return observed, (observed * length) / (n_first * n_second)


def analyse(name, sequence, cps_scores=None):
    """Compute composition statistics for one coding sequence."""
    sequence = sequence.upper().replace("U", "T")
    length = len(sequence)
    if length == 0:
        raise ValueError(f"{name}: empty sequence")

    counts = {base: sequence.count(base) for base in BASES}
    n_codons = length // 3
    third_positions = sequence[2:n_codons * 3:3]

    cpg_n, cpg_oe = dinucleotide_oe(sequence, "C", "G")
    upa_n, upa_oe = dinucleotide_oe(sequence, "T", "A")

    result = {
        "sequence": name,
        "length_nt": length,
        "codons": n_codons,
        "A_pct": counts["A"] / length * 100,
        "C_pct": counts["C"] / length * 100,
        "G_pct": counts["G"] / length * 100,
        "U_pct": counts["T"] / length * 100,
        "GC_pct": (counts["G"] + counts["C"]) / length * 100,
        "GC3_pct": ((third_positions.count("G") + third_positions.count("C"))
                    / len(third_positions) * 100) if third_positions else float("nan"),
        "CpG_n": cpg_n,
        "CpG_per_kb": cpg_n / length * 1000,
        "CpG_OE": cpg_oe,
        "UpA_n": upa_n,
        "UpA_per_kb": upa_n / length * 1000,
        "UpA_OE": upa_oe,
    }

    if cps_scores is not None:
        codons = [sequence[i:i + 3] for i in range(0, n_codons * 3, 3)]
        pairs = [codons[i] + codons[i + 1] for i in range(len(codons) - 1)]
        if pairs:
            result["CPB"] = sum(cps_scores.get(p, 0.0) for p in pairs) / len(pairs)

    if length % 3:
        print(f"warning: {name} is {length} nt, not a multiple of three; "
              "GC3 and CPB assume the reading frame starts at position 1")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", nargs="+", help="one or more FASTA files")
    parser.add_argument("--cps", help="codon pair score table, e.g. data/CPS_human.csv")
    parser.add_argument("--output", default="composition.csv", help="output CSV path")
    args = parser.parse_args()

    cps_scores = load_cps(args.cps) if args.cps else None

    rows = []
    for path in args.fasta:
        for record in SeqIO.parse(path, "fasta"):
            rows.append(analyse(record.id, str(record.seq), cps_scores))

    if not rows:
        raise SystemExit("no sequences found")

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v)
                             for k, v in row.items()})

    columns = ["sequence", "length_nt", "GC_pct", "GC3_pct", "U_pct",
               "CpG_n", "CpG_OE", "UpA_n", "UpA_OE"]
    if cps_scores is not None:
        columns.append("CPB")
    print("\t".join(columns))
    for row in rows:
        print("\t".join(f"{row[c]:.3f}" if isinstance(row[c], float) else str(row[c])
                        for c in columns))
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
