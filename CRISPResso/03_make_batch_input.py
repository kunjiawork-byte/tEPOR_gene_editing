#!/usr/bin/env python3
"""Build the original CRISPRessoBatch table with an explicit target/read mode."""
import argparse
from pathlib import Path
import csv

def read_fasta(path):
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif line.strip():
            records[name] += line.strip()
    return records

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-dir", required=True, type=Path)
parser.add_argument("--target", required=True, choices=["HBG_CBE", "EPOR_sg1", "EPOR_sg2", "BCL11A_CBE", "BCL11A_Cas9"])
parser.add_argument("--reads", required=True, choices=["r1", "paired"])
parser.add_argument("--output", type=Path, default=Path("CRISPResso_batch_input.txt"))
args = parser.parse_args()
resources = Path(__file__).resolve().parent
amplicons = read_fasta(resources / "amplicons.fasta")
guides = read_fasta(resources / "guides.fasta")
targets = {
    "HBG_CBE": ("HBG2", "sgHBG-1"),
    "EPOR_sg1": ("EPOR", "EPOR-sg1"),
    "EPOR_sg2": ("EPOR", "EPOR-sg2"),
    "BCL11A_CBE": ("BCL11A", "BCL11A-1620 (BCL11A CBE)"),
    "BCL11A_Cas9": ("BCL11A", "BCL11A-1617 (BCL11A Cas9)"),
}
amplicon, guide = targets[args.target]
r1_files = sorted(args.input_dir.resolve().glob("trimmed-*.R1.fastq.gz"))
if not r1_files:
    parser.error("No trimmed-*.R1.fastq.gz files found")
if args.output.resolve() in [resources / "amplicons.fasta", resources / "guides.fasta"]:
    parser.error("Output must not overwrite a sequence resource")
columns = ["name", "amplicon_seq", "guide_seq", "fastq_r1"]
if args.reads == "paired":
    columns.append("fastq_r2")
count = 0
with args.output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(columns)
    for r1 in r1_files:
        base = r1.name[:-len(".R1.fastq.gz")]
        r2 = r1.with_name(base + ".R2.fastq.gz")
        if args.reads == "paired" and not r2.is_file():
            print(f"Skipping {r1.name}: missing {r2.name}")
            continue
        row = [base[len("trimmed-"):], amplicons[amplicon], guides[guide], str(r1)]
        if args.reads == "paired":
            row.append(str(r2))
        writer.writerow(row)
        count += 1
print(f"Wrote {count} samples to {args.output}")
