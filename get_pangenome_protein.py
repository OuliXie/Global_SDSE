#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm

# Script for translating nucleotide sequence from pangenome to amino acid
# Takes multisequence alignments and/or fasta containing sequences of a single gene
# Translates sequence, excludes sequences with internal stop codon
# then takes longest sequences and writes it to a representative pan-genome aa fasta file

def trim_seq(record):
    # remove gaps and trim sequence to multiple of 3
    no_gap = record.seq.ungap("-")
    remainder = len(no_gap) % 3
    if remainder == 0:
        return SeqRecord(seq = Seq(no_gap), id = record.id, description = "")
    else:
        return SeqRecord(seq = Seq(no_gap[:-remainder]), id = record.id, description = "")

def translate_fasta(fasta):
    # read in fasta
    prefix = os.path.basename(fasta).split(".")[0]
    translated_record = []

    # parse in records and translate sequence to first stop codon
    for record in SeqIO.parse(fasta, "fasta"):
        trimmed = trim_seq(record)
        aa = trimmed.seq.translate(table=11).strip("*")
        aa_record = SeqRecord(seq = aa, id = prefix, description = "")
        # only append record if no internal stop codon (excludes pseudogenes)
        if aa_record.seq.find("*") == -1:
            translated_record.append(aa_record)

    # find longest record
    longest = max(translated_record, key=lambda x: len(x))

    return longest

def main():

    description = "Generates representative pan-genome amino acid fasta"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_directory",
        help="Directory of multisquence alignments or per-gene multifastas",
        required=True,
        type=os.path.abspath,
        nargs="+",
    )

    args = parser.parse_args()

    # get fasta file paths
    directory = args.input_directory[0]
    input_files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    protein_seqs = []

    for fasta in tqdm(input_files):
        protein_seqs.append(translate_fasta(fasta))

    SeqIO.write(protein_seqs, "pan_genome_reference_protein.fa", "fasta")

    return


if __name__ == "__main__":
    main()
