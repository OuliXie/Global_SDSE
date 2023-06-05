#!/usr/bin/env python

import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

# Script for sorting annotated_gene_presence_absence
# Uses associated gff to fix issues with re found genes and out of order locus tags from PGAP


def parse_gff(file):
    # Adapted from post_run_gff_output.py
    # Read in and split off FASTA portion
    raw_file = open(file, 'r')
    # Read in and split off FASTA portion
    lines = raw_file.read().replace(',', '')
    split = lines.split('##FASTA')[0]
    # Remove headers
    body = []
    for line in split.splitlines():
        if "##" not in line:
            body.append(line)
    # Parse locus tag only in order
    parsed_gff = []
    for gff_line in body:
        initial_split = gff_line.split("\t")
        attribute_split = " ".join(initial_split[8:]).split(";")
        locus_tag = ""
        for attribute in attribute_split:
            if "locus_tag" in attribute:
                locus_tag = attribute.split("=")[1]
        parsed_gff.append(locus_tag)

    return parsed_gff


def order_extract(genome, gpa, gff_path, out_dir):
    # Read in gff
    gff_locus = parse_gff(gff_path + "/" + genome + ".gff")
    # Select columns from annotated_gene_presence_absence.csv
    gpa = gpa[["Gene", "recombinase", "T4SS", "phage", "Non-unique Gene name", "Annotation", "No. isolates",
                       "No. sequences", genome]]
    # Temp column containing only the first locus_tag for CDS which were found to belong to the same gene
    unique_tag = [x.split(";")[0] for x in gpa[genome].to_list()]
    gpa = gpa.assign(locus_unique=unique_tag)
    # Check which locus tags have been merged in gff
    gff_unique = [x for x in gff_locus if x in unique_tag]
    # Subset annotated_gene_presence_absence.csv
    gpa = gpa[gpa.loc[:, "locus_unique"].isin(gff_unique)]
    gpa = gpa.set_index("locus_unique")
    # Reorder based on gff order
    gpa = gpa.reindex(gff_unique)

    # Remove phage annotations if gene is also a recombinase
    gpa["phage"] = np.where(gpa["recombinase"] != "", "", gpa["phage"])

    gpa.to_csv(out_dir + "/" + genome + "_gene_presence_absence.csv", index=False)

    return


def main():
    description = "Orders annotated_gene_presence_absence.csv by each genome"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="annotated_gene_presence_absence.csv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-g",
        "--gffs",
        dest="gffs",
        help="Path to post-Panaroo corrected gffs",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=False,
        type=os.path.abspath,
        default="ordered_gene_presence_absence"
    )

    args = parser.parse_args()

    # Make dedicated directory for outputs
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Read in annotated_gene_presence_absence
    gpa = pd.read_csv(args.input_file, low_memory=False).fillna("")

    # Make copy in output folder
    gpa.to_csv(args.output_dir + "/" + "annotated_gene_presence_absence.csv", index=False)

    # Get list of genomes
    exclude = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
               "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
               "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc",
               "Avg group size nuc", "recombinase", "T4SS", "phage"]
    genomes = gpa.loc[:, ~gpa.columns.isin(exclude)].columns

    for genome in tqdm(genomes):
        order_extract(genome, gpa, args.gffs, args.output_dir)

    # Write sequence_id.txt file
    pd.DataFrame(genomes.to_list()).to_csv("sequence_id.txt", index=False, header=False)


if __name__ == "__main__":
    main()
