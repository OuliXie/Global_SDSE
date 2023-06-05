#!/usr/bin/env python

import os
import argparse
import pandas as pd
from tqdm import tqdm


def assign_core(summary_df, core):
    # using Panaroo classification of core here which is based on % threshold without rounding
    total = max(summary_df["No. isolates"])

    summary_df["Classification"] = ["Core" if (x / total) >= core else
                                    "Non_MGE" for x in summary_df["No. isolates"]]

    return summary_df


def add_mge(file, summary_df):
    # read in accessory segment
    segment = pd.read_csv(file)

    # find summarised class of MGE for segment
    mge_class = segment.columns[-1]

    # remove hotspot designation
    mge_class = mge_class.split(";Hotspot")[0]

    summary_df.loc[list(segment["Gene"]), mge_class] += 1

    return summary_df


def final_assign(summary_df, core):
    # create temporary new column to find genes which have had at least 1 classification
    summary_df["sum_class"] = summary_df.loc[:, "Phage":"Non_MGE"].sum(axis=1)
    # create temporary new column to count number of MGE-related occurrences
    summary_df["sum_mge"] = summary_df.loc[:, "Phage":"Degraded"].sum(axis=1)
    # create temporary new column to count number of non-degraded MGE-related occurrences
    summary_df["sum_intact"] = summary_df.loc[:, "Phage":"IS;ME;Integron"].sum(axis=1)

    # If sum of any MGE more often than non-MGE, will be classified as MGE
    # Assign max ignoring degraded because want MGE to override degraded
    summary_df.loc[(summary_df.sum_intact > 0) & (summary_df.sum_mge >= summary_df.Non_MGE), "Classification"] = \
        summary_df.loc[(summary_df.sum_intact > 0) & (summary_df.sum_mge >= summary_df.Non_MGE),
                       "Phage": "IS;ME;Integron"].idxmax(axis=1)
    # If only in degraded or non-MGE, take the greater of the two
    summary_df.loc[(summary_df.sum_mge > 0) & (summary_df.sum_intact == 0), "Classification"] = \
        summary_df.loc[(summary_df.sum_mge > 0) & (summary_df.sum_intact == 0),
                       "Degraded": "Non_MGE"].idxmax(axis=1)
    # If no classifications for gene, then because it is only present on accessory only fragments
    # Some of these may be part of small plasmids depending on accessory only fragment size cutoff selected
    total = max(summary_df["No. isolates"])
    threshold = int(total * core)
    summary_df.loc[(summary_df.sum_class == 0) & (summary_df["No. isolates"] < threshold), "Classification"] = \
        "Unclassified"

    return summary_df


def main():
    description = "Summarises if gene is core or has been detected to be part of a MGE"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="List of classified MGE csvs for each sequence between a single core-core pair",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-g",
        "--gpa",
        dest="gpa",
        help="Annotated gene presence absence csv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-c",
        "--core",
        dest="core",
        help="Core threshold",
        required=True,
        type=float,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=False,
        type=os.path.abspath,
    )

    args = parser.parse_args()

    # create list of input csvs paths
    csvs = []
    with open(args.input_file) as f:
        for line in f:
            line = line.strip()
            csvs.append(line)

    # read in and clean annotated_gene_presence_absence.csv
    summary_df = pd.read_csv(args.gpa, low_memory=False)
    summary_df.set_index("Gene", inplace=True)
    summary_df = summary_df.loc[:, ~summary_df.columns.isin(["Avg sequences per isolate", "Avg sequences per isolate",
                                                             "Genome Fragment", "Order within Fragment",
                                                             "Accessory Fragment", "Accessory Order with Fragment",
                                                             "QC", "Min group size nuc", "Max group size nuc",
                                                             "Avg group size nuc", "recombinase", "T4SS", "phage"])]
    summary_df["Phage"], summary_df["ICE"], summary_df["ME"], summary_df["Phage_like"], summary_df["IS"], \
        summary_df["Integron"], summary_df["IS;Phage"], summary_df["IS;Phage_like"], summary_df["IS;ICE"], \
        summary_df["IS;ME"], summary_df["IS;Integron"], summary_df["Phage;Integron"], \
        summary_df["Phage_like;Integron"], summary_df["ICE;Integron"], summary_df["ME;Integron"], \
        summary_df["IS;Phage;Integron"], summary_df["IS;Phage_like;Integron"], summary_df["IS;ICE;Integron"], \
        summary_df["IS;ME;Integron"], summary_df["Degraded"], summary_df["Non_MGE"], summary_df["Unclassified"], \
        summary_df["Classification"] = [0] * 23

    # assign core genes with Non_MGE as placeholder for other genes
    # note using Panaroo classification of core here which is based on % without rounding
    # otherwise will output a greater number of core genes than Panaroo
    # may result in a small number of genes which are classified as core by Corekaburra as non_MGE here
    assign_core(summary_df, args.core)

    # increment MGE class for each gene in a sequence and core-core combination
    print("Summarising gene types")
    for file in tqdm(csvs):
        add_mge(file, summary_df)

    # pick final classification of MGE association for each gene based on frequency of their association with each
    # MGE class
    final_assign(summary_df, args.core)

    # remove temporary column
    summary_df = summary_df.drop(columns=["sum_class", "sum_mge", "sum_intact"])

    if args.output_dir is None:
        summary_df.to_csv("classified_genes.csv")
    else:
        summary_df.to_csv(args.output_dir + "/" + "classified_genes.csv")


if __name__ == "__main__":
    main()
