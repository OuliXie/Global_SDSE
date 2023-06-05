#!/usr/bin/env python

import os
import argparse
import pandas as pd
import numpy as np

# Script for use within per_pair_accessory.sh
# Classifies MGE type based on rules from Khedkar et al. NAR 2022
# Will also assign classification to accessory genes in the segment
# Automatically resolves segments with multiple recombinases detected


def get_rec(rules):
    rec_IS_list = []
    rec_phage_list = []
    rec_CE_list = []
    rec_IS_phage_list = []
    rec_phage_CE_list = []
    # get IS_Tn only recombinases
    for recombinase, classification in rules.items():
        if classification["IS_Tn"] == 1 and classification["Phage"] == 0 and classification["CE"] == 0:
            rec_IS_list.append(recombinase)
        # get phage only recombinases
        elif classification["Phage"] == 1 and classification["IS_Tn"] == 0 and classification["CE"] == 0:
            rec_phage_list.append(recombinase)
        # get CE only recombinases
        elif classification["CE"] == 1 and classification["IS_Tn"] == 0 and classification["Phage"] == 0:
            rec_CE_list.append(recombinase)
        # get IS_Tn|phage recombinases
        elif classification["IS_Tn"] == 1 and classification["Phage"] == 1 and classification["CE"] == 0:
            rec_IS_phage_list.append(recombinase)
        # get phage|CE recombinases
        elif classification["Phage"] == 1 and classification["CE"] == 1 and classification["IS_Tn"] == 0:
            rec_phage_CE_list.append(recombinase)

    return rec_IS_list, rec_phage_list, rec_CE_list, rec_IS_phage_list, rec_phage_CE_list


def process_csv(sequence):
    # reads in csv file and adds MGE classification columns
    segment = pd.read_csv(sequence)
    segment = segment.rename(columns={"recombinase": "Recombinase", "phage": "Ph_gene"})

    # add blank columns for MGE types
    segment["IS"], segment["Phage"], segment["Phage_like"], \
        segment["ICE"], segment["ME"], segment["Integron"], segment["Degraded"], \
        segment["Hotspot"] = [0, 0, 0, 0, 0, 0, 0, 0]

    return segment


def classify_single(r, rules, phage_genes, T4SS_genes, mge_df):
    # find recombinase category
    for category, recombinases in rules.items():
        for i in recombinases:
            if r == i:
                rec_gene = category

    # debugging
    try:
        rec_gene
    except NameError:
        print(r)
        print(rules)
        exit(1)

    # increment columns for recombinases belonging to only 1 category of MGE
    if rec_gene == "IS_tn":
        mge_df["IS"] += 1
    elif rec_gene == "integron":
        mge_df["Integron"] += 1
    elif rec_gene == "phage":
        # must have 2 or more phage structural genes in segment to be called phage
        if phage_genes >= 2:
            mge_df["Phage"] += 1
            # if both ICE genes also present, add hotspot designation
            if any(x in ["virb4", "I_traU"] for x in T4SS_genes) and any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes):
                mge_df["Hotspot"] += 1
        else:
            mge_df["Phage_like"] += 1
            # if both ICE genes also present, add hotspot designation
            if any(x in ["virb4", "I_traU"] for x in T4SS_genes) and any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes):
                mge_df["Hotspot"] += 1
    elif rec_gene == "CE":
        # must have at least 1 ATPas and 1 coupling protein within segment to be called ICE
        if any(x in ["virb4", "I_traU"] for x in T4SS_genes) and any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes):
            mge_df["ICE"] += 1
            # if called ICE but >=2 phage genes also present, then label as hotspot as possible nested element
            if phage_genes >= 2:
                mge_df["Hotspot"] += 1
        else:
            mge_df["ME"] += 1
            # if called ME but >=2 phage genes also present, then label as hotspot as possible nested element
            if phage_genes >= 2:
                mge_df["Hotspot"] += 1
    elif rec_gene == "IS_tn_phage":
        if phage_genes >= 2:
            mge_df["Phage"] += 1
            # if both ICE genes also present, add hotspot designation
            if any(x in ["virb4", "I_traU"] for x in T4SS_genes) and any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes):
                mge_df["Hotspot"] += 1
        else:
            mge_df["IS"] += 1
    elif rec_gene == "IS_tn_CE":
        if any(x in ["virb4", "I_traU"] for x in T4SS_genes) and any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes):
            mge_df["ICE"] += 1
            # if called ICE but >=2 phage genes also present, then label as hotspot as possible nested element
            if phage_genes >= 2:
                mge_df["Hotspot"] += 1
        else:
            mge_df["IS"] += 1
    elif rec_gene == "phage_CE":
        # for recombinases which can be found in phage or CE
        # if at least 2 phage structural genes and does not contain both essential ICE genes, then call phage
        if phage_genes >= 2 and \
                (False == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                           any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            mge_df["Phage"] += 1
        # if at least 2 phage structural genes AND contains both essential ICE genes, then can't classify and call ME
        # and add hotspot designation
        elif phage_genes >= 2 and \
                (True == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                          any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            mge_df["ME"] += 1
            mge_df["Hotspot"] += 1
        # if less than 2 phage structural genes AND contains both essential ICE genes, then classify as ICE
        elif phage_genes < 2 and \
                (True == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                          any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            mge_df["ICE"] += 1
        # if less than 2 phage structural genes and doesn't contain both essential ICE genes, then classify as ME
        elif phage_genes < 2 and \
                (False == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                           any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            mge_df["ME"] += 1

    return mge_df


def classify(segment, recombinase_dict):
    mge_df = segment

    # looks for recombinase, phage/ICE genes, and classifies elements
    rec_genes = segment["Recombinase"].dropna().to_list()
    phage_genes = int(np.sum(segment["Ph_gene"]))
    T4SS_genes = segment["T4SS"].dropna().to_list()

    # Step 1: look if recombinase gene present
    if len(rec_genes) == 0:
        # next look if there are any phage/ICE genes present in segment:
        if phage_genes > 0 or len(T4SS_genes) > 0:
            mge_count = {"IS": 0, "Phage": 0, "Phage_like": 0, "ICE": 0, "ME": 0,
                         "Integron": 0, "Degraded": 1, "Hotspot": 0}
            mge_df["Degraded"] = 1
            return mge_df, mge_count
        else:
            mge_count = {"IS": 0, "Phage": 0, "Phage_like": 0, "ICE": 0, "ME": 0,
                         "Integron": 0, "Degraded": 0, "Hotspot": 0}
            return mge_df, mge_count
    # Step 2: check if single recombinase present in segment
    elif len(rec_genes) == 1:
        mge_df = classify_single(rec_genes[0], recombinase_dict, phage_genes, T4SS_genes, mge_df)
    # Step 3: classify if multiple recombinases present
    elif len(rec_genes) >= 2:
        for r in rec_genes:
            classify_single(r, recombinase_dict, phage_genes, T4SS_genes, mge_df)
        # Add hotspot designation if >= 4 recombinases in the segment as higher likelihood of nested element
        if len(rec_genes) >= 4:
            mge_df["Hotspot"] += 1

    # Step 4: resolve MGE type for segments with multiple recombinase genes
    # set maximum of phage, phage_like, ICE and ME categories to 1
    # allow multiple IS and integrons in nested elements
    mge_count = {"IS": mge_df["IS"].unique()[0], "Phage": mge_df["Phage"].unique()[0],
                 "Phage_like": mge_df["Phage_like"].unique()[0], "ICE": mge_df["ICE"].unique()[0],
                 "ME": mge_df["ME"].unique()[0], "Integron": mge_df["Integron"].unique()[0],
                 "Hotspot": mge_df["Hotspot"].unique()[0]}
    if mge_df["Phage"].unique()[0] >= 1:
        mge_count["Phage"] = 1
        mge_df["Phage"] = 1
    if mge_df["Phage_like"].unique()[0] >= 1:
        mge_count["Phage_like"] = 1
        mge_df["Phage_like"] = 1
    if mge_df["ICE"].unique()[0] >= 1:
        mge_count["ICE"] = 1
        mge_df["ICE"] = 1
    if mge_df["ME"].unique()[0] >= 1:
        mge_count["ME"] = 1
        mge_df["ME"] = 1
    if mge_df["Hotspot"].unique()[0] >= 1:
        mge_count["Hotspot"] = 1
        mge_df["Hotspot"] = 1

    if len(rec_genes) > 1:
        if mge_count["ME"] >= 1 and phage_genes >= 2 and \
                (True == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                          any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            # clear phage/phage_like/ICE classifications as this only occurs if phage and ICE structural genes present
            # and a recombinase which occurs in both phage and CE is present - therefore not classifiable beyond ME
            mge_df["Phage"], mge_df["Phage_like"], mge_df["ICE"] = [0, 0, 0]
            # add hotspot designation (although this should have already been triggered during individual rules)
            mge_df["Hotspot"] = 1
            mge_count.update({"Phage": 0, "Phage_like": 0, "ICE": 0, "Hotspot": 1})
        elif mge_count["Phage"] >= 1 and mge_count["ICE"] >= 1:
            # clear phage/phage_like/ICE classifications and reclassify as ME
            # this only occurs if phage and ICE structural genes present and separate phage-specific
            # and CE-specific recombinases are present - cannot classify beyond ME
            mge_df["Phage"], mge_df["Phage_like"], mge_df["ICE"] = [0, 0, 0]
            mge_count.update({"Phage": 0, "Phage_like": 0, "ICE": 0})
            # add hotspot designation (although this should have already been triggered during individual rules)
            mge_df["ME"], mge_df["Hotspot"] = [1, 1]
            mge_count.update({"ME": 1, "Hotspot": 1})
        elif mge_count["Phage"] >= 1 and mge_count["ME"] >= 1 and \
                (False == (any(x in ["virb4", "I_traU"] for x in T4SS_genes) and
                           any(x in ["t4cp1", "t4cp2"] for x in T4SS_genes))):
            # clear ME and classify as Phage
            mge_df["ME"] = 0
            mge_count.update({"ME": 0})
            # add hotspot designation (although this should have already been triggered during individual rules)
            mge_df["Hotspot"] = 1
            mge_count.update({"Hotspot": 1})
        elif mge_count["Phage_like"] >= 1 and mge_count["ICE"] >= 1:
            # clear Phage_like and classify as ICE
            mge_df["Phage_like"] = 0
            mge_count.update({"Phage_like": 0})
            # add hotspot designation (although this should have already been triggered during individual rules)
            mge_df["Hotspot"] = 1
            mge_count.update({"Hotspot": 1})
        elif mge_count["Phage_like"] >= 1 and mge_count["ME"] >= 1:
            # clear Phage_like and classify as ME as occurs when both phage and CE recombinases present and no
            # phage or ICE structural genes
            mge_df["Phage_like"] = 0
            mge_count.update({"Phage_like": 0})
            # add hotspot designation
            mge_df["Hotspot"] = 1
            mge_count.update({"Hotspot": 1})

    # data check
    if (mge_count["Phage"] + mge_count["Phage_like"] + mge_count["ICE"] + mge_count["ME"]) > 1:
        print("Error in MGE counts - multiple phage/phage-like/ICE/ME counted")
        exit(1)

    return mge_df, mge_count


def summarise_mge(mge_df, mge_count, ignore):
    # gives a summarised classification of the MGE segment
    # e.g., if phage present in segment with 1 or more IS, will classify as IS;Phage

    mge_df_summarised = mge_df.drop(["IS", "Phage", "Phage_like", "ICE", "ME", "Integron", "Degraded", "Hotspot"],
                                    axis=1)
    summarised_class = []
    for k, v in mge_count.items():
        if v >= 1:
            summarised_class.append(k)

    if len(summarised_class) >= 1:
        summarised_class = ';'.join(summarised_class)
        mge_df_summarised[summarised_class] = 1
    else:
        if ignore:
            # Gives "unclassified" category intended for when sequence break in segment and no MGE detected
            # So don't classify accessory genes which are likely part of a MGE but are separated from recombinase
            mge_df_summarised["Unclassified"] = 1
        else:
            mge_df_summarised["Non_MGE"] = 1

    return mge_df_summarised


def main():
    description = "Classifies MGE element from accessory segment"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="Accessory fragment csvs extracted from gene_presence_absence",
        required=True,
        type=os.path.abspath,
        nargs="+"
    )

    io_opts.add_argument(
        "-r",
        "--rules",
        dest="rules",
        help="Path to recombinase rules tsv",
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
    )

    parser.add_argument(
        "-b",
        action="store_true",
        help="Ignore gene classification if no MGE found - use only when sequence break in gene pair",
    )

    args = parser.parse_args()

    # read in recombinase rules file as dictionary
    rules = pd.read_csv(args.rules, sep='\t')
    rules = rules.set_index("rec_subfamily").to_dict("index")

    # list recombinases belonging to only one class of MGE
    recombinase_dict = {}
    recombinase_dict["IS_tn"], recombinase_dict["phage"], recombinase_dict["CE"], \
        recombinase_dict["IS_tn_phage"], recombinase_dict["phage_CE"] = get_rec(rules)
    recombinase_dict["integron"] = ["Integron"]
    recombinase_dict["IS_tn_CE"] = ["c2_n1ser"]
    recombinase_dict["Cellular"] = ["Arch1", "Candidate", "Cyan", "Xer"]

    segments = []
    # read in each segment
    for file in args.input_file:
        segments.append(process_csv(file))

    # classify MGE for segment including resolving multiple recombinases
    mge_processed = [classify(x, recombinase_dict) for x in segments]

    # collapse categories to a single nested classification e.g., IS;Phage
    summarised_mge_df = []
    for seg in mge_processed:
        summarised_mge_df.append(summarise_mge(seg[0], seg[1], args.b))

    # output
    if args.output_dir is None:
        for seg in mge_processed:
            gff = seg[0].columns[8]
            seg[0].to_csv(gff + "_mge_full.csv", index=False)
        for df in summarised_mge_df:
            gff = df.columns[8]
            df.to_csv(gff + "_mge_summarised.csv", index=False)
    else:
        for seg in mge_processed:
            gff = seg[0].columns[8]
            seg[0].to_csv(args.output_dir + "/" + gff + "_mge_full.csv", index=False)
        for df in summarised_mge_df:
            gff = df.columns[8]
            df.to_csv(args.output_dir + "/" + gff + "_mge_summarised.csv", index=False)


if __name__ == "__main__":
    main()
