#!/usr/bin/env python3
# coding: utf-8

"""Module to plot phages binning result from checkv output summary.

Core function to plot:
    - build_vmags_summary
    - pie_bins_size_distribution
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def build_vmags_summary(checkv_summary):
    """Build vMAGs quality summary table.

    Parameters:
    -----------
    checkv_summary : pandas.core.frame.DataFrame
        Table with the informations about the final phage bins.
    """
    # Add a qualitive quality column in bin_summary with the quality of the MAGs
    # checkv_summary["MAG_quality"] = "ND"
    # Change NA value to 0
    mask = checkv_summary.completeness == "NA"
    checkv_summary.loc[mask, "completeness"] = 0
    checkv_summary = checkv_summary.loc[checkv_summary.provirus == "No", :]
    checkv_summary = checkv_summary.reset_index()
    # Build a small table with the sum for each quality category.
    mags_summary = pd.DataFrame(
        [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]],
        columns=["bins", "size"],
    )
    for i in range(len(checkv_summary)):
        completness = float(checkv_summary.loc[i, "completeness"])
        size = int(checkv_summary.loc[i, "contig_length"])
        if completness >= 90:
            mags_summary.loc[0, "bins"] += 1
            mags_summary.loc[0, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">90"
        if completness >= 70:
            mags_summary.loc[1, "bins"] += 1
            mags_summary.loc[1, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">70"
        elif completness >= 50:
            mags_summary.loc[2, "bins"] += 1
            mags_summary.loc[2, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">50"
        elif completness >= 30:
            mags_summary.loc[3, "bins"] += 1
            mags_summary.loc[3, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">30"
        else:
            mags_summary.loc[4, "bins"] += 1
            mags_summary.loc[4, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = "<30"
    return mags_summary


def pie_bins_size_distribution(checkv_summary, out_file):
    """Function to plot a camembert of the fraction of size corresponding
    to their completness.

    Parameters:
    -----------
    checkv_summary : pandas.core.frame.DataFrame
        Table with the informations about the final phage bins.
    out_file : str
        Path where to save the figure.
    """
    mags_summary = build_vmags_summary(checkv_summary)
    # Plot the camembert of the size ratio.
    labels = [
        ">90%: {0} - {1}Mb".format(
            int(mags_summary.loc[0, "bins"]),
            round(mags_summary.loc[0, "size"] / 1000000, 2),
        ),
        ">70%: {0} - {1}Mb".format(
            int(mags_summary.loc[1, "bins"]),
            round(mags_summary.loc[1, "size"] / 1000000, 2),
        ),
        ">50%: {0} - {1}Mb".format(
            int(mags_summary.loc[2, "bins"]),
            round(mags_summary.loc[2, "size"] / 1000000, 2),
        ),
        ">30%: {0} - {1}Mb".format(
            int(mags_summary.loc[3, "bins"]),
            round(mags_summary.loc[3, "size"] / 1000000, 2),
        ),
        "Others contigs: {0} - {1}Mb".format(
            int(mags_summary.loc[4, "bins"]),
            round(mags_summary.loc[4, "size"] / 1000000, 2),
        ),
    ]
    total_size = np.sum(list(mags_summary["size"]))
    fig, ax = plt.subplots()
    plt.pie(
        mags_summary["size"],
        colors=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.legend(labels, bbox_to_anchor=(0.9, 0.0, 1.0, 1.0), loc="upper right")
    plt.text(
        -1.5,
        -1.2,
        "Total size of the assembly: {0}Mb".format(
            round(total_size / 1000000, 2)
        ),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.35,
        "Percentage of phage 50% complete: {0}%".format(
            round(sum(mags_summary["size"][0:3]) / total_size * 100, 2)
        ),
        fontdict=None,
    )
    plt.title("Size proportion of phage bins depending on their completness.")
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")


def barplot_bins_size(labels, list_checkv_summary, out_file):
    """Plot the size distribution of the bins quality of different experiments
    or methods.

    Parameters:
    -----------
    labels : list os str
        List of the labels of the different experiments/methods.
    list_checkv_summary : list of pandas.DataFrame
        List of the checkV quality summary tables.
    out_file : str
        Path to the output file.
    """
    # Create an empty dataframe.
    stacked_data = pd.DataFrame(
        [[0] * 5] * len(labels),
        columns=[">90", ">70", ">50", ">30", "<30"],
        index=labels,
    )

    # Extract vMAGs summary and add it to the table.
    for i, checkv_summary in enumerate(list_checkv_summary):
        mags_summary = build_vmags_summary(checkv_summary)
        stacked_data.loc[labels[i], :] = list(mags_summary["size"])

    # Plot the table.
    stacked_data = stacked_data.apply(lambda x: x * 100 / sum(x), axis=1)
    stacked_data.plot(
        kind="bar",
        stacked=True,
        color=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.title("Phages bins quality distribution")
    plt.xlabel("Experiment")
    plt.ylabel("Percentage bins size (%)")
    plt.legend(loc="upper right", bbox_to_anchor=(0.2, 0.0, 1.0, 1.0))
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")


def barplot_bins_number(labels, list_checkv_summary, out_file):
    """Plot the number distribution of the bins quality of different experiments
    or methods.

    Parameters:
    -----------
    labels : list os str
        List of the labels of the different experiments/methods.
    list_checkv_summary : list of pandas.DataFrame
        List of the checkV quality summary tables.
    out_file : str
        Path to the output file.
    """
    # Create an empty dataframe.
    stacked_data = pd.DataFrame(
        [[0] * 5] * len(labels),
        columns=[">90", ">70", ">50", ">30", "<30"],
        index=labels,
    )

    # Extract vMAGs summary and add it to the table.
    for i, checkv_summary in enumerate(list_checkv_summary):
        mags_summary = build_vmags_summary(checkv_summary)
        stacked_data.loc[labels[i], :] = list(mags_summary["bins"])

    # Plot the table.
    stacked_data = stacked_data.apply(lambda x: x * 100 / sum(x), axis=1)
    stacked_data.plot(
        kind="bar",
        stacked=True,
        color=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.title("Phages bins quality distribution")
    plt.xlabel("Experiment")
    plt.ylabel("Percentage bins size (%)")
    plt.legend(loc="upper right", bbox_to_anchor=(0.2, 0.0, 1.0, 1.0))
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
