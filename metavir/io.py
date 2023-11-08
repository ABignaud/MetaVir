#!/usr/bin/env python3
# coding: utf-8

"""Core tools to build I/O for MetaVir

This mdoule contains all core I/O functions:
    - generate_temp_dir
    - import_anvio_binning
    - import_contig_data_phages
    - import_network
    - import_phages_contigs
    - write_phage_data
"""


import networkx as nx
import os
import pandas as pd
from os.path import join, exists
from random import getrandbits


def generate_temp_dir(path):
    """Temporary directory generation

    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Generates a temporary file with a random name at the input path.

    Parameters:
    -----------
    path : str
        The path at which the temporary directory will be created.

    Returns:
    --------
    str
        The path of the newly created temporary directory.
    """
    exist = True
    while exist:
        # Keep trying random directory names if they already exist
        directory = str(hex(getrandbits(32)))[2:]
        full_path = join(path, directory)
        if not exists(full_path):
            exist = False
    try:
        os.makedirs(full_path)
    except PermissionError:
        raise PermissionError(
            "The temporary directory cannot be created in {}. "
            "Make sure you have write permission.".format(path)
        )
    return full_path


def import_anvio_binning(binning_file):
    """Import Anvio binning file.

    Parameters:
    -----------
    binning_file : str
        Path to the binning file from anvio to import.

    Returns:
    --------
    dict:
        Dictionnary with contig name as keys and bin name as value.
    """
    binning_result = {}
    with open(binning_file, "r") as binning:
        for line in binning:
            line = line.split()
            # Remove the split name add by anvio
            binning_result[line[0].split("_split")[0]] = line[1]
    return binning_result


def import_contig_data_phages(contig_data_file, binning_result, phages_list):
    """Import contigs data.

    Parameters:
    -----------
    contig_data_file : str
        Path to the contigs data file from MetaTOR.
    binning_result : dict
        Dictionnary with contig name as keys and bin name as value.
    phages_list : list
        List of the phages contigs names.

    Returns:
    --------
    pandas.DataFrame:
        Table with the contig information given with more column with the given
        anvio binning and the phage annotation.
    list
        List of the phage contigs ID.
    """
    contig_data = pd.read_csv(contig_data_file, sep="\t")
    contig_data["Binned"] = False
    contig_data["Final_bin"] = "ND"
    contig_data["MGE"] = False
    phages_list_id = []
    for i in contig_data.index:
        if contig_data.loc[i, "Name"] in phages_list:
            contig_data.loc[i, "MGE"] = True
            phages_list_id.append(contig_data.index[i])
        try:
            contig_data.loc[i, "Final_bin"] = binning_result[
                contig_data.loc[i, "Name"]
            ]
            contig_data.loc[i, "Binned"] = True
        except KeyError:
            continue
    return contig_data, phages_list_id


def import_network(network_file):
    """Import MetaTOR network file.

    Parameters:
    -----------
    network_file : str
        Path to the network file to import.

    Returns:
    --------
    networkx.classes.graph.Graph:
        Network as networkx class.
    """
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )
    return network


def import_phages_contigs(phages_file):
    """Import list of phages contigs.

    Parameters:
    -----------
    phages_file : str
        Path to the phages file which contains the list of phages contigs. One
        contig per line.

    Returns:
    --------
    list:
        List of the phages contigs names.
    """
    phages_list = []
    with open(phages_file, "r") as phages:
        for contig in phages:
            phages_list.append(contig.split()[0])
    return phages_list


def write_phage_data(phage_data, out_file):
    """Write phage binning information.

    Parameters:
    -----------
    phages_data : pandas.DataFrame
        Table with the phage binning information.
    out_file : str
        Path to write the phage data information.
    """

    # Drop the phage id column which is the concatenation of the two previous
    # ones.
    try:
        phage_data.drop("tmp", inplace=True, axis=1)
    except KeyError:
        pass

    # Write the data frame
    phage_data.to_csv(out_file, sep="\t", index=False, float_format="%.2f")
