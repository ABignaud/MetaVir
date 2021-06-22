#!/usr/bin/env python3
# coding: utf-8

import networkx as nx


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
    dict:
        Dictionnary with the contig id as keys and with the values of the contig
        name the associated bin name, and either if the contig is binned and if
        it's a phage contig. The name of the keys are "id", "bin", "binned",
        "phage".
    list
        List of the phages ID.
    """
    contig_data = dict()
    phages_list_id = []
    with open(contig_data_file, "r") as file:
        for line in file:
            # Skip header.
            if line.startswith("ID"):
                continue
            line = line.split()
            try:
                bin_name = binning_result[line[1]]
                binned = True
            except KeyError:
                binned = False
            if line[1] in phages_list:
                phage = True
                phages_list_id.append(line[0])
            else:
                phage = False
            contig_data[line[0]] = {
                "name": line[1],
                "bin": bin_name,
                "binned": binned,
                "phage": phage,
            }
    return contig_data, phages_list_id
