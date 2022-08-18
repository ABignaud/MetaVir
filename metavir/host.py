#!/usr/bin/env python3
# coding: utf-8

"""Module with the phages bacterial host detection functions. 

It detects phages bacterial host MAGs from the output binning of MetaTOR and the
metaHiC network. The module extract the subnetwork and then identify the bins 
which interaction with the pahge to rank potential bacterial host.

A core class Subentwork is build to handle the phage contigs interaction.

Main function to call the class and build the ouput is host_detection.
"""


import networkx
import pandas as pd
from metavir.log import logger


class Subnetwork:
    """Class to handle subnetwork information of one bin."""

    def __init__(self, subnetwork):
        """Initialize nodes list and their weights (number of hits with the
        associated contigs normalized by metaTOR).

        Parameters:
        -----------
        subnetwork : networkx.classes.reportviews.EdgeDataView
            Subnetwork with all edges of a given contig.
        """
        self.nodes = []
        self.weights = []
        for edge in subnetwork:
            self.id = edge[0]
            self.nodes.append(edge[1])
            self.weights.append(edge[2])
        self.len = len(self.weights)
        self.scored = False

    def setScore(self):
        """Set scores for each associated contigs. The score is just a ratio of
        the sum of the weigths.
        """
        self.score = [x / sum(self.weights) for x in self.weights]

    def getMaxScore(self):
        """Return highest score and the associated contig.

        Returns:
        --------
        str:
            Contig name with the highest score (maximum connectivity with the
            given contig).
        float:
            Highest score value (contig score).
        """
        return self.nodes[self.score.index(max(self.score))], max(self.score)

    def setBinScore(self, contig_data):
        """Set scores for each associated bins.

        Parameters:
        -----------
        contig_data : pandas.DataFrame:
            Table with the contig information given with more column with the given
            anvio binning and the phage annotation.
        """
        self.bins = dict()
        for i in range(self.len):
            node = self.nodes[i] - 1
            score = self.score[i]
            # Uncomment to add virus as bins.
            # if contig_data[node]["virus"]:
            #     try:
            #         self.bins[contig_data[node]["name"]]["score"] += score
            #     except KeyError:
            #         self.bins[contig_data[node]["name"]] = {"score" : score}
            if (
                contig_data.loc[node, "Binned"]
                and not contig_data.loc[node, "Phage"]
            ):
                try:
                    self.bins[contig_data.loc[node, "Final_bin"]][
                        "score"
                    ] += score
                except KeyError:
                    self.bins[contig_data.loc[node, "Final_bin"]] = {
                        "score": score
                    }
        self.scored = True

    def getMaxBinScore(self, contig_data=None):
        """Return highest score and the associated bin.

        Returns:
        --------
        str:
            Bin name with the highest score (maximum connectivity with the given
            contig).
        float:
            Highest score value (bin score).
        """
        if not self.scored:
            self.setBinScore(contig_data)
        max_score = 0
        max_bin = "-"
        for bin_name in self.bins:
            score = self.bins[bin_name]["score"]
            if score > max_score:
                max_score = score
                max_bin = bin_name
        return max_bin, max_score

    def getBinScore(self):
        """Return the count of connected bin. The threshold is set to 0.2 to be
        considered as connected.

        Returns:
        --------
        int:
            Count of connected bins.
        """
        c = 0
        for bin_name in self.bins:
            score = self.bins[bin_name]["score"]
            if score >= 0.2:
                c += 1
        return c


def host_detection(network, contig_data, phages_list, phages_list_id, outfile):
    """Main function to detect the host.

    Parameters:
    -----------
    network : networkx.classes.graph.Graph
        MetaTOR network of the HiC data.
    contig_data : dict
        Dictionnary with the contig name as keys and with the values of the
        contig id the associated bin name, and either if the contig is binned
        and if it's a phage contig. The name of the keys are "id", "bin",
        "binned", "phage".
    phages_list : list
        List of the phages contigs names.
    phages_list_id : list
        List of the phages contigs IDs.
    outfile : str
        Path  to write the ouput table.

    Returns:
    --------
    dict:
        Dictionnary with the phage contig as keys and with the associated host
        as value.
    """
    # Compute the score with the subnetwork and return bins in each categories
    # and build the associated table.
    phage_data = pd.DataFrame(contig_data.loc[phages_list_id, :])
    phage_data["Host"] = "ND"
    A, B, C = 0, 0, 0
    for contig_id in phages_list_id:
        network_id = contig_id + 1
        contig = contig_data.loc[contig_id, "Name"]
        # Do not compute subnetwork with no HiC contacts.
        if contig_data.loc[contig_id, "Hit"] > 0:
            # Manage the case of the self interacting contigs.
            try:
                subnetwork = network.edges(network_id, data="weight")
                sub = Subnetwork(subnetwork)
                sub.setScore()
                sub.setBinScore(contig_data)
                bin_name, score = sub.getMaxBinScore()
                count = sub.getBinScore()
                if count == 0:
                    C += 1
                    phage_data.loc[contig_id, "Host"] = "No associated bin. ID: " + str(
                        C
                    )
                elif count == 1:
                    A += 1
                    phage_data.loc[contig_id, "Host"] = bin_name
                else:
                    B += 1
                    phage_data.loc[
                        contig_id, "host"
                    ] = "More than one associated bin. ID: " + str(B)
            except networkx.exception.NetworkXError:
                C += 1
                phage_data.loc[contig_id, "Host"] = "No associated bin. ID: " + str(
                    C
                )
        else:
            C += 1
            phage_data.loc[contig_id, "Host"] = "No associated bin. ID: " + str(
                C
            )

    logger.info("{0} phages associated with one bin.".format(A))
    logger.info("{0} phages associated with more than one bin.".format(B))
    logger.info("{0} phages associated with no bin.".format(C))

    # Save the host table.
    phage_data.to_csv(outfile, sep="\t")
    return phage_data
