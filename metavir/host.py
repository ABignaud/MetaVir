#!/usr/bin/env python3
# coding: utf-8

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
        self.binned = False

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
        contig_data : dict
            Dictionnary with the contig name as keys and with the values of the
            bin names, and either if the contig is binned and if it's a phage
            contig. The name of the keys are "bin", "binned", "phage".
        """
        self.bins = dict()
        for i in range(self.len):
            node = str(self.nodes[i])
            score = self.score[i]
            # Uncomment to add virus as bins.
            # if contig_data[node]["virus"]:
            #     try:
            #         self.bins[contig_data[node]["name"]]["score"] += score
            #     except KeyError:
            #         self.bins[contig_data[node]["name"]] = {"score" : score}
            if contig_data[node]["binned"] and not contig_data[node]["phage"]:
                try:
                    self.bins[contig_data[node]["bin"]]["score"] += score
                except KeyError:
                    self.bins[contig_data[node]["bin"]] = {"score": score}
        self.binned = True

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
        if not self.binned:
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
    host_data = pd.DataFrame(
        ["-"] * len(phages_list_id), index=phages_list, columns=["host"]
    )
    A, B, C = 0, 0, 0
    for contig_id in phages_list_id:
        network_id = int(contig_id)
        contig = contig_data[contig_id]["name"]
        subnetwork = network.edges(network_id, data="weight")
        sub = Subnetwork(subnetwork)
        sub.setScore()
        sub.setBinScore(contig_data)
        bin_name, score = sub.getMaxBinScore()
        count = sub.getBinScore()
        mask = host_data.index == contig_data[contig_id]["name"]
        if count == 0:
            C += 1
            host_data.loc[mask, "host"] = "No associated bin. ID: " + str(C)
        elif count == 1:
            A += 1
            host_data.loc[mask, "host"] = bin_name
        else:
            B += 1
            host_data.loc[
                mask, "host"
            ] = "More than one associated bin. ID: " + str(B)

    logger.info("{0} phages associated with one bin.".format(A))
    logger.info("{0} phages associated with more than one bin.".format(B))
    logger.info("{0} phages associated with no bin.".format(C))

    # Save the host table.
    host_data.to_csv(outfile, sep="\t")
    return host_data
