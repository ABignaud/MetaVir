#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for MetaVir

This module contains all classes related to MetaVir commands:
    - host
    - binning

Note
----
Structure based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
: https://github.com/rgreinho/docopt-subcommands-example
abignaud, 20201118

Raises
------
NotImplementedError
    Will be raised if AbstractCommand is called for some reason instead of one
    of its children.
"""

from docopt import docopt
from os.path import join
import metavir.binning as mtb
import metavir.host as mth
import metavir.io as mio
import os
import shutil
from metavir.log import logger


class AbstractCommand:
    """Abstract base command class

    Base class for the commands from which other metaVir commands derive.
    """

    def __init__(self, command_args, global_args):
        """Initialize the commands"""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError

    def check_output_path(self, path, force=False):
        """Throws error if the output file exists. Create required file tree
        otherwise.
        """
        # Get complete output filename and prevent overwriting unless force is
        # enabled
        if not force and os.path.exists(path):
            raise IOError(
                "Output file already exists. Use --force to overwrite"
            )
        if dirname(path):
            os.makedirs(dirname(path), exist_ok=True)


class Host(AbstractCommand):
    """Detect host of phage annotated contigs.

    It will return an output file with the phages information from MetaTOR
    binning with the added information from the anvio binning and the detected
    bacterial MAG host of the phages.

    usage:
        host --network=FILE --binning=FILE --phages=FILE --contig-data=FILE
        [--outfile=FILE] [--threshold=0.1]

    options:
        -b, --binning=FILE      Path to the anvio binning file.
        -c, --contig-data=FILE  Path to the MetaTOR contig data file.
        -n, --network=FILE      Path to the network file.
        -o, --outfile=FILE      Path where to write the output file.
        -p, --phages=FILE       Path to the file with phages contigs list.
        -t, --threshold=FLOAT   Threshold to consider an association with a MAG.
                                [Default: 0.1]
    """

    def execute(self):
        # Defined the output file if none are given
        if not self.args["--outfile"]:
            self.args["--outfile"] = "./phages_data_host.tsv"

        # Import the files
        binning_result = mio.import_anvio_binning(self.args["--binning"])
        phages_list = mio.import_phages_contigs(self.args["--phages"])
        contig_data, phages_list_id = mio.import_contig_data_phages(
            self.args["--contig-data"], binning_result, phages_list
        )
        network = mio.import_network(self.args["--network"])

        # Run the host detection
        mth.host_detection(
            network,
            contig_data,
            phages_list,
            phages_list_id,
            self.args["--outfile"],
            self.args["--threshold"],
        )


class Binning(AbstractCommand):
    """Bin phages contigs.

    Bin the phage contigs in phages MAGs using the bacterial host detection from
    MetaVir and the metagenomic binning base on sequences and coverage from
    metabat2.

    The results are then checked using checkV and some plots are displayed to
    viusalize them. It will return the updated phages data, the phages fasta,
    the detailed ouputs of checKV and the plots.

    The phage fasta will contain one entry by phage MAG, with 180 "N" spacers
    between contigs.

    usage:
        binning --network=FILE --binning=FILE --phages=FILE --contigs-data=FILE --fasta=FILE
        [--checkv-db=DIR] [--depth=FILE] [--method=pairs] [--no-clean-up]
        [--outdir=DIR] [--pairs=STR] [--plot] [--random] [--threads=1]
        [--tmpdir=DIR] [--threshold-bin=0.8] [--threshold-asso=0.1]

    options:
        -b, --binning=FILE      Path to the anvio binning file.
        -c, --contigs-data=FILE  Path to the MetaTOR contig data file.
        --checkv-db=DIR         Directory where the checkV database is stored.
                                By default the CHECKVDB environment variable is
                                used.
        -d, --depth=FILE        Path to the depth file from metabat2 script:
                                jgi_summarize_bam_contig_depths.
        -f, --fasta=FILE        Path to the fasta file with tha phage contigs
                                sequences.
        -m, --method=STR        Method for the binning. Either 'metabat' or
                                'pairs' [Default: pairs].
        -n, --network=FILE      Path to the network file.
        -N, --no-clean-up       If enabled, remove the temporary files.
        -o, --outdir=DIR        Path to the output directory where the output
                                will be written. Default current directory.
        -q, --pairs=STR         Path of the pairs file separated by a comma.
        -p, --phages=FILE       Path to the file with phages contigs list.
        -P, --plot              If enable, make summary plots.
        -r, --random            If enable, make a random binning.
        -s, --threshold-bin=FLOAT       Threshold to use for binning. 
                                        [Default: 0.8]
        -S, --threshold-asso=FLOAT      Threshold to use for association. If 
                                several MAGs have value higher than this ratio 
                                of total contatcs several association are 
                                considered. [Default: 0.1]
        -t, --threads=INT       Number of threads to use for checkV.
                                [Default: 1]
        -T, --tmpdir=DIR        Path to temporary directory. [Default: ./tmp]
    """

    def execute(self):
        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            self.args["--tmpdir"] = "./tmp"
        tmp_dir = mio.generate_temp_dir(self.args["--tmpdir"])
        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)

        # Set remove tmp for checkV.
        if not self.args["--no-clean-up"]:
            remove_tmp = True
        else:
            remove_tmp = False

        # Set checkV database path
        if not self.args["--checkv-db"]:
            self.args["--checkv-db"] = os.getenv("CHECKVD")

        # Sanity check
        if self.args["--method"] == "pairs" and not self.args["--pairs"]:
            logger.error("Pair file is necessary if method is pairs.")
            raise ValueError

        if self.args["--method"] == "metabat" and not self.args["--depth"]:
            logger.error("Depth file is necessary if method is metabat.")
            raise ValueError

        pairs_files = self.args["--pairs"]
        if pairs_files:
            pairs_files = pairs_files.split(",")

        # Import the files
        binning_result = mio.import_anvio_binning(self.args["--binning"])
        phages_list = mio.import_phages_contigs(self.args["--phages"])
        contigs_data, phages_list_id = mio.import_contig_data_phages(
            self.args["--contigs-data"], binning_result, phages_list
        )
        network = mio.import_network(self.args["--network"])

        # Run the phages binning
        mtb.phage_binning(
            checkv_db=self.args["--checkv-db"],
            depth_file=self.args["--depth"],
            fasta_phages_contigs=self.args["--fasta"],
            network=network,
            contigs_data=contigs_data,
            phages_list_id=phages_list_id,
            out_dir=self.args["--outdir"],
            pairs_files=pairs_files,
            tmp_dir=tmp_dir,
            threshold_bin=float(self.args["--threshold-bin"]),
            threshold_asso=float(self.args["--threshold-asso"]),
            association=True,
            plot=self.args["--plot"],
            remove_tmp=remove_tmp,
            threads=int(self.args["--threads"]),
            method=self.args["--method"],
            random=self.args["--random"],
        )

        # Delete the temporary folder.
        if remove_tmp:
            shutil.rmtree(tmp_dir)
            # Delete pyfastx index:
            os.remove(self.args["--fasta"] + ".fxi")
