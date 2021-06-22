#!/usr/bin/env python3
# coding: utf-8

from docopt import docopt
import metavir.host as mth
import metavir.io as mio


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

    usage:
        host --network=STR --binning=STR --phages=STR --contig-data=STR
        [--outfile=STR]

    options:
        -b, --binning=STR       Path to the anvio binning file.
        -c, --contig-data=STR   Path to the MetaTOR contig data file.
        -n, --network=STR       Path to the network file.
        -o, --outfile=STR       Path where to write the output file.
        -p, --phages=STR        Path to the file with phages contigs list.
    """

    def execute(self):
        # Defined the output file if none are given
        if not self.args["--outfile"]:
            self.args["--outfile"] = "./phages_host.tsv"

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
        )
