#!/usr/bin/env python3
# coding: utf-8

from docopt import docopt
import metavir.binning as mtb
import metavir.host as mth
import metavir.io as mio
import os
import shutil


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


class Binning(AbstractCommand):
    """Bin phages contigs.

    usage:
        binning --depth=STR --fasta=STR --host=STR [--no-clean-up]
        [--outdir=STR] [--threads=1] [--tmpdir=STR]

    options:
        -d, --depth=STR     Path to the depth file from metabat2 script:
                            jgi_summarize_bam_contig_depths.
        -f, --fasta=STR     Path to the fasta file with tha phage contigs
                            sequences.
        -h, --host=STR      Path to the bacterial host associated to the phages
                            contigs generated with metavir host.
        -N, --no-clean-up   If enabled, remove the temporary files.
        -o, --outdir=STR    Path to the output directory where the output will
                            will be written. Default current directory.
        -t, --threads=INT   Number of threads to use for checkV. [Default: 1]
        -T, --tmpdir=STR       Path to temporary directory. [Default: ./tmp]
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

        # Run the phages binning
        mtb.phage_binning(
            self.args["--depth"],
            self.args["--fasta"],
            self.args["--host"],
            self.args["--outdir"],
            remove_tmp,
            int(self.args["--threads"]),
            tmp_dir,
        )

        # Delete pyfastx index:
        os.remove(self.args["--fasta"] + ".fxi")
        # Delete the temporary folder.
        if remove_tmp:
            shutil.rmtree(tmp_dir)
