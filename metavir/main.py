
#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# abignaud, 20210621

"""
MetaVir pipeline for generating and manipulating phage MAGs.
usage:
    metator [-hv] <command> [<args>...]
options:
    -h, --help              shows the help
    -v, --version           shows the version

The subcommands are:
    host     Generate metaHiC contigs network from fastq reads or bam files
                and normalize it.
    binning   Partition metaHiC network using Louvain or Leiden algorithm.
    pipeline    Use all the others command to give a binned output in one
                command line.
    validation  Validates bins using CheckM and make a recursive partition to
                try to decontaminate them.
"""

from docopt import docopt
from docopt import DocoptExit
import metator.commands as commands
from metator.version import __version__


def main():
    args = docopt(__doc__, version=__version__, options_first=True)
    # Retrieve the command to execute.
    command_name = args.pop("<command>").capitalize()

    # Retrieve the command arguments.
    command_args = args.pop("<args>")
    if command_args is None:
        command_args = {}
    # After 'popping' '<command>' and '<args>', what is left in the args
    # dictionary are the global arguments.

    # Retrieve the class from the 'commands' module.
    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command.")
        raise DocoptExit()
    # Create an instance of the command.
    command = command_class(command_args, args)
    # Execute the command.
    command.execute()


if __name__ == "__main__":
    main()