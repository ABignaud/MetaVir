#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# abignaud, 20210621

"""
MetaVir pipeline for generating and manipulating phage MAGs.

usage:
    metavir [-hv] <command> [<args>...]
    
options:
    -h, --help              shows the help
    -v, --version           shows the version

The subcommands are:
    host        Detect bacterial host from a metaHiC network binned by metaTOR
                given a annotated phages list.
    binning     Build phages MAGs based on metagenomic binning using metabat2
                and the host detection from the metaHiC data.
    pipeline    Use both others command and run them sequentially.
"""

from docopt import docopt
from docopt import DocoptExit
import metavir.commands as commands
from metavir.version import __version__


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
