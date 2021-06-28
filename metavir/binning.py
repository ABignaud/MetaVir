#!/usr/bin/env python3
# coding: utf-8

"""Module with the phages binning functions. 

It generates phages bins from the detected bacterial host of the contigs and 
from the metabat2 clusterization based on the sequence and shotgun coverage 
information. 

Core function to partition phages contigs:
    - build_phage_depth
    - generate_phages_bins
    - generate_phages_fasta
    - run_checkv
    - run_metabat
"""


import checkv
import metavir.figures as mtf
import metavir.io as mio
import pandas as pd
import subprocess as sp
from metavir.log import logger
from os.path import join
import shutil


def build_phage_depth(contigs_file, depth_file, phages_data, phage_depth_file):
    """Build phage depth form the whole assembly depth file from metabat script.

    Parameters:
    -----------
    contigs_file : str
        Path to the temporary file containing the list of the phage contigs name
        in the same order as the depth file will be written.
    depth_file : str
        Path to the wholed depth file from Metabat2 script.
    phages_data : pandas.DataFrame
        Table with the phage contig names as index and the detected bacterial
        bins as column.
    phage_depth_file : str
        Path to write the depth file with only the phage contigs depth.
    """

    # Import the whole depth file as dataframe.
    whole_depth = pd.read_csv(depth_file, sep="\t")

    # Extract contigs name list.
    phage_list = list(phages_data.Name)

    # Extract line of the phage contigs
    mask = []
    for i in whole_depth.contigName:
        if i in phage_list:
            mask.append(True)
        else:
            mask.append(False)
    phage_depth = whole_depth.loc[mask]

    # Write the contigs list as the same order as the depth file.
    with open(contigs_file, "w") as f:
        for contig_name in phage_depth.contigName:
            f.write("%s\n" % contig_name)

    # Write phages depth file to use metabat2.
    phage_depth.to_csv(phage_depth_file, sep="\t", index=False)


def generate_phage_bins(phages_data):
    """Generates the binning of the pahges contigs based on both HiC
    information (host detection) and the coverage and sequences (metabat2
    binning)

    Parameters:
    -----------
    phages_data : pandas.DataFrame
        Table with the contigs name as index and with information from both
        host detected from metavir host and cluster form metabat2.

    Returns:
    --------
    pandas.DataFrame:
        Input table with the phage bin id column added.
    dict:
        Dictionnary with the phage bin id as key and the list of the contigs
        name as value.
    """

    # Creates an unique ID for each future bin.
    phages_data["tmp"] = (
        phages_data.Host + "___" + list(map(str, phages_data.Metabat_bin))
    )

    # Create a new column with the bin id information added.
    bins_ids = {}
    phage_bins = {}
    phages_data["MetaVir_bin"] = 0
    bin_id = 0
    for contig in phages_data.index:
        phage_id = phages_data.loc[contig, "tmp"]
        # Test if the phage id have been already seen.
        try:
            bin_id_old = bins_ids[phage_id]
            phages_data.loc[contig, "MetaVir_bin"] = bin_id_old
            phage_bins[bin_id_old].append(phages_data.loc[contig, "Name"])
        # Increment the bin id if it's the first time the phage id have been
        # seen.
        except KeyError:
            bin_id += 1
            bins_ids[phage_id] = bin_id
            phages_data.loc[contig, "MetaVir_bin"] = bin_id
            phage_bins[bin_id] = [phages_data.loc[contig, "Name"]]
    return phages_data, phage_bins


def generate_phages_fasta(fasta, phage_bins, out_file, tmp_dir):
    """Generate the fasta file with one sequence entry for each bin. The
    sequences are generated like that to be accepted by checkV as one virus. In
    the sequences 180 "N" spacers are added between two contigs.

    Parameters:
    -----------
    fasta : str
        Path to the fasta file with the original phages contigs sequences.
    phage_bins : dict
        A dictionnary with the id of the phage bins as keys and the list
        of name of their contigs as values.
    out_file : str
        Path to the output file where the fasta of all the phages bins will be
        written.
    tmp_dir : str
        Path to the temporary directory to write the temporary contigs list
        files.
    """

    nb_bins = 0
    # For each bin create a list of the contigs and extract them from the
    # fasta to create a new fasta file with only the bin.
    with open(out_file, "w") as out:
        for bin_id in phage_bins:
            # Extract the list of the contigs from the contigs data file.
            list_contigs_name = phage_bins[bin_id]
            nb_bins += 1
            # Create a temporary fasta file.
            contigs_file = join(tmp_dir, "MetaVIR_{0}.txt".format(bin_id))
            temp_file = join(tmp_dir, "MetaVIR_{0}.fa".format(bin_id))
            with open(contigs_file, "w") as f:
                for contig_name in list_contigs_name:
                    f.write("%s\n" % contig_name)
            cmd = "pyfastx extract {0} -l {1} > {2}".format(
                fasta, contigs_file, temp_file
            )
            process = sp.Popen(cmd, shell=True)
            process.communicate()

            # Concatenated the fasta in one sequence entry for checkV with 180
            # "N" spacer between contigs.
            with open(temp_file, "r") as tmp:
                start = True
                for line in tmp:
                    if line.startswith(">"):
                        if start:
                            start = False
                            if bin_id == 1:
                                out.write(
                                    ("%s\n" % ">MetaVIR_{0}".format(bin_id))
                                )
                            else:
                                out.write(
                                    ("\n%s\n" % ">MetaVIR_{0}".format(bin_id))
                                )
                        else:
                            out.write(
                                "N" * 60
                                + "\n"
                                + "N" * 60
                                + "\n"
                                + "N" * 60
                                + "\n"
                            )
                    else:
                        out.write(line)
    logger.info("{0} bins have been extracted".format(nb_bins))


def run_checkv(checkv_db, fasta, out_dir, remove_tmp, threads):
    """Function to launch end to end workflow from checkV.

    Parameters:
    -----------
    checkv_db : str
        Path to the directory of the reference database.
    fasta : str
        Path to the fasta of phages sequences to check.
    out_dir : str
        Path to the checkV output directory where the results of checkV will be
        written.
    remove_tmp : bool
        If True, remove temporary files from checkV.
    threads : int
        Number of threads to use to launch checkV.
    """

    # Defined checkV arguments.
    checkv_args = {
        "db": checkv_db,
        "input": fasta,
        "output": out_dir,
        "quiet": True,
        "remove_tmp": False,
        "restart": True,
        "threads": threads,
    }
    # Run checkV.
    checkv.modules.end_to_end.main(checkv_args)

    # Remove temporary directory if required. This is done separately as it
    # raises an error in checkV.
    if remove_tmp:
        shutil.rmtree(join(out_dir, "tmp"))


def run_metabat(
    contigs_file, input_fasta, outfile, phage_depth_file, temp_fasta
):
    """Function to launch metabat binning which is based on sequence and
    coverage information.

    Paremeters:
    -----------
    contigs_file : str
        Path to the file with the list of the phages contigs in the same order
        as the depth file.
    input_fasta : str
        Path to the fasta containing the phages sequences. It ocould have more
        sequences.
    outfile : str
        Path to write the clustering results of metabat.
    phage_depth_file : str
        Path to the depth information of the phages file.
    temp_fasta : str
        Path to write a temporary fasta with the phages sequences in the same
        order as the depth file.

    Returns:
    --------
    pandas.DataFrame:
        Table with the phage contigs name as index and clustering result column.
    """

    # Extract fasta to have sequences at the same order as the depth file.
    cmd = "pyfastx extract {0} -l {1} > {2}".format(
        input_fasta, contigs_file, temp_fasta
    )
    process = sp.Popen(cmd, shell=True)
    process.communicate()

    # Run metabat2 without the bin output with no limit of bin size and save
    # cluster information in the output file.
    cmd = "metabat2 -i {0} -a {1} -o {2} -s 0 --saveCls --noBinOut".format(
        temp_fasta, phage_depth_file, outfile
    )
    process = sp.Popen(cmd, shell=True)
    process.communicate()

    # Import metabat result as a pandas dataframe and return it.
    metabat = pd.read_csv(
        outfile, sep="\t", index_col=False, names=["Name", "Metabat_bin"]
    )
    return metabat


def phage_binning(
    checkv_db,
    depth_file,
    fasta_phages_contigs,
    phages_data_file,
    out_dir,
    remove_tmp,
    threads,
    tmp_dir,
):
    """Main function to bin phages contigs.

    Generates a fasta file where each entry is a phage bin, with 180bp "N" as
    spacers between contigs.

    Parameters:
    -----------
    checkv_db : str
        Path to the directory of the reference database.
    depth_file : str
        Path to depth file of the whole metagenome from metabat script.
    fasta_phages_contigs : str
        Path to the fasta containing the phages sequences. It could contain
        other sequences.
    phages_data_file : str
        Path to the output file from metavir host detection workflow.
    out_dir : str
        Path to the directory where to write the output data.
    remove_tmp : bool
        If eneabled, remove temporary files of checkV.
    threads : int
        Number of threads to use for checkV.
    tmp_dir : str
        Path to temporary directory for intermediate files.
    """

    # Create output and temporary files.
    phage_depth_file = join(tmp_dir, "phage_depth.txt")
    contigs_file = join(tmp_dir, "phage_contigs.txt")
    temp_fasta = join(tmp_dir, "phages.fa")
    metabat_output = join(tmp_dir, "metabat_phages_binning.tsv")
    phage_data_file = join(out_dir, "phages_data_final.tsv")
    fasta_phages_bins = join(out_dir, "phages_binned.fa")
    checkv_dir_contigs = join(out_dir, "checkV_contigs")
    checkv_dir_bins = join(out_dir, "checkV_bins")
    figure_file_pie = join(out_dir, "pie_phage_bins_size_distribution.png")
    figure_file_bar_size = join(
        out_dir, "barplot_phage_bins_size_distribution.png"
    )
    # figure_file_bar_nb = join(
    #     out_dir, "barplot_phage_bins_numbers_distribution.png"
    # )

    # Import host data from the previous step.
    phages_data = pd.read_csv(phages_data_file, sep="\t", index_col=0)

    # Launch metabat binning.
    build_phage_depth(contigs_file, depth_file, phages_data, phage_depth_file)
    metabat = run_metabat(
        contigs_file,
        fasta_phages_contigs,
        metabat_output,
        phage_depth_file,
        temp_fasta,
    )

    # Make binning based on both metabat binning and host detection.
    phages_data = phages_data.merge(metabat)
    phages_data, phage_bins = generate_phage_bins(phages_data)
    mio.write_phage_data(phages_data, phage_data_file)

    # Generate fasta for checkV quality check.
    generate_phages_fasta(
        fasta_phages_contigs, phage_bins, fasta_phages_bins, tmp_dir
    )

    # Run checkV on phage contigs and bins.
    run_checkv(checkv_db, temp_fasta, checkv_dir_contigs, remove_tmp, threads)
    run_checkv(
        checkv_db, fasta_phages_bins, checkv_dir_bins, remove_tmp, threads
    )

    # Plot figures
    checkv_summary_contigs = pd.read_csv(
        join(checkv_dir_contigs, "quality_summary.tsv"), sep="\t"
    )
    checkv_summary_bins = pd.read_csv(
        join(checkv_dir_bins, "quality_summary.tsv"), sep="\t"
    )
    mtf.pie_bins_size_distribution(checkv_summary_bins, figure_file_pie)
    mtf.barplot_bins_size(
        ["Contigs", "Bins"],
        [checkv_summary_contigs, checkv_summary_bins],
        figure_file_bar_size,
    )
    # mtf.barplot_bins_number(
    #     ["Contigs", "Bins"],
    #     [checkv_summary_contigs, checkv_summary_bins],
    #     figure_file_bar_nb,
    # )
