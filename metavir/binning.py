#!/usr/bin/env python3
# coding: utf-8

"""Module with the phages binning functions. 

It generates phages bins from the detected bacterial host of the contigs and 
from the metabat2 clusterization based on the sequence and shotgun coverage 
information. 

Core function to partition phages contigs:
    - build_matrix
    - build_phage_depth
    - generate_bin_summary
    - generate_phages_bins_metabat
    - generate_phages_bins_pairs
    - generate_phages_fasta
    - phage_binning
    - resolve_matrix
    - run_checkv
    - run_metabat
    - shuffle_phage_bins
    - update_phage_data
"""


import checkv
import metavir.figures as mtf
import metavir.io as mio
import numpy as np
import pandas as pd
import pypairix
import subprocess as sp
from metavir.log import logger
from os.path import join
import shutil
from typing import List, Tuple


def build_matrix(contigs, contigs_size, pairs_files):
    """Function to extract the pairs from a set of contigs from pairs files. Run
    faster if the files are indexed. The contacts are stored in a raw matrix of
    contacts between the list of given contigs. The contacts beloz 1kb are
    removed.

    Parameters:
    -----------
    contigs : List of str
        List of the phage contigs names uses in the alignment.
    contigs_size : List of int
        List of the size in bp of the contigs (same order as the contigs).
    pairs_files : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.

    Return:
    -------
    np.array:
        Raw matrix of contacts between the given contigs.
    """

    # Initiation
    npairs = 0
    n = len(contigs)
    mat = np.zeros((n, n))
    # Write one pair file for all the ones given.
    for pairs_file in pairs_files:
        # Check if the pairix index exist
        try:
            pairs_data = pypairix.open(pairs_file)
            pypairix_index = True
        except pypairix.PairixError:
            logger.warning("No pairix index found. Iterates on the pairs.")
            pypairix_index = False
        # Need a sorted (chr1 chr2 pos1 pos2) pair file indexed with pairix.
        if pypairix_index:
            for i, contig in enumerate(contigs):
                # Only need to retrieve the upper triangle.
                for j in range(i, len(contigs)):
                    pairs_lines = pairs_data.query2D(
                        contig,
                        0,
                        contigs_size[i],
                        contigs[j],
                        0,
                        contigs_size[j],
                        1,
                    )
                    for p in pairs_lines:
                        npairs += 1
                        # The threshold of 1000 is to remove the close range
                        # contacts.
                        if i == j:
                            if np.abs(int(p[2]) - int(p[4])) > 1000:
                                mat[i, i] += 1
                        else:
                            mat[i, j] += 1
        # else Iterates on the input pairs file (take much longer than with
        # the index).
        else:
            with open(pairs_file, "r") as input_pairs:
                for pairs_line in input_pairs:
                    # Ignore header lines.
                    if pairs_line.startswith("#"):
                        continue
                    # Split the line on the tabulation and check if both contigs
                    # are in the bin.
                    pairs = pairs_line.split("\t")
                    if pairs[1] in contigs and pairs[3] in contigs:
                        npairs += 1
                        i = contigs.index(pairs[1])
                        j = contigs.index(pairs[3])
                        # The threshold of 1000 is to remove the close range
                        # contacts.
                        if i == j:
                            if np.abs(int(pairs[2]) - int(pairs[4])) > 1000:
                                mat[i, i] += 1
                        else:
                            mat[i, j] += 1
    logger.info(f"{npairs} pairs extracted.")
    return mat


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


def generate_bin_summary(
    contig_data: "pd.DataFrame", phage_bins: dict, outfile: str
) -> "pd.DataFrame":
    """Function to generate and write the phage binning summary.

    Parameters
    ----------
    contig_data : pd.DataFrame
        Table with contig information from metator.
    phage_bins : dict
        Dictionnary with the phage bin id as key and the list of the contigs
        name as value.
    outfile : str
        Output file where to write the summary table. 

    Return
    ------
    pd.DataFrame :
        Summary table of the vrial bins. 
    """
    # Create output empty DataFrame
    cols = {
        'BinName': pd.Series(dtype='str'),
        'BinLength': pd.Series(dtype='int'),
        'GC': pd.Series(dtype='float'),
        'Score': pd.Series(dtype='float'),
        'ContigsNumber': pd.Series(dtype='int'),
        'MetagenomicBin': pd.Series(dtype='str'),
        'Hit': pd.Series(dtype='int'),
        'Contigs': pd.Series(dtype='str'),
    }
    summary = pd.DataFrame(cols, index=phage_bins.keys())

    # Iterates on the bins to fill the table.
    for bin_id in phage_bins.keys():
        summary.loc[bin_id, 'BinName'] = f'MetaVir_{bin_id:05d}'
        contigs = phage_bins[bin_id]['Contigs']
        summary.loc[bin_id, 'ContigsNumber'] = len(contigs)
        summary.loc[bin_id, 'Contigs'] = ','.join(contigs)
        summary.loc[bin_id, 'Score'] =  phage_bins[bin_id]['Score']
        summary.loc[bin_id, 'MetagenomicBin'] =  phage_bins[bin_id]['Bin']
        length, hit, gc = 0, 0, 0
        for contig in contigs:
            length += contig_data.loc[contig, 'Size']
            hit += contig_data.loc[contig, 'Hit']
            gc += contig_data.loc[contig, 'GC_content'] * contig_data.loc[contig, 'Size']
        summary.loc[bin_id, 'BinLength'] = length
        summary.loc[bin_id, 'Hit'] = hit
        summary.loc[bin_id, 'GC'] = gc / length
    
    # Write the summary.
    summary.to_csv(
        outfile, 
        sep='\t',
        na_rep='NA',
        float_format='%.2f',
        header=True,
        index=False,
    )
    return summary


def generate_phage_bins_metabat(phages_data):
    """Generates the binning of the phages contigs based on both HiC
    information (host detection) and the coverage and sequences information
    (metabat2 binning).

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
            phage_bins[bin_id_old]['Contig'].append(phages_data.loc[contig, "Name"])
        # Increment the bin id if it's the first time the phage id have been
        # seen.
        except KeyError:
            bin_id += 1
            bins_ids[phage_id] = bin_id
            phages_data.loc[contig, "MetaVir_bin"] = bin_id
            phage_bins[bin_id]['Contig'] = [phages_data.loc[contig, "Name"]]
            phage_bins[bin_id]['Score'] = np.nan 
    return phages_data, phage_bins


def generate_phage_bins_pairs(phages_data, pairs_files):
    """Generates the binning of the phages contigs based on both HiC
    information (host detection) and the coverage and sequences information
    (metabat2 binning).

    Parameters:
    -----------
    phages_data : pandas.DataFrame
        Table with the contigs name as index and with information from both
        host detected from metavir host and cluster form metabat2.
    pairs_file : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.

    Returns:
    --------
    pandas.DataFrame:
        Input table with the phage bin id column added.
    dict:
        Dictionnary with the phage bin id as key and the list of the contigs
        name as value.
    """

    # Extract contigs and contigs size.
    contigs = list(phages_data.Name)
    contigs_size = list(phages_data.Size)

    # Build matrix
    mat = build_matrix(contigs, contigs_size, pairs_files)

    # Associates contigs
    bins = resolve_matrix(mat)

    # Update phages_data and generates bins.
    phages_data, phage_bins = update_phage_data(phages_data, bins)

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
            list_contigs_name = phage_bins[bin_id]['Contigs']
            nb_bins += 1
            # Create a temporary fasta file.
            contigs_file = join(tmp_dir, f"MetaVIR_{bin_id:05d}.txt")
            temp_file = join(tmp_dir, f"MetaVIR_{bin_id:05d}.fa")
            with open(contigs_file, "w") as f:
                for contig_name in list_contigs_name:
                    f.write("%s\n" % contig_name)
            cmd = f"pyfastx extract {fasta} -l {contigs_file} > {temp_file}"
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
                            out.write((f">MetaVIR_{bin_id:05d}\n"))
                        else:
                            out.write(
                                "N" * 200
                            )
                    else:
                        out.write(line)
    logger.info(f"{nb_bins} bins have been extracted")


def phage_binning(
    checkv_db,
    depth_file,
    fasta_phages_contigs,
    out_dir,
    pairs_files,
    phages_data_file,
    plot,
    remove_tmp,
    threads,
    tmp_dir,
    method="pairs",
    random=False,
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
    method : str
        Method for the phage binning. Either 'pairs', 'metabat'.
    out_dir : str
        Path to the directory where to write the output data.
    pairs_files : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.
    phages_data_file : str
        Path to the output file from metavir host detection workflow.
    plot : bool
        If True make some summary plots.
    random : bool
        If enabled, make a andom shuffling of the bins.
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
    phage_data_file = join(out_dir, "phages_bin_summary.tsv")
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

    if method == "pairs":

        with open(contigs_file, "w") as f:
            for contig_name in list(phages_data.Name):
                f.write("%s\n" % contig_name)
        # Extract fasta to have sequences at the same order as the depth file.
        cmd = "pyfastx extract {0} -l {1} > {2}".format(
            fasta_phages_contigs, contigs_file, temp_fasta
        )
        process = sp.Popen(cmd, shell=True)
        process.communicate()
        phages_data, phage_bins = generate_phage_bins_pairs(
            phages_data, pairs_files
        )

    if method == "metabat":
        # Launch metabat binning.
        build_phage_depth(
            contigs_file, depth_file, phages_data, phage_depth_file
        )
        metabat = run_metabat(
            contigs_file,
            fasta_phages_contigs,
            metabat_output,
            phage_depth_file,
            temp_fasta,
        )

        # Make binning based on both metabat binning and host detection.
        phages_data = phages_data.merge(metabat)
        phages_data, phage_bins = generate_phage_bins_metabat(phages_data)

    if random:
        # Shuffle to simulate random bins. Uncomment to do it
        phages_data, phage_bins = shuffle_phage_bins(phages_data)

    # Generate fasta for checkV quality check.
    generate_phages_fasta(
        fasta_phages_contigs, phage_bins, fasta_phages_bins, tmp_dir
    )

    for bin_id in phage_bins:
        if association:
            phage_bins[bin_id]["Bin"] = mtb.asociate(phage_bins)
        else:
            phage_bins[bin_id]["Bin"] = np.nan

    summary = generate_bin_summary(phages_data, phage_bins, phage_data_file)

    # Run checkV on phage contigs and bins.
    if plot:
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


def resolve_matrix(mat: "np.ndarray", threshold: float = 1.0) -> List(Tuple):
    """Main function to bin phages contigs.

    From the marix of contacts associates the contigs with a lot of
    interactions. To normalize the contacts we divide the contacts by the
    geometric mean of the contacts in intra. To avoid the noise of the close
    range contacts we have remove the contacts below 1000bp. When two contigs
    are binned, they are fused together.

    Parameters:
    -----------
    mat : np.array
        Matrix of the raw contacts between the contigs. Upper triangle and the
        contacts in intra below 1000bp are not kept.
    threshold : float
        Threshold of score to bin contigs. [Default: 1]

    Returns:
    List of tuple:
        List of the Tuple of the associated bins.
    """

    bins = []
    n = len(mat)
    # Save a copy of the initial matrix to keep the intra initial values.
    mat0 = np.copy(mat)

    # Normalize values
    for i in range(n):
        for j in range(i + 1, n):
            if (mat[i, i] > 0) and (mat[j, j] > 0):
                mat[i, j] = mat[i, j] / np.sqrt(mat[i, i] * mat[j, j])
            else:
                mat[i, j] = 0

    # Remove the intra contacts.
    for i in range(n):
        mat[i, i] = 0

    # While there is an association bigger than 1 associates the contigs.
    maxi = np.max(mat)
    while maxi > threshold:
        # Handle the case where we have multiple points at the maximum value.
        try:
            i, j = map(int, np.where(mat == maxi))
        except TypeError:
            i, j = map(
                int, [np.where(mat == maxi)[0][0], np.where(mat == maxi)[1][0]]
            )

        # Compute the sum of the count to create the merge vector.
        a = mat0[i, :] + mat0[:, i] + mat0[j, :] + mat0[:, j]
        # Compute a new intra count.
        intra = mat0[i, i] + mat0[j, j] + mat0[i, j]
        mat0[i, i] = intra
        # Normalize the new vector.
        for k in range(n):
            if k < i:
                mat0[k, i] = a[k]
                if mat0[k, k] > 0:
                    mat[k, i] = a[k] / np.sqrt(intra * mat0[k, k])
                else:
                    mat[k, i] = 0
            elif k > i:
                mat0[i, k] = a[k]
                if mat0[k, k] > 0:
                    mat[i, k] = a[k] / np.sqrt(intra * mat0[k, k])
                else:
                    mat[k, i] = 0
        # Remove the contig j as it has been merged with the contig i.
        mat[j, :] = np.zeros(n)
        mat[:, j] = np.zeros(n)
        mat0[j, :] = np.zeros(n)
        mat0[:, j] = np.zeros(n)

        bins.append([i, j, maxi])
        maxi = np.max(mat)

    return bins


def shuffle_phage_bins(phages_data):
    """Function to shuffle id to imitate a random binning with the same bins
    distribution as the one created by MetaVir.

    Parameters:
    -----------
    phages_data : pandas.DataFrame
        Table with the contigs name as index and with information from both
        host detected from metavir host and cluster form metabat2 and with phage
        bins id.

    Returns:
    --------
    pandas.DataFrame:
        Input table with the phage bin id column randomly shuffled.
    dict:
        Dictionnary with the phage bin id as key and the list of the contigs
        name as value from the shuffle contigs bins.
    """

    # Shuffle the ids of the dataframe
    phages_bin_ids = phages_data.MetaVir_bin
    shuffle_ids = np.random.permutation(phages_bin_ids)
    phages_data["shuffle"] = shuffle_ids

    phages_bins = {}
    # Shuffle phages bins according to the shuffle dataframe.
    for index in phages_data.index:
        contig = phages_data.loc[index, "Name"]
        phage_bin_id = phages_data.loc[index, "shuffle"]
        try:
            phages_bins[phage_bin_id]['Contig'] .append(contig)
        except KeyError:
            phages_bins[phage_bin_id]['Contig']  = [contig]
            phages_bins[phage_bin_id]['Score'] = np.nan
    return phages_data, phages_bins


def update_phage_data(phages_data, bins):
    """Function to update the paheg bins data.

    Parameters
    ----------
    phages_data : pd.DataFrame
        Table wit the phage contig information.
    bins : list of tuple
        List of pairs of contigs with their score of association.

    Returns
    -------
    pandas.DataFrame :
        Table wit the phage contig information updated.
    dictionnnary :
        Dictionnary of the phage bins.
    """
    # Initiation
    phages_data["MetaVir_bin"] = 0
    phages_data["MetaVir_Score"] = 0
    bin_id = 0
    phage_bins = {}
    phages_data.set_index(np.arange(len(phages_data)), inplace=True)
    for contig_tuple in bins:
        i, j, score = contig_tuple
        # If no existing bin, creates one.
        if (phages_data.loc[i, "MetaVir_bin"] == 0) and (
            phages_data.loc[j, "MetaVir_bin"] == 0
        ):
            bin_id += 1
            current_bin = bin_id
            current_score = score
        # If one existing bin, append it.
        elif (phages_data.loc[i, "MetaVir_bin"] == 0) or (
            phages_data.loc[j, "MetaVir_bin"] == 0
        ):
            current_bin = max(
                phages_data.loc[i, "MetaVir_bin"],
                phages_data.loc[j, "MetaVir_bin"],
            )
            current_score = min(score, max(
                phages_data.loc[i, "MetaVir_Score"],
                phages_data.loc[j, "MetaVir_Score"],
            ))
        # Complicated case as we have to fuse two existing bins. The bin j
        # is not reused.
        else:
            current_bin = phages_data.loc[i, "MetaVir_bin"]
            bin_j = phages_data.loc[j, "MetaVir_bin"]
            if current_bin > bin_j:
                current_bin -= 1
            bin_id -= 1
            for k in range(len(phages_data)):
                if phages_data.loc[k, "MetaVir_bin"] == bin_j:
                    phages_data.loc[k, "MetaVir_bin"] = current_bin
                elif phages_data.loc[k, "MetaVir_bin"] > bin_j:
                    phages_data.loc[k, "MetaVir_bin"] -= 1
            current_score = min(
                score,
                phages_data.loc[i, "MetaVir_Score"],
                phages_data.loc[j, "MetaVir_Score"],
            )   
        phages_data.loc[i, "MetaVir_bin"] = current_bin
        phages_data.loc[j, "MetaVir_bin"] = current_bin
        phages_data.loc[i, "MetaVir_Score"] = current_score
        phages_data.loc[j, "MetaVir_Score"] = current_score

    # Update unbinned contigs and build phage bins.
    for k in phages_data.index:
        if phages_data.loc[k, "MetaVir_bin"] == 0:
            bin_id += 1
            phages_data.loc[k, "MetaVir_bin"] = bin_id
            phage_bins[bin_id] = {
                'Contigs': [phages_data.loc[k, "Name"]]
                'Score': np.nan
            }
        else:
            curr_bin = phages_data.loc[k, "MetaVir_bin"]
            try:
                phage_bins[curr_bin]['Contig'].append(
                    phages_data.loc[k, "Name"]
                )
                phage_bins[curr_bin]['Score'] = min(
                    phage_bins[curr_bin]['Score'], 
                    phages_data.loc[k, "MetaVir_Score"]
                )
            except KeyError:
                phage_bins[curr_bin] = {
                'Contigs': [phages_data.loc[k, "Name"]]
                'Score': phages_data.loc[k, "MetaVir_Score"]
            }
    return phages_data, phage_bins
