# MetaVir

Metagenomics phage bacterial host detection and phage binning pipeline based on metaHiC data and short reads metagenomic assembly.

The goal of this pipeline is to detect phage bacterial host and to build phage MAGs based on the contig host detected, and a classical binning of the contigs based on the covergae and the sequences.

## Installation

### Requirements

* Python 3.6 or later is required.
* The following librairies are required but will be automatically installed with the pip installation: `checkv`, `docopt`, `networkx`, `numpy`, `pandas`, `pyfastx`.
* The following software and database should be installed separately if you used pip installation:
  * [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)
  * [checkV](https://bitbucket.org/berkeleylab/checkv/src/master/) database

### Using pip

```sh
pip3 install metavir
```

Or to use the latest version

```sh
pip3 install -e git+https://github.com/ABignaud/MetaVir.git@master#egg=metator
```

### Using docker container

```sh
git clone https://github.com/ABignaud/MetaVir.git
cd MetaVir
docker build --tag metavir .
docker run metavir {hsot|binning} [parameters]
```

A docker image will be soon available.

## Usage

```sh
metavir {host|binning} [parameters]
```

There are two main steps in the metaVir pipeline, which must be run in the following order:

* `host` : Detect bacterial host in the network of a given annotated phages contigs list.
* `binning` : Build phages MAGs based on metagenomic binning using [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) and the host detection form the metaHiC data.

There are a number of other, optional, miscellaneous actions:

* `version` : display current version number.
* `help` : display help message.

### Output files

#### host

* **host.tsv:** Table with two columns the phage contigs name and the bacterial MAG host.

#### binning

* **metabat_phages_binning.tsv:** Clustering output from [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)
* **phages_data.tsv:** Clustering output summary from the binning.
* **phages_binned.fa:** Fasta file with the sequences of the binnes phage MAGs. Each entry represent one phage MAG, and contigs are delimited by 180bp "N" spacers.
* **checkV_contigs:** Directory with [checkV](https://bitbucket.org/berkeleylab/checkv/src/master/) output of phage contigs.
* **checkV_bins:** Directory with [checkV](https://bitbucket.org/berkeleylab/checkv/src/master/) output of phage bins.
* **Some plots:**
  * barplot_phage_bins_size_distribution.png
  * pie_phage_bins_size_distribution.png

## References

[MetaHiC phage-bacteria infection network reveals active cycling phages of the healthy human gut](https://elifesciences.org/articles/60608),  Martial Marbouty, Agnès Thierry, Gaël A Millot, Romain Koszul, 2021

## Contact

### Authors

* amaury.bignaud@pasteur.fr
* martial.marbouty@pasteur.fr
* romain.koszul@pasteur.fr

### Research lab

[Spatial Regulation of Genomes](https://research.pasteur.fr/en/team/spatial-regulation-of-genomes/) (Institut Pasteur, Paris)
