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
git clone https://github.com/ABignaud/MetaVir.git
pip3 install -e ./MetaVir
```

Installation of metabat2 and checkV database are necessary to use the binning module. To install metabat2, follow the instructions [here](https://bitbucket.org/berkeleylab/metabat/src/master/) and to install and/or update checkv database, instructions are available [here](https://bitbucket.org/berkeleylab/checkv/src/master/). To use MetaVir, it is mandatory to set the environnement variable CHECKVDB manually or using their workflows.

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

* `host` : Detect bacterial host from a metaHiC network binned by metaTOR given a annotated phages list.
* `binning` : Build phages MAGs based on metagenomic binning using [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) and the host detection form the metaHiC data.

There are a number of other, optional, miscellaneous actions:

* `pipeline` : Run both steps sequentially.
* `version` : display current version number.
* `help` : display help message.

### Output files

#### host

* **host.tsv:** Table with two columns the phage contigs name and the bacterial MAG host.

#### binning

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
