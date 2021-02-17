# MpGAP (General multi-platform genome assembly pipeline)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3997375.svg)](https://doi.org/10.5281/zenodo.3445485) ![](https://img.shields.io/github/v/release/fmalmeida/MpGAP) [![Build Status](https://travis-ci.org/fmalmeida/MpGAP.svg?branch=master)](https://travis-ci.org/fmalmeida/MpGAP) ![](https://img.shields.io/badge/dependencies-docker-informational) [![Documentation Status](https://readthedocs.org/projects/mpgap/badge/?version=latest)](https://mpgap.readthedocs.io/en/latest/?badge=latest) ![](https://img.shields.io/badge/Nextflow-v20.01-yellowgreen)

MpGAP is an easy to use nextflow docker-based pipeline that adopts well known software for genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes. This pipeline wraps up the following software:

* [Canu](https://github.com/marbl/canu)
* [Flye](https://github.com/fenderglass/Flye)
* [Unicycler](https://github.com/rrwick/Unicycler)
* [Spades](https://github.com/ablab/spades)
* [Nanopolish](https://github.com/jts/nanopolish)
* [Medaka](https://github.com/nanoporetech/medaka)
* [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
* [Pilon](https://github.com/broadinstitute/pilon)
* [QUAST](https://github.com/ablab/quast)

## Further reading

This pipeline has two complementary pipelines (also written in nextflow) for [NGS preprocessing](https://github.com/fmalmeida/mpgap) and [annotation](https://github.com/fmalmeida/bacannot) that can give the user a complete workflow for bacterial genomics analyses.

## Table of contents

* [Requirements](https://github.com/fmalmeida/mpgap#requirements)
* [Installation](https://github.com/fmalmeida/mpgap#installation)
* [Documentation](https://github.com/fmalmeida/mpgap#documentation)
  * [Hybrid strategies explanation](https://github.com/fmalmeida/mpgap#explanation-of-hybrid-strategies)
  * [Full usage](https://github.com/fmalmeida/mpgap#usage)
  * [Command line examples](https://github.com/fmalmeida/mpgap#command-line-usage-examples)
  * [Configuration File](https://github.com/fmalmeida/mpgap#using-the-configuration-file)
  * [Interactive and graphical execution](https://github.com/fmalmeida/mpgap#interactive-graphical-configuration-and-execution)
* [Known issues](https://github.com/fmalmeida/mpgap#known-issues)

## Requirements

This pipeline has only two dependencies: [Docker](https://www.docker.com) and [Nextflow](https://github.com/nextflow-io/nextflow).

* Unix-like operating system (Linux, macOS, etc)
* Java 8 (or higher)
* Nextflow (version 20.01 or higher)
* Docker
  * Image: `fmalmeida/mpgap`

## Installation

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).
    * You can give this [in-house script](https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh) a try.
    * After installed, you need to download the required Docker image

          docker pull fmalmeida/mpgap

> Each release is accompanied by a Dockerfile in the docker folder. When using releases older releases, users can create the correct image using
the Dockerfile that goes alongside with the release (Remember to give the image the correct name, as it is in dockerhub and the nextflow script).
The latest release will always have its docker image in dockerhub.

2. Install Nextflow (version 20.01 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow run fmalmeida/mpgap --help

> Users can get let the pipeline always updated with: `nextflow pull fmalmeida/MpGAP`

## Documentation

### Explanation of hybrid strategies

Hybrid assemblies can be produced using one of two available strategies:

#### Strategy 1

By using Unicycler and/or SPAdes hybrid assembly modes. For instance, it will use Unicycler hybrid mode which will first assemble a high quality assembly graph with Illumina data and then it will use long reads to bridge the gaps. More information about Unicycler Hybrid mode can be found [here](https://github.com/rrwick/Unicycler#method-hybrid-assembly).

#### Strategy 2

By polishing a long reads only assembly with Illumina reads. For that, users will have to set `--illumina_polish_longreads_contigs` to true. This will tell the pipeline to produce a long reads only assembly (with canu, flye or unicycler) and polish it with Pilon (for unpaired reads) or with [Unicycler-polish program](https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md) (for paired end reads).

> Note that, `--illumina_polish_longreads_contigs` parameter is an alternative workflow, when used, it will execute ONLY strategy 2 and not both strategies. When false, only strategy 1 will be executed.

> Remember to select the desired assemblers to run with `--try_canu`, `--try_flye` and/or `--try_unicycler`

#### Example:

        nextflow run fmalmeida/mpgap --outdir output --threads 5 --assembly_type hybrid --try_unicycler --try_flye --try_canu
        --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq" --shortreads_single "dataset_1/sampled/illumina_single.fastq"
        --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq" --illumina_polish_longreads_contigs

### Usage

* Complete command line explanation of parameters:
    + `nextflow run fmalmeida/mpgap --help`
* See usage examples in the command line:
    + `nextflow run fmalmeida/mpgap --examples`
* However, users are encouraged to read the [complete online documentation](https://mpgap.readthedocs.io/en/latest/?badge=latest).

### Command line usage examples

Command line executions are exemplified [in the manual](https://mpgap.readthedocs.io/en/latest/examples.html).

### Using the configuration file

All parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is used the pipeline is executed as `nextflow run fmalmeida/mpgap -c ./configuration-file`. Your configuration file is what will tell the pipeline which type of data you have, and which processes to execute. Therefore, it needs to be correctly configured.

To create a configuration file in your working directory:

* For Hybrid assemblies:

      nextflow run fmalmeida/mpgap --get_hybrid_config

* For Long reads only assemblies:

      nextflow run fmalmeida/mpgap --get_lreads_config

* For illumina only assemblies:

      nextflow run fmalmeida/mpgap --get_sreads_config

### Interactive graphical configuration and execution

Users can trigger a graphical and interactive pipeline configuration and execution by using [nf-core launch](https://nf-co.re/launch) utility.

#### Install nf-core

```bash
# Install nf-core
pip install nf-core>=1.10
```

#### launch the pipeline

nf-core launch will start an interactive form in your web browser or command line so you can configure the pipeline step by step and start the execution of the pipeline in the end.

```bash
# Launch the pipeline
nf-core launch fmalmeida/mpgap
```

It will result in the following:

<p align="center">
<img src="./images/nf-core-asking.png" width="500px"/>
</p>

<p align="center">
<img src="./images/nf-core-gui.png" width="400px"/>
</p>

#### nextflow tower

This pipeline also accepts that users track its execution of processes via [nextflow tower](https://tower.nf/). For that users will have to use the parameters `--use_tower` and `--tower_token`.

## Known issues

1. Whenever using unicycler with unpaired reads, an odd platform-specific SPAdes-related crash seems do randomly happen as it can be seen in the issue discussed at https://github.com/rrwick/Unicycler/issues/188. As a workaround, Ryan says to use the `--no_correct` parameter which solves the issue and does not have a negative impact on assembly quality. Therefore, whenever using Illumina unpaired reads, this parameter is automatically used by the pipeline.
2. Users are advised to do not perform hybrid annotations with multiple samples (giving reads as glob patterns, such as ? and *). This would result in unexpected outcomes since the pipeline is not trained to understand/search pairs when multiple samples are given. It is possible that the reads get mixed, or that reads are analysed separately.
   + These glob patterns are effective when running the pipeline with multiple samples that have the same type of input reads. For example, running a long reads assembly for multiple samples with nanopore data (using `--longreads *.fastq --lr_type nanopore`). This is possible because, as said before, the pipeline is not trained to understand/search pairs when multiple samples are given, which means, that each fastq found with the glob pattern will be individually assembled using the parameters set.

## Citation

To cite this pipeline users can use our Zenodo tag or directly via the github url.

Users are encouraged to cite the programs used in this pipeline whenever they are used:

* [Canu](https://github.com/marbl/canu)
* [Flye](https://github.com/fenderglass/Flye)
* [Unicycler](https://github.com/rrwick/Unicycler)
* [Spades](https://github.com/ablab/spades)
* [Nanopolish](https://github.com/jts/nanopolish)
* [QUAST](https://github.com/ablab/quast)
* [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
* [Pilon](https://github.com/broadinstitute/pilon).
