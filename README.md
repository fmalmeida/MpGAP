# MpGAP (General multi-platform genome assembly pipeline)

![](https://travis-ci.com/fmalmeida/MpGAP.svg?branch=master) [![DOI](https://zenodo.org/badge/200904121.svg)](https://zenodo.org/badge/latestdoi/200904121) [![Documentation Status](https://readthedocs.org/projects/mpgap/badge/?version=latest)](https://mpgap.readthedocs.io/en/latest/?badge=latest) ![](https://img.shields.io/docker/cloud/build/fmalmeida/mpgap) ![](https://img.shields.io/badge/Nextflow-v20.01-yellowgreen)



MpGAP is a nextflow docker-based pipeline that wraps up [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler), [Spades](https://github.com/ablab/spades), [Nanopolish](https://github.com/jts/nanopolish), [Medaka](https://github.com/nanoporetech/medaka), [QUAST](https://github.com/ablab/quast), [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus) and [Pilon](https://github.com/broadinstitute/pilon).

This is an easy to use pipeline that adopts well known software for genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes.

This pipeline has only two dependencies: [Docker](https://www.docker.com) and [Nextflow](https://github.com/nextflow-io/nextflow).

Also, check out my other two complementary pipelines for [preprocessing NGS raw data](https://github.com/fmalmeida/NGS-preprocess) and to annotate microbial genomes.

## Table of contents

* [Requirements](https://github.com/fmalmeida/MpGAP#requirements)
* [Quickstart](https://github.com/fmalmeida/MpGAP#quickstart)
* [Documentation](https://github.com/fmalmeida/MpGAP#documentation)
  * [Workflow explanation](https://github.com/fmalmeida/MpGAP#workflow)
  * [Full usage](https://github.com/fmalmeida/MpGAP#usage)
  * [Configuration File](https://github.com/fmalmeida/MpGAP#using-the-configuration-file)

## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Java 8
* Docker
  * `fmalmeida/mpgap`

## Quickstart

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).
    * You can give this [in-house script](https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh) a try.
    * After installed, you need to download the required Docker image

          docker pull fmalmeida/mpgap

2. Install Nextflow (version 0.24.x or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow fmalmeida/MpGAP --help

## Documentation

### Workflow

This pipeline was designed to accept data from Illumina, Pacbio and Oxford Nanopore Technologies sequence platforms. It is possible to perform short or long reads only assemblies as well as hybrid assemblies.

#### Short reads only

Short reads only assemblies can be produced with Unicycler and/or SPAdes assemblers. Since short reads have an extremely high basecall accuracy, no further steps are done.

Users can assemble Illumina paired or single end reads. [Here](https://github.com/fmalmeida/MpGAP#usage-examples) you can find a few usage examples.

#### Long reads only

Long reads only assemblies can be produced with Canu, Flye and/or Unicycler assemblers. Users can assemble pacbio or nanopore reads.

Addionally, users can perform a polishing step since long reads only assemblies generally suffer from basecall accuracy. For that, it is only necessary to set path to a directory containing raw nanopore FAST5 data to polish with Nanopolish or, for pacbio, to set path to \*.bax.h5 or \*.subreads.bam files to use VarianCaller to polish the assembly.

[Here](https://github.com/fmalmeida/MpGAP#usage-examples) you can find a few usage examples.

#### Hybrid

It is possible to perform two types of hybrid assemblies.

1. By using Unicycler and SPAdes hybrid assembly modes. For instance, will use Unicycler hybrid mode which will first assemble a high quality assembly graph with Illumina data and then it will use long reads to bridge the gaps. More information about Unicycler Hybrid mode can be found [here](https://github.com/rrwick/Unicycler#method-hybrid-assembly).
2. By polishing a long reads only assembly with Illumina reads. For that, users will have to set **illumina_polish_longreads_contigs** to true. This will tell the pipeline to produce a long reads only assembly and polish it with Pilon (for unpaired reads) or with [Unicycler-polish program](https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md) (for paired end reads).

Note that, **illumina_polish_longreads_contigs** parameter is an addition, when true, it will execute both strategies 1 and 2. Then false, only strategy 1 will be executed.

[Here](https://mpgap.readthedocs.io/en/latest/examples.html#examples) you can find a few usage examples.

## Usage

Users are encouraged to read the [documentation](https://mpgap.readthedocs.io/en/latest/index.html).

### Using the configuration file

All parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is used the pipeline is executed as `nextflow run fmalmeida/MpGAP -c ./configuration-file`

Your configuration file is what will tell the pipeline which type of data you have, and which processes to execute. Therefore, it needs to be correctly configured.

Create a configuration file in your working directory:

* For Hybrid assemblies:

      nextflow run fmalmeida/MpGAP --get_hybrid_config

* For Long reads only assemblies:

      nextflow run fmalmeida/MpGAP --get_lreads_config

* For illumina only assemblies:

      nextflow run fmalmeida/MpGAP --get_sreads_config

* To download the YAML file used to pass additional parameters to assemblers:

      nextflow run fmalmeida/MpGAP --get_yaml

## Citation

Users are encouraged to cite the programs used in this pipeline whenever they are used. They are: [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler), [Spades](https://github.com/ablab/spades), [Nanopolish](https://github.com/jts/nanopolish), [QUAST](https://github.com/ablab/quast), [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus) and [Pilon](https://github.com/broadinstitute/pilon).
