[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3997375.svg)](https://doi.org/10.5281/zenodo.3445485) [![Releases](https://img.shields.io/github/v/release/fmalmeida/mpgap)](https://github.com/fmalmeida/mpgap/releases) [![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://mpgap.readthedocs.io/en/latest/?badge=latest) [![Dockerhub](https://img.shields.io/badge/Docker-fmalmeida/mpgap-informational)](https://hub.docker.com/r/fmalmeida/mpgap) [![Docker build](https://img.shields.io/docker/cloud/build/fmalmeida/mpgap)](https://hub.docker.com/r/fmalmeida/mpgap) ![Docker Pulls](https://img.shields.io/docker/pulls/fmalmeida/mpgap) [![Nextflow version](https://img.shields.io/badge/Nextflow%20>=-v20.07-important)](https://www.nextflow.io/docs/latest/getstarted.html) [![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/mpgap/blob/master/LICENSE)

<p align="center">

  <h1 align="center">MpGAP pipeline</h2>

  <p align="center">
    <h3 align="center">A generic multi-platform genome assembly pipeline</h3>
    <br />
    <a href="https://mpgap.readthedocs.io/en/latest/index.html"><strong>See the documentation »</strong></a>
    <br />
    <br />
    <a href="https://github.com/fmalmeida/MpGAP/issues">Report Bug</a>
    ·
    <a href="https://github.com/fmalmeida/MpGAP/issues">Request Feature</a>
  </p>
</p>

## About

MpGAP is an easy to use nextflow docker-based pipeline that adopts well known software for genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes. This pipeline wraps up the following software:

* [Canu](https://github.com/marbl/canu)
* [Flye](https://github.com/fenderglass/Flye)
* [Raven](https://github.com/lbcb-sci/raven)
* [Haslr](https://github.com/vpc-ccg/haslr)
* [Unicycler](https://github.com/rrwick/Unicycler)
* [Spades](https://github.com/ablab/spades)
* [Shovill](https://github.com/tseemann/shovill)
* [Nanopolish](https://github.com/jts/nanopolish)
* [Medaka](https://github.com/nanoporetech/medaka)
* [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
* [Pilon](https://github.com/broadinstitute/pilon)
* [QUAST](https://github.com/ablab/quast)
* [MultiQC](https://multiqc.info/)

## Further reading

This pipeline has two complementary pipelines (also written in nextflow) for [NGS preprocessing](https://github.com/fmalmeida/mpgap) and [annotation](https://github.com/fmalmeida/bacannot) that can give the user a complete workflow for bacterial genomics analyses.

## Requirements

This pipeline has only two dependencies: [Docker](https://www.docker.com) and [Nextflow](https://github.com/nextflow-io/nextflow).

* Unix-like operating system (Linux, macOS, etc)
  + Windows users maybe can execute it using the linux subsystem for windows as shown in:
    + https://docs.microsoft.com/pt-br/windows/wsl/install-win10
    + https://www.nextflow.io/docs/latest/getstarted.html
    + https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux
* Java 8 (or higher)
* Nextflow (version 20.01 or higher)
* Docker
  * Image: `fmalmeida/mpgap`

## Installation

1. If you don't have it already install [Docker](https://docs.docker.com/) in your computer.
    * After installed, you need to download the required Docker images

          docker pull fmalmeida/mpgap

> Each release is accompanied by a Dockerfile in the docker folder. When using releases older releases, users can create the correct image using
the Dockerfile that goes alongside with the release (Remember to give the image the correct name, as it is in dockerhub and the nextflow script).
The latest release will always have its docker image in dockerhub.

2. Install Nextflow (version 20.01 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow run fmalmeida/mpgap --help

> Users can let the pipeline always updated with: `nextflow pull fmalmeida/mpgap`

## Documentation

### Explanation of hybrid strategies

Hybrid assemblies can be produced using one of two available strategies:

#### Strategy 1

By using Unicycler, Haslr and/or SPAdes hybrid assembly modes. For instance, it can use the Unicycler hybrid mode which will first assemble a high quality assembly graph with Illumina data and then it will use long reads to bridge the gaps. More information about Unicycler Hybrid mode can be found [here](https://github.com/rrwick/Unicycler#method-hybrid-assembly).

#### Strategy 2

By polishing a long reads only assembly with Illumina reads. For that, users will have to set `--strategy_2` to true. This will tell the pipeline to produce a long reads only assembly (with canu, flye, raven or unicycler) and polish it with Pilon (for unpaired reads) or with [Unicycler-polish program](https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md) (for paired end reads).

> Note that, `--strategy_2` parameter is an alternative workflow, when used, it will execute ONLY strategy 2 and not both strategies. When false, only strategy 1 will be executed.

#### Example:

        nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
        --shortreads_single "path-to/illumina_unpaired.fastq" --lr_type 'nanopore' --longreads "path-to/ont_reads.fastq" --strategy_2

### Usage

<a href="https://mpgap.readthedocs.io/en/latest/index.html"><strong>Users are advised to read the complete documentation »</strong></a>

* Complete command line explanation of parameters:
    + `nextflow run fmalmeida/mpgap --help`
* See usage examples in the command line:
    + `nextflow run fmalmeida/mpgap --examples`

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

## Known issues

1. Whenever using unicycler with unpaired reads, an odd platform-specific SPAdes-related crash seems do randomly happen as it can be seen in the issue discussed at https://github.com/rrwick/Unicycler/issues/188.
  + As a workaround, Ryan says to use the `--no_correct` parameter which solves the issue and does not have a negative impact on assembly quality. Therefore, whenever using Illumina unpaired reads, this parameter is automatically used by the pipeline.
2. Whenever running the pipeline for multiple samples at once using glob patterns such as '*' and '?', users are advised to do not perform hybrid assemblies, nor combining both paired and unpaired short reads in short reads only assemblies. Because the pipeline is not yet trained to properly search for the correct pairs, and since nextflow channels are random, we cannot ensure that the combination of data used in these to assembly types will be right. The pipeline treats each input file as a unique sample, and it will execute it individually.
  + To date, the use of glob patterns only works properly with long reads only assembly, or short reads only assemblies using either paired or unpaired reads, not both at the same time. For example:
    + `nextflow run [...] --longreads 'my_data/*.fastq' --lr_type 'nanopore' --outdir my_results`
    + The pipeline will load and assembly each fastq in the `my_data` folder and assemble it, writing the results for each read in a sub-folder with the reads basename in the `my_results` output folder.
    + `nextflow run [...] --shortreads_single 'my_data/*.fastq' --outdir my_results`
    + The pipeline will load and assembly each fastq in the `my_data` folder and assemble it, writing the results for each read in a sub-folder with the reads basename in the `my_results` output folder.

> However, we are currently working in a proper way to execute the hybrid and combination of short reads in assemblies for multiple samples at once so that users can properly execute it without confusion. But it will come in v2.3.

3. Sometimes, shovill assembler can fail and cause the pipeline to fail due to problems in estimating the genome size. This, is actually super simple to solve! Instead of letting the shovill assembler estimate the genome size, you can pass the information to it and prevent its fail:
    + `--shovill_additional_parameters '--gsize 3m'`

## Citation

To cite this pipeline users can use our Zenodo tag or directly via the github url.

Users are encouraged to cite the programs used in this pipeline whenever they are used:

* [Canu](https://github.com/marbl/canu)
* [Flye](https://github.com/fenderglass/Flye)
* [Raven](https://github.com/lbcb-sci/raven)
* [Haslr](https://github.com/vpc-ccg/haslr)
* [Unicycler](https://github.com/rrwick/Unicycler)
* [Spades](https://github.com/ablab/spades)
* [Shovill](https://github.com/tseemann/shovill)
* [Nanopolish](https://github.com/jts/nanopolish)
* [Medaka](https://github.com/nanoporetech/medaka)
* [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
* [Pilon](https://github.com/broadinstitute/pilon)
* [QUAST](https://github.com/ablab/quast)
* [MultiQC](https://multiqc.info/)
