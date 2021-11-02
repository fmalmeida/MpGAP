<img src="images/lOGO_3.png" width="300px">

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

MpGAP is an easy to use nextflow docker-based pipeline that adopts well known software for _de novo_ genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes. This pipeline wraps up the following software:

|| **Source** |
|:- | :- |
| **Assemblers** | [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Raven](https://github.com/lbcb-sci/raven), [Shasta](https://github.com/chanzuckerberg/shasta), [wtdbg2](https://github.com/ruanjue/wtdbg2), [Haslr](https://github.com/vpc-ccg/haslr), [Unicycler](https://github.com/rrwick/Unicycler), [Spades](https://github.com/ablab/spades), [Shovill](https://github.com/tseemann/shovill) |
| **Polishers** | [Nanopolish](https://github.com/jts/nanopolish), [Medaka](https://github.com/nanoporetech/medaka), [gcpp](https://github.com/PacificBiosciences/gcpp), [Pilon](https://github.com/broadinstitute/pilon) |
| **Quality check** | [QUAST](https://github.com/ablab/quast), [MultiQC](https://multiqc.info/) |

### Release notes

Are you curious about changes between releases? See the [changelog](markdown/CHANGELOG.md).

* I **strongly**, **vividly**, **mightily** recommend the usage of the latest versions hosted in master branch, which is nextflow's default.
    + The latest will always have support, bug fixes and generally maitain the same processes (I mainly add things instead of removing) that also were in previous versions.
    + But, if you **really** want to execute an earlier release, please [see the instructions for that](markdown/earlier_releases_instructions.md).
* Versions below 2.0 are no longer supported.

### Further reading

This pipeline has two complementary pipelines (also written in nextflow) for [NGS preprocessing](https://github.com/fmalmeida/ngs-preprocess) and [prokaryotic genome annotation](https://github.com/fmalmeida/bacannot) that can give the user a complete workflow for bacterial genomics analyses.

## Requirements

This pipeline has only two dependencies: [Docker](https://www.docker.com) and [Nextflow](https://github.com/nextflow-io/nextflow).

* Unix-like operating system (Linux, macOS, etc)
  + Windows users maybe can execute it using the linux subsystem for windows as shown in:
    + https://nextflow.io/blog/2021/setup-nextflow-on-windows.html
    + https://docs.microsoft.com/pt-br/windows/wsl/install-win10
    + https://www.nextflow.io/docs/latest/getstarted.html
    + https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux
* Java 8 (or higher)
* Nextflow (version 20.01 or higher)
* Docker
  * Image: `fmalmeida/mpgap:v3.0`

## Installation

1. If you don't have it already install [Docker](https://docs.docker.com/) in your computer.
    * After installed, you need to download the required Docker images

          docker pull fmalmeida/mpgap:v3.0

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

        nextflow run fmalmeida/mpgap --output output --threads 5 --input "samplesheet.yml" --strategy_2

#### Tip

Users can perform both strategy 1 and strategy 2 hybrid assemblies, at the same time, for a sample by configuring the sample, in the samplesheet, with the key: `hybrid_strategy: both`. For more information please read: https://mpgap.readthedocs.io/en/latest/examples.html

### Usage

<a href="https://mpgap.readthedocs.io/en/latest/index.html"><strong>For understading pipeline usage and configuration, users must read the complete online documentation »</strong></a>

#### Warnings

* Remember to **always** write input paths inside double quotes.
* When using paired end reads it is **required** that input reads are set with the "{1,2}" pattern. For example: "SRR6307304_{1,2}.fastq". This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq"
* When running hybrid assemblies, mixing short read types or performing multiple assemblies it is advised to use the YAML samplesheet for multi-samples workflows.

### Using the configuration file

All parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is used the pipeline is executed as `nextflow run fmalmeida/mpgap -c ./configuration-file`. Your configuration file is what will tell the pipeline which type of data you have, and which processes to execute. Therefore, it needs to be correctly configured.

* To create a configuration file in your working directory:
  
      nextflow run fmalmeida/mpgap --get_config

### Interactive graphical configuration and execution

#### Via NF tower launchpad (good for cloud env execution)

Nextflow has an awesome feature called [NF tower](https://tower.nf). It allows that users quickly customise and set-up the execution and configuration of cloud enviroments to execute any nextflow pipeline from nf-core, github (this one included), bitbucket, etc. By having a compliant JSON schema for pipeline configuration it means that the configuration of parameters in NF tower will be easier because the system will render an input form.

Checkout more about this feature at: https://seqera.io/blog/orgs-and-launchpad/

<p align="center">
<img src="https://j.gifs.com/GRnqm7.gif" width="500px"/>
</p>

#### Via nf-core launch (good for local execution)

Users can trigger a graphical and interactive pipeline configuration and execution by using [nf-core launch](https://nf-co.re/launch) utility. nf-core launch will start an interactive form in your web browser or command line so you can configure the pipeline step by step and start the execution of the pipeline in the end.

```bash
# Install nf-core
pip install nf-core

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
  + As a workaround, Ryan says to use the `--no_correct` parameter which solves the issue and does not have a negative impact on assembly quality.
  + Therefore, if you run into this error when using unpaired data you can activate this workaroud with `--unicycler_additional_parameters "--no_correct"`.
2. Sometimes, shovill assembler can fail and cause the pipeline to fail due to problems in estimating the genome size. This, is actually super simple to solve! Instead of letting the shovill assembler estimate the genome size, you can pass the information to it and prevent its fail:
    + `--shovill_additional_parameters '--gsize 3m'`

## Citation

To cite this pipeline users can use our Zenodo tag or directly via the github url. Users are encouraged to cite the programs used in this pipeline whenever they are used.

Please, do not forget to cite the software that were used whenever you use its outputs. See [the list of tools](markdown/list_of_tools.md).
