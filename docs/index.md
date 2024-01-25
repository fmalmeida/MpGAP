# Welcome to <u>ngs-preprocess</u> pipeline documentation

<img src="./lab_logo.png" width="300px">

[![F1000 Paper](https://img.shields.io/badge/Citation%20F1000-10.12688/f1000research.139488.1-orange)](https://doi.org/10.12688/f1000research.139488.1)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/fmalmeida/ngs-preprocess?include_prereleases&label=Latest%20release)](https://github.com/fmalmeida/ngs-preprocess/releases)
[![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://ngs-preprocess.readthedocs.io/en/latest/?badge=latest)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/ngs-preprocess/blob/master/LICENSE)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40fmarquesalmeida-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/fmarquesalmeida)
[![Zenodo Archive](https://img.shields.io/badge/Zenodo-Archive-blue)](https://doi.org/10.5281/zenodo.3451405)

## About

[ngs-preprocess](https://github.com/fmalmeida/ngs-preprocess) is a pipeline designed to provide an easy-to-use framework for preprocessing sequencing reads from Illumina, Pacbio and Oxford Nanopore platforms. It is developed with [Nextflow](https://www.nextflow.io/docs/latest/index.html) and [Docker](https://www.docker.com/).

## Workflow

The pipeline wraps up the following tools and analyses:

| Software | Analysis |
| :------- | :------- |
| [sra-tools](https://github.com/ncbi/sra-tools) & [entrez-direct](https://anaconda.org/bioconda/entrez-direct) | Interaction with SRA database for fetching fastqs and metadata |
| [fastp](https://github.com/OpenGene/fastp) | Fast all-in-one preprocessing for FastQ files |
| [porechop](https://github.com/rrwick/Porechop)** | ONT reads trimming and demultiplexing |
| [pycoQC](https://github.com/tleonardi/pycoQC) | ONT reads QC |
| [NanoPack](https://github.com/wdecoster/nanopack) | Long reads QC and filter |
| [bax2bam](https://anaconda.org/bioconda/bax2bam) | Convert PacBio bax files to bam |
| [bam2fastx](https://github.com/PacificBiosciences/pbtk#bam2fastx) | Extract reads from PacBio bam files |
| [lima](https://github.com/PacificBiosciences/barcoding) | PacBio reads demultiplexing |
| [pacbio ccs](https://ccs.how/) | Generate PacBio Highly Accurate Single-Molecule Consensus Reads |

!!! info "About porechop"

    Although discontinued since 2018, porechop is included as a legacy compatibility for old nanopore runs, old sequencing kit libraries and old sequencer versions.
    
    However, the newest versions of MinKNOW is able to output trimmed and demultiplexed fastq data, meaning this step is not required anymore.

    Finally, it is also okay to not remove adapters from reads as some assemblers may be aware and even benefit of the sequences.

!!! note "Quickstart"

    A [quickstart](quickstart.md#) is available so you can quickly get the gist of the pipeline's capabilities.

## Usage

The pipeline's common usage is very simple as shown below:

```bash
# usual command-line
nextflow run fmalmeida/ngs-preprocess \
    --sra_ids "list_of_sra.txt" \
    --lreads_min_length 750 \
    --output "./preprocessed_data" \
    ...
```

!!! quote

    Some parameters are required, some are not. Please read the pipeline's manual reference to understand each parameter.

## Citation

In order to cite this pipeline, please refer to:

> Almeida FMd, Campos TAd and Pappas Jr GJ. Scalable and versatile container-based pipelines for de novo genome assembly and bacterial annotation. F1000Research 2023, 12:1205 (<https://doi.org/10.12688/f1000research.139488.1>)

## Support contact

Whenever a doubt arise feel free to contact me via the [github issues](https://github.com/fmalmeida/ngs-preprocess/issues).
