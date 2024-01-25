# Welcome to <u>MpGAP</u> pipeline documentation

<img src="./lab_logo.png" width="300px">

[![F1000 Paper](https://img.shields.io/badge/Citation%20F1000-10.12688/f1000research.139488.1-orange)](https://doi.org/10.12688/f1000research.139488.1)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/fmalmeida/mpgap?include_prereleases&label=Latest%20release)](https://github.com/fmalmeida/mpgap/releases)
[![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://mpgap.readthedocs.io/en/latest/?badge=latest)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/mpgap/blob/master/LICENSE)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40fmarquesalmeida-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/fmarquesalmeida)
[![Zenodo Archive](https://img.shields.io/badge/Zenodo-Archive-blue)](https://doi.org/10.5281/zenodo.3445485)

## About

[MpGAP](https://github.com/fmalmeida/mpgap) is a pipeline developed with [Nextflow](https://www.nextflow.io/docs/latest/index.html) and [Docker](https://www.docker.com/). It was designed to provide an easy-to-use framework for *de novo* genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes.

## Workflow

The pipeline wraps up the following tools and analyses:

| Software | Analysis |
| :------- | :------- |
|  [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler), [Raven](https://github.com/lbcb-sci/raven), [Shasta](https://github.com/chanzuckerberg/shasta) and [wtdbg2](https://github.com/ruanjue/wtdbg2) | Long reads assembly |
| [Haslr](https://github.com/vpc-ccg/haslr), [Unicycler](https://github.com/rrwick/Unicycler) and [SPAdes](https://github.com/ablab/spades) | Hybrid assembly |
| [Shovill](https://github.com/tseemann/shovill), [Unicycler](https://github.com/rrwick/Unicycler), [Megahit](https://github.com/voutcn/megahit) and [SPAdes](https://github.com/ablab/spades) | Short reads assembly |
| [Nanopolish](https://github.com/jts/nanopolish), [Medaka](https://github.com/nanoporetech/medaka), [gcpp](https://github.com/PacificBiosciences/gcpp) and [Pilon](https://github.com/broadinstitute/pilon) | Assembly polishing |
| [Quast](https://github.com/ablab/quast) and [MultiQC](https://multiqc.info/) | Assembly QC |

!!! note "Quickstart"

    A [quickstart](quickstart.md#) is available so you can quickly get the gist of the pipeline's capabilities.

## Usage

The pipeline's common usage is very simple as shown below:

```bash
# usual command-line
nextflow run fmalmeida/mpgap \
    -profile docker \
    --output ./results \
    --tracedir ./results/pipeline_info \
    --input input.yml \
    --max_cpus 20 \
    --max_memory '40.GB' \
    ...
```

!!! quote

    Some parameters are required, some are not. Please read the pipeline's manual reference to understand each parameter.

## Citation

In order to cite this pipeline, please refer to:

> Almeida FMd, Campos TAd and Pappas Jr GJ. Scalable and versatile container-based pipelines for de novo genome assembly and bacterial annotation. F1000Research 2023, 12:1205 (<https://doi.org/10.12688/f1000research.139488.1>)

## Support contact

Whenever a doubt arise feel free to contact me via the [github issues](https://github.com/fmalmeida/mpgap/issues).
