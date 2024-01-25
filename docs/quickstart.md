# Quickstart

As an use case, we will use 30X of one of the *Escherichia coli* sequencing data (Biosample: [SAMN10819847](https://www.ncbi.nlm.nih.gov/biosample/10819847))
that is available from a recent study that compared the use of different long read technologies in hybrid assembly of 137 bacterial genomes [[1](https://doi.org/10.1099/mgen.0.000294)].

## Get the data

We have made this subsampled dataset available in [Figshare](https://figshare.com/articles/dataset/Illumina_pacbio_and_ont_sequencing_reads/14036585).


```bash
# Download data from figshare
wget -O reads.zip https://ndownloader.figshare.com/articles/14036585/versions/4

# Unzip
unzip reads.zip
```

Now we have the necessary data to perform the quickstart.

!!! note "Where my outputs go?"

    The pipeline will always create different subdirectories for the resulting assemblies based on your sample IDs, and selected assemblies and assembly strategies selected. All inside the selected `params.output` directory.

## Preparing the input samplesheet

The pipeline reads the input files from a samplesheet in YAML format. A list of available YAML keys to be used in the samplesheet and how to properly create it is available in the [samplesheet reference page](samplesheet.md#).

Here, taking advantage of the `hybrid_strategy` YAML key, we will create a samplesheet entry for the input reads that performs a hybrid assembly in both strategies 1 and 2.

!!! note "The assembly strategies"

  If this key is not used, the pipeline will run the default strategy (1), which can be changed with the parameter `--hybrid_strategy`. For more information on the hybrid assembly strategies please see the [manual reference page](manual.md#).

A proper samplesheet for this data will look like this:

```yaml
# this is a YAML file
# samplesheet file of e. coli 30X reads
# input entry will perform both hybrid strategies
samplesheet:
  - id: ecoli_30X
    illumina:
      - SRR8482585_30X_1.fastq.gz
      - SRR8482585_30X_2.fastq.gz
    nanopore: SRX5299443_30X.fastq.gz
    hybrid_strategy: both
    genome_size: 4m
```

Copy it's content and save it in a file called `samplesheet.yml`, and now, we are able to run the pipeline (check it below).

```bash
# Run the pipeline
nextflow run fmalmeida/mpgap \
  --output _ASSEMBLY \
  --max_cpus 5 \
  --skip_spades \
  --input "samplesheet.yml" \
  --unicycler_additional_parameters ' --mode conservative ' \
  -profile <docker/singularity/conda>
```

!!! tip

    Additional parameters to assemblers can be given with `--{assembler}_additional_parameters`.
    Moreover, specific software can be turned off with the parameters `--skip_{assembler}`.

## About hybrid strategy 2 and long reads polishing

Additionally, for hybrid strategy 2, users can also execute a long reads polishing step in their assemblies prior to the polishing with short reads.

The long reads polishers available are:

* [Medaka](https://github.com/nanoporetech/medaka) and [Nanopolish](https://github.com/jts/nanopolish) for nanopore data;
* [gcpp](https://github.com/PacificBiosciences/gcpp) for pacbio data.

To use them, users must either select a medaka model or pass to the pipeline  the ONT fast5 directory or the pacbio bam file. This will make de pipeline work in the following order: 

1. long reads assembly
2. polishing with long reads models
3. polishing with short reads with Pilon

Please see the [samplesheet](samplesheet.md#) and [manual](manual.md#) reference pages for more information.

## Using test profile


Users can also used a pre-configured test profile which will automatically load a list of SRA run ids for download.

```bash
# short-reads
nextflow run fmalmeida/mpgap -profile test,sreads,<docker/singularity>

# long-reads
nextflow run fmalmeida/mpgap -profile test,lreads,<ont/pacbio>,<docker/singularity>

# hybrid
nextflow run fmalmeida/mpgap -profile test,hybrid,<ont/pacbio>,<docker/singularity>
```

## Afterwards

Now you can used these datasets to, for example, annotate a genome. For this, check out the [Bacannot](https://bacannot.readthedocs.io/en/latest/index.html) pipeline that we've developed for such task.
