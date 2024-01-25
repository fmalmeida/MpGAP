# Manual

```bash
# Get help in the command line
nextflow run fmalmeida/ngs-preprocess --help
```

!!! tip

    All these parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution cleaner and more readable. See a [config example](config.md#).

## Input description

* path to fastq files containing sequencing reads
* path to Pacbio .bam or .h5 files containing raw data
* path containing list of SRA IDs

!!! warning "Watch your input"

    Users must **never** use hard or symbolic links. This will make nextflow fail.

    Whenever using REGEX for a pattern match, for example "illumina/SRR9847694_{1,2}.fastq.gz" or "illumina/SRR*.fastq.gz", it MUST ALWAYS be inside double quotes.

    **Remember:** the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will process each read (or pair) that match this pattern separately.

## Output options

| <div style="width:180px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--output`  | :material-check: | NA       | Directory to store output files |

## Max job request

| <div style="width:120px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--max_cpus`  | :material-close: | 4  | Max number of threads a job can use across attempts |
| `--max_memory` | :material-close: | 6.GB | Max amount of memory a job can use across attempts |
| `--max_time` | :material-close: | 40.h | Max amount of time a job can take to run

## SRA IDs as input

As of version v2.5, users can also select data directly from SRA. One just need to provide a txt file containing SRA run ids, one per line, e.g. [Example](https://github.com/fmalmeida/test_datasets/blob/main/sra_ids.txt).

| <div style="width:160px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--sra_ids`      | :material-close: | NA | Path to txt file containing list of SRA run IDs |

## Short reads input

| <div style="width:220px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--shortreads` | :material-check: | NA | String Pattern to find short reads. Example: "SRR6307304_{1,2}.fastq" |
| `--shortreads_type` | :material-check: | NA | (single \| paired). Tells whether input is unpaired or paired end |
| `--fastp_average_quality` | :material-close: | 20 | Fastp will filter out reads with mean quality less than this |
| `--fastp_correct_pairs` | :material-close: | false | If set, tells Fastp to try to correct paired end reads. Only works for paired end reads |
| `--fastp_merge_pairs` | :material-close: | false | If set, tells Fastp to try to merge read pairs |
| `--fastp_additional_parameters` | :material-close: | false | Pass on any additional parameter to Fastp. The tool's parameters are described in their [manual](https://github.com/OpenGene/fastp) |

## Long reads input

| <div style="width:220px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--lreads_min_length` | :material-close: | 500 | Length min. threshold for filtering long reads (ONT or Pacbio) |
| `--lreads_min_quality` | :material-close: | 5 | Quality min. threshold for filtering long reads (ONT or Pacbio) |
| `--nanopore_fastq` | :material-check: | NA | Sets path to nanopore fastq files. Pre-processes basecalled long reads |
| `--nanopore_is_barcoded` | :material-close: | false | Tells whether your data (Nanopore or Pacbio) is barcoded or not. It will split barcodes into single files. Users with legacy pacbio data need to first produce a new barcoded_subreads.bam file |
| `--nanopore_sequencing_summary` | :material-close: | NA | Path to nanopore 'sequencing_summary.txt'. Using this will make the pipeline render a sequencing statistics report using pycoQC. pycoQC reports will be saved using the files basename, so please, use meaningful basename, such as: sample1.txt, sample2.txt, etc. Preferentially, using the same basename as the fastq |
| `--pacbio_bam` | :material-close: | NA | Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ. Always keep subreads.bam and its relative subreads.bam.pbi files in the same directory |
| `--pacbio_h5` | :material-close: | NA | Path to directory containing legacy bas.h5 data file (1 per directory). It will be used to extract reads in FASTQ file. All its related files (e.g. bax.h5 files) must be in the same directory |
| `--pacbio_barcodes` | :material-close: | NA | Path to xml/fasta file containing barcode information. It will split barcodes into single files. Will be used for all pacbio inputs, h5 or bam |
| `--pacbio_barcode_design` | :material-close: | same | Select the combination of barcodes for demultiplexing. Options: same, different, any |
| `--pacbio_get_hifi` | :material-close: | false | Whether or not to try to compute CCS reads. Will be used for all pacbio inputs, h5 or bam |

All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a [config](config.md#) example.


## Examples

For a better understanding of the usage we provided a feel examples. See some [examples](examples.md#).
