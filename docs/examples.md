# CLI usage Examples

!!! warning "Watch your REGEX"

    The pipeline does not concatenate the reads. Whenever you use a pattern such as \* with unpaired reads the pipeline will process each read separately.

## Illumina paired end reads.

This command will select all the read pairs that match the pattern "path-to/SRR*_{1,2}.fastq.gz" and process each pair separately.

```bash
nextflow run fmalmeida/ngs-preprocess \
  --max_cpus 3 \
  --output illumina_paired \
  --shortreads "path-to/SRR*_{1,2}.fastq.gz" \
  --shortreads_type "paired" \
  --fastp_merge_pairs
```

!!! note

    Since `--shortreads` will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.
    
    When using paired end reads it is required that inputs are set with the "{1,2}" pattern. For example: "SRR6307304_{1,2}.fastq". This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq"

    `--fastp_merge_pairs` triggers the Fastp module to merge read pairs.

## Illumina single end reads

This command will select all the reads that match the pattern "path-to/SRR*.fastq.gz" and process each one separately.

```bash
nextflow run fmalmeida/ngs-preprocess \
  --max_cpus 3 \
  --output illumina_single \
  --shortreads "path-to/SRR*.fastq.gz" \
  --shortreads_type "single" \
  --fastp_additional_parameters " --trim_front1 5 --trim_tail1 5 "
```

!!! note

    In this example, we pass on an additional parameter (`--trim_front1 5 --trim_tail1 5`) to Fastp so it trims the reads using a fixed number of bases from the head and tail of reads.
    
    If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

## ONT reads (fastq)

This command will select all the reads that match the pattern "path-to/SRR*.fastq.gz" and process each one separately.

```bash
nextflow run fmalmeida/ngs-preprocess \
  --max_cpus 3 \
  --output ONT \
  --nanopore_fastq "path-to/SRR*.fastq.gz" \
  --lreads_min_length 1000
```

!!! note

  The parameter `--lreads_min_length` applies a minimum read length threshold to filter the reads.

## Pacbio raw (subreads.bam) reads

This command will select all the reads that match the pattern "path-to/m140905_*.subreads.bam" and process each one separately.

```bash
nextflow run fmalmeida/ngs-preprocess \
  --max_cpus 3 \
  --output pacbio_subreads \
  --pacbio_bam "path-to/m140905_*.subreads.bam" \
  --pacbio_get_hifi \
  -with-report
```

!!! note

    The parameter `--pacbio_get_hifi` will make the pipeline **try** to produce the high fidelity pacbio ccs reads.

    `-with-report` will generate nextflow execution reports.
  
    If multiple reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

## Pacbio raw (legacy .bas.h5 to subreads.bam) reads

```bash
nextflow run fmalmeida/ngs-preprocess \
  --pacbio_h5 E01_1/Analysis_Results/ \
  --output E01_1/Analysis_Results/preprocessed \
  --max_cpus 3
```

!!! note

    This example refers to the SMRT Cell data files available at: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly. The path `E01_1/Analysis_Results/` is the directory where the legacy \*.bas.h5 and \*.bax.h5 files are located. The pipeline will load the bas files available in the directory.
    
    Pacbio bas.h5 file and its related bax.h5 files MUST be in the same directory

## Running with a nf-core interactive graphical interface

```bash
nf-core launch fmalmeida/ngs-preprocess
```


## Running with a configuration file

```bash
nextflow run fmalmeida/ngs-preprocess -c nextflow.config
```
