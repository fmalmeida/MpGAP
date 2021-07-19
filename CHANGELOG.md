# MpGAP pipeline changelog

The tracking for changes started in v2.

## v2.3

1. Since pacbio GenomicConsensus and Arrow have reached its end of life, these software have been replaced by their new polisher [gcpp](https://github.com/PacificBiosciences/gcpp).

2. I have identified a small error when using more than one pacbio BAM for polishing. The pipeline, was not loading them all into a nextflow channel list and it was causing the input channel tuple in the pacbio polisher module to have wrong size and wrongly match its variables.
    + The parameter to load pacbio bams has been changed from `--pacbio_all_bam_path` to `--pacbio_bams`.

3. The dockerfile have been changed so it now properly downloads the latest version of software when built. Some of the software wer being download with wget but were being set to a specific freezed version.

4. The nextflow config file have been altered so it now properly calls for a versioned Docker image. Docker images will now have version tags. E.g. `fmalmeida/mpgap:v2.3`.

## v2.2

The main changes in the pipeline are summarized below. The majority of the changes were made to make the pipeline easier to rookies.

### Implementation from issue #7

As stated in issue \#7, we have implemented an option (whenever available) for long reads assemblers to treat the input long reads as **corrected long reads**.

When using the `--corrected_lreads` parameter, the following assemblers will be affected:

* [Canu](https://github.com/marbl/canu)
    + Will trigger the `-corrected` parameter
* [Flye](https://github.com/fenderglass/Flye)
    + Will pass the long reads as `--pacbio-corr` or `--nano-corr`
* [Raven](https://github.com/lbcb-sci/raven)
    + Will trigger the `--weaken` parameter

> Be cautious when using this parameter. If your reads are not corrected, and you use this parameter, you will probably do not generate any contig.

> Additionally, assembly names have been standardised to `{assembler}_assembly`.

### Output directories organization

The output files and folders that have been changed to have a easier and cleaner organization, to provide more readable output files. Now, each assembly is written as a sub-folder under the main output directory (`--outdir`). The sub-folders are written using the input reads basenames:

* Short reads basenames for short reads only assemblies
* Long reads basenames for hybrid and long reads only assemblies

However, whenever running the pipeline for multiple samples at once using glob patterns such as '*' and '?', users are advised to do not perform hybrid assemblies, nor combining both paired and unpaired short reads in short reads only assemblies. Because the pipeline is not yet trained to properly search for the correct pairs, and since nextflow channels are random, we cannot ensure that the combination of data used in these to assembly types will be right. The pipeline treats each input file as a unique sample, and it will execute it individually.

* To date, the use of glob patterns only works properly with long reads only assembly, or short reads only assemblies using either paired or unpaired reads, not both at the same time. For example:
  + `nextflow run [...] --longreads 'my_data/*.fastq' --lr_type 'nanopore' --outdir my_results`
  + The pipeline will load and assembly each fastq in the `my_data` folder and assemble it, writing the results for each read in a sub-folder with the reads basename in the `my_results` output folder.
  + `nextflow run [...] --shortreads_single 'my_data/*.fastq' --outdir my_results`
  + The pipeline will load and assembly each fastq in the `my_data` folder and assemble it, writing the results for each read in a sub-folder with the reads basename in the `my_results` output folder.

> However, we are currently working in a proper way to execute the hybrid and combination of short reads in assemblies for multiple samples at once so that users can properly execute it without confusion. But it will come in v2.3.

### Parameters changed

1. We have removed the `--assembly_type` parameter. It is not necessary anymore. Instead, the pipeline will check the input parameters given, and based on the combination given it will choose between the assembly modes.
2. The `--illumina_polish_longreads_contigs` parameter have been changed to `--strategy_2`.
  + It still does the same thing, it has only been renamed.
3. Turning on/off assemblies. The `--try_*` parameters have been removed. The pipeline now is capable of understanding the assembly mode that will be executed, and to properly select all the assemblers that are capable of performing the assembly mode required.
  + Now, instead of selecting the assemblers to run, by default, the pipeline will run all the assemblers that match the assembly mode. So now, users can turn off some specific assemblers by using the `--skip_*` parameters.

### New assemblers

Three more assemblers have been added to the pipeline:

1. [Shovill](https://github.com/tseemann/shovill)
  + For paired short reads assemblies
2. [Haslr](https://github.com/vpc-ccg/haslr)
  + For hybrid assemblies
3. [Raven](https://github.com/lbcb-sci/raven)
  + For long reads assemblies

### New QC summaries

Finally, [MultiQC](https://multiqc.info/) has been added to the pipeline to provide a rapid and nice comparison between the generated assemblies.

## v2.1

This version have no additions to the pipeline workflow. It has additions in the modes of configuring and executing the pipeline, which are highlighted below.

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/mpgap
```

Checkout the paremeters `--use_tower` and `--tower_token` to activate pipeline execution in nextflow tower.
