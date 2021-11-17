# MpGAP pipeline changelog

The tracking for changes started in v2.

## v3.1

This release addresses the issue https://github.com/fmalmeida/MpGAP/issues/17.

Now the pipeline mainly relies on packages built through conda into a main environment, making it possible to use it with docker, conda or singularity by properly using the `-profile` parameter.

Please read more about it at: https://github.com/fmalmeida/MpGAP#selecting-between-profiles

## v3.0.1

### hotfix

Super small fix to properly load YAML file when using the pipeline with cloud computing environments such as AWS/S3-bucket:

```bash
# from
parameter_yaml = new FileInputStream(new File(params.input))
# to
parameter_yaml = file(params.input).readLines().join("\n")
```

## v3.0

### hotfix

* Since shasta release v0.8 (Oct/2021), shasta now expects users to select a pre-set configuration file.
    + This has been included as a parameter `--shasta_config` that sets a default to `Nanopore-Oct2021`
      + Please read the [shasta configuration manual](https://chanzuckerberg.github.io/shasta/Configurations.html) page to know the available models.
      + It can also be passed from inside the YAML file in a sample-specific manner.
      + Please read more about it in the online documentation: [Samplesheet configuration](https://mpgap.readthedocs.io/en/latest/samplesheet.html) and [Parameters manual](https://mpgap.readthedocs.io/en/latest/manual.html)

### input configuration

* A YAML samplesheet file has been implemented in order to properly use the pipeline, in all workflow types, for multiple samples at once.
    + Because of that, we had to remove the possibility to pass the input reads via the command line and now, the files input data files, must always be set inside the YAML samplesheet, even if analysing only one sample.
    + Read more at: https://mpgap.readthedocs.io/en/latest/samplesheet.html
    + Check the template samplesheet at: https://github.com/fmalmeida/mpgap/blob/master/example_samplesheet.yaml
    + The samplesheet is given with the parameter `--input`
* Due to the implementation above, the folowing parameters are now deprecated, since they are now set inside the YAML file:
    + `--longreads`
    + `--lr_type`
    + `--pacbio_bam`
    + `--nanopolish_fast5Path`
    + `--shortreads_paired`
    + `--shortreads_single`
* A major change has also ocurred with the `wtdbg2_techonology` parameter
    + Now, by default, the pipeline will check wheter long reads inputs (for each sample) are nanopore or pacbio
    + If they are nanopore, the wtdbg2 techonology parameter is automatically set to `ont`
    + If they are pacbio, the wtdbg2 techonology parameter is automaically set to `sq`
    + This wtdbg2 parameter has the following options: "ont" for Nanopore reads, "rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS reads.
    + It can also be passed from inside the YAML file in a sample-specific manner.
    + Please read more about it in the online documentation: [Samplesheet configuration](https://mpgap.readthedocs.io/en/latest/samplesheet.html) and [Parameters manual](https://mpgap.readthedocs.io/en/latest/manual.html)

### nomenclature change

* In order to make it simple and natural, two changes ocurred in input/output parameters
    + The `--outdir` parameter is now `--output`
    + The `--medaka_sequencing_model` parameter is now `--medaka_model`
    + The `--corrected_lreads` parameter is now `--corrected_long_reads`.
      + It can also be passed from inside the YAML file in a sample-specific manner.
      + Please read more about it in the online documentation: [Samplesheet configuration](https://mpgap.readthedocs.io/en/latest/samplesheet.html) and [Parameters manual](https://mpgap.readthedocs.io/en/latest/manual.html)
    + The `--genomeSize` parameter is now `--genome_size`.
      + It can also be passed from inside the YAML file in a sample-specific manner.
      + Please read more about it in the online documentation: [Samplesheet configuration](https://mpgap.readthedocs.io/en/latest/samplesheet.html) and [Parameters manual](https://mpgap.readthedocs.io/en/latest/manual.html)
    + The `--strategy_2` parameter is now `--hybrid_strategy` which expects a value indicating the strategies to perform.
      + It still defaults to strategy 1
      + It can also be passed from inside the YAML file in a sample-specific manner.
      + Please read more about it in the online documentation: [Samplesheet configuration](https://mpgap.readthedocs.io/en/latest/samplesheet.html) and [Parameters manual](https://mpgap.readthedocs.io/en/latest/manual.html)
  
### comments

* Since this changes are major changes, the pipeline main version has changed and it is now in v3.0
    + The docker image is `fmalmeida/mpgap:v3.0`.

## v2.3.1

This patch release is related to the issue [#19](https://github.com/fmalmeida/MpGAP/issues/19) which raises attention that shovill was not being used to its fully extent. Shovill was just being used with spades as its core.

Now this has been fixed and shovill now will produce three assemblies for comparisons, with both spades, megahit and skesa as its core. However this limits the functionality of the  `--shovill_additional_parameters`. Since these three assemblers are already done with shovill, users cannot pass anymore the shovill `--assembler` parameter with `--shovill_additional_parameters`. Which means, something as `--shovill_additional_parameters " --assembler skesa "` will raise an error.

Also, the CLI help messages and in the config were made a bit more clear to avoid confusions.

## v2.3

1. Since pacbio GenomicConsensus and Arrow have reached its end of life, these software have been replaced by their new polisher [gcpp](https://github.com/PacificBiosciences/gcpp).

2. I have identified a small error when using more than one pacbio BAM for polishing. The pipeline, was not loading them all into a nextflow channel list and it was causing the input channel tuple in the pacbio polisher module to have wrong size and wrongly match its variables.
    + The parameter to load pacbio bams has been changed from `--pacbio_all_bam_path` to `--pacbio_bams`.

3. The dockerfile have been changed so it now properly downloads the latest version of software when built. Some of the software were being download with wget but were being set to a specific freezed version.

4. The nextflow config file have been altered so it now properly calls for a versioned Docker image. Docker images will now have version tags. E.g. `fmalmeida/mpgap:v2.3`.

5. Two other longreads assemblers have been added: wtdbg2 and shasta. Please read the docs to check the parameters related to them.

    + One of this assemblers, wtdbg2, has an additional parameter (`--wtdbg2_technology`) that must be properly set when assembling pacbio reads with it.

6. A new parameter called `--prefix` has been added so users can create custom prefixes for their samples. By default (if not used), the pipeline will create a custom prefix using the input reads names.

7. A new document called ASSEMBLY_SUMMARY.txt is now given. Before the pipeline only produced a HTML report with multiqc, however, HTML can not be read in the terminal. Thus, we now also produce this text file.

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
