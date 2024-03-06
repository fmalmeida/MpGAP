# Output files

Here, using the results produced in the [Non-bacterial dataset section](non_bacteria.md#), we give users a glimpse over the main outputs produced by MpGAP. The command used in there wrote the results under the `genome_assembly` directory.

!!! note

    Please take note that the pipeline uses the directory set with the `--output` parameter as a storage place in which it will create a folder for the final results, separated by sample, technology and assembly strategy.

## Directory tree

After a successful execution, you will have something like this:

```bash
# Directory tree from the running dir
genome_assembly
├── aspergillus_fumigatus           # directory containing the assembly results for a given sample these are written with the 'id' value. In our example we have only one, but if input data samplesheet had more samples we would have one sub-directory for each.
│   └── longreads_only              # results for long reads only assembly. A sub-directory is created for results of each assembly strategy to allow you running multiple strategies at once
│       ├── 00_quality_assessment   # QC reports
│       ├── canu                    # Canu assembly
│       ├── flye                    # Flye assembly
│       ├── medaka_polished_contigs # Assemblies of all assemblers polished with medaka
│       ├── raven                   # Raven assembly
│       ├── shasta                  # Shasta assembly
│       └── wtdbg2                  # Shasta assembly
├── bacannot_samplesheet.yml        # a template input ready for bacannot pipeline
├── final_assemblies                # A folder contatining a copy of all the assemblies generated, raw and polished
│   ├── aspergillus_fumigatus_canu_assembly.fasta
│   ├── aspergillus_fumigatus_canu_medaka_consensus.fa
│   ├── aspergillus_fumigatus_flye_assembly.fasta
│   ├── aspergillus_fumigatus_flye_medaka_consensus.fa
│   ├── < ... > etc.
├── input.yml                       # Copy of given input samplesheet for data provenance
└── pipeline_info                   # directory containing the nextflow execution reports
    ├── mpgap_report_2023-12-28_12-25-18.html
    ├── mpgap_timeline_2023-12-28_12-25-18.html
    └── mpgap_tracing_2023-12-28_12-25-18.txt
```

## The pre-formatted Bacannot input samplesheet

Once finished, the pipeline also generates a file called `bacannot_samplesheet.yml` (showed below). Basically this samplesheet defines all the **minimum** definitions in order to annotate these generated genomes using the [Bacannot](https://bacannot.readthedocs.io/en/latest/) pipeline.

```yaml
samplesheet:
  - id: aspergillus_fumigatus_shasta
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_shasta_assembly.fasta
  - id: aspergillus_fumigatus_flye
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_flye_assembly.fasta
  - id: aspergillus_fumigatus_raven
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_raven_assembly.fasta
  - id: aspergillus_fumigatus_wtdbg2
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_wtdbg2_assembly.fasta
  - id: aspergillus_fumigatus_canu
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_canu_assembly.fasta
  - id: aspergillus_fumigatus_shasta_medaka
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_shasta_medaka_consensus.fa
  - id: aspergillus_fumigatus_canu_medaka
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_canu_medaka_consensus.fa
  - id: aspergillus_fumigatus_raven_medaka
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_raven_medaka_consensus.fa
  - id: aspergillus_fumigatus_flye_medaka
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_flye_medaka_consensus.fa
  - id: aspergillus_fumigatus_wtdbg2_medaka
    assembly: /array1/falmeida/git_repos/MpGAP/testing/results/final_assemblies/aspergillus_fumigatus_wtdbg2_medaka_consensus.fa
```

!!! note

    One must keep in mind that, this template samplesheet contains only the **bare minimum** to launch bacannot but many other customizations are possible. For example, one can also set, for each input genome, a different resfinder panel if not wanting to run the same for all. And many other things. Therefore, users can/must use this output as a template for easily customization of the bacannot pipeline input to readily use the results of the mpgap pipeline.

    For more information, please refer to the [Bacannot](https://bacannot.readthedocs.io/en/latest/) documentation.

## Example of QC outputs

Here I am going to display just a very few examples of results produced, focusing on the QC, as the main result is a normal assembly, performed by each assembler.

**Summary of Assembly Statistics in TXT format**

Open it [here](../assets/ASSEMBLY_SUMMARY.txt).

**MultiQC Report - HTML**

Open it [here](../assets/multiqc_report.html).

**Quast Report of Flye assembly - HTML**

Open it [here](../assets/flye_medaka/report.html).
