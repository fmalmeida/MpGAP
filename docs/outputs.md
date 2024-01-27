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
├── input.yml                       # Copy of given input samplesheet for data provenance
└── pipeline_info                   # directory containing the nextflow execution reports
    ├── mpgap_report_2023-12-28_12-25-18.html
    ├── mpgap_timeline_2023-12-28_12-25-18.html
    └── mpgap_tracing_2023-12-28_12-25-18.txt
```

## Example of QC outputs

Here I am going to display just a very few examples of results produced, focusing on the QC, as the main result is a normal assembly, performed by each assembler. Just for reference, I am also displaying the drawing of the assembly graph with Bandage, for the results of the Flye assembler

**Bandage visualization of Flye assembly graph**

<center>
  <img src="../assets/LengthvsQualityScatterPlot_dot.png" width="85%">
</center>

**MultiQC Report - HTML**

Open it [here](../assets/NanoPlot-report.html).

**Quast Report of Flye assembly - HTML**

Open it [here](../assets/NanoPlot-report.html).
