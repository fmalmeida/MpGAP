# Output files

Here, using the results produced in the [Non-bacterial dataset section](non_bacteria.md#), we give users a glimpse over the main outputs produced by bacannot. The command used in the quickstart wrote the results under the `preprocessed_reads` directory.

!!! note

    Please take note that the pipeline uses the directory set with the `--output` parameter as a storage place in which it will create a folder for the final pre-processed reads and for the intermediate files, separated by sequencing technology.

## Directory tree

After a successful execution, you will have something like this:

```bash

# Directory tree from the running dir
preprocessed_reads
# directory containing the final results of the data cleaning
├── final_output                   
│   └── nanopore
│       └── SRR23337893.filtered.fq.gz
# directory containing the nextflow execution reports
├── pipeline_info
│   ├── ngs_preprocess_report_2023-11-18_10-07-36.html
│   ├── ngs_preprocess_timeline_2023-11-18_10-07-36.html
│   ├── ngs_preprocess_tracing_2023-11-18_10-07-36.txt
# directory containing the intermediate files produced by the tools used during pre-processing, and, QC
├── preprocessing_outputs
│   └── nanopore
│       ├── porechop
│       └── QC
# directory containing the intermediate files when downloading data from SRA
└── SRA_FETCH
    ├── FASTQ
    │   └── SRR23337893_data
    └── SRR23337893_sra_runInfo.csv
```

## Example of QC outputs

Here I am going to display just a very few examples of results produced, focusing on the QC, as the main result is a cleaned FASTQ file.

**Length versus Quality Scatterplot**

<center>
  <img src="../assets/LengthvsQualityScatterPlot_dot.png" width="85%">
</center>

**NanoPlot Report HTML**

Open it [here](../assets/NanoPlot-report.html).

**NanoStats Report TXT**

Open it [here](../assets/NanoStats.txt).
