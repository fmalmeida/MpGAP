//
// This file holds several functions specific to the main.nf workflow in the fmalmeida/ngs-preprocess pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.5281/zenodo.3445485\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/fmalmeida/ngs-preprocess#citation"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --max_cpus 3 --output illumina_paired --shortreads \"path-to/SRR*_{1,2}.fastq.gz\" --shortreads_type \"paired\" --fastp_merge_pairs" + "\n  " + "nextflow run ${workflow.manifest.name} --max_cpus 3 --output ONT --nanopore_fastq \"path-to/SRR*.fastq.gz\" --lreads_min_length 1000" + "\n" + "\n  " + "More CLI examples at: https://ngs-preprocess.readthedocs.io/en/latest/examples.html#"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, false)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(false)
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, false)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(false)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Download template config
        if (params.get_config) {
            new File("MPGAP.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/conf/defaults.config").getText())
            log.info """

            MPGAP.config file saved in working directory
            After configuration, run:
              nextflow run fmalmeida/mpgap -c ./MPGAP.config
            Nice code

            """.stripIndent()
            System.exit(0)
        }

        // Download template samplesheet
        if (params.get_samplesheet) {
            new File("MPGAP.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/example_samplesheet.yml").getText())
            log.info """

            Samplesheet (MPGAP_samplesheet.yml) file saved in working directory
            Nice code!

            """.stripIndent()
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        // if (params.enable_conda) {
        //     Utils.checkCondaChannels(log)
        // }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

    }

}
