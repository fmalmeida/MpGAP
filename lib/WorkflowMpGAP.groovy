//
// This file holds several functions specific to the the fmalmeida/mpgap pipeline
//

class WorkflowMpGAP {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.get_config && !params.get_samplesheet && !params.help) {
            if (!params.input) {
            log.error "ERROR!\nA major error has occurred!\n\t==> A samplesheet has not been provided. Please, provide a samplesheet to run the analysis.\n\t Online documentation is available at: https://mpgap.readthedocs.io/en/latest/\nPlease, read the docs.\nCheers."
            System.exit(1)
            }
        }

        if (params.hybrid_strategy.toString() != "1" && params.hybrid_strategy.toString() != "2" && params.hybrid_strategy.toString() != "both") {
            log.error "ERROR!\nA major error has occurred!\n\t==>  Parameter --hybrid_strategy must be either 1, 2 or both.\n\t Online documentation is available at: https://mpgap.readthedocs.io/en/latest/\nPlease, read the docs.\nCheers."
            System.exit(1)
        }

        if (params.corrected_longreads && params.high_quality_longreads) {
            log.error "ERROR!\nA major error has occurred!\n\t==>  Parameters --corrected_longreads and --high_quality_longreads were used at the same time. These activate assembler configurations for reads of different quality levels. Cannot be used at the same time ( uncorrected < corrected < high_quality ).\nCheers."
            System.exit(1)
        }
    }

}
