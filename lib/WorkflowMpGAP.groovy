//
// This file holds several functions specific to the the fmalmeida/mpgap pipeline
//

class WorkflowMpGAP {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        params.hybrid_strategy = params.hybrid_strategy.toString()
        if (!params.get_config && !params.get_samplesheet && !params.help) {
            if (!params.input) {
            log.error "ERROR!\nA major error has occurred!\n\t==> A samplesheet has not been provided. Please, provide a samplesheet to run the analysis.\n\t Online documentation is available at: https://mpgap.readthedocs.io/en/latest/\nPlease, read the docs.\nCheers."
            System.exit(1)
            }
        }
    }

}
